import numpy as np
from skimage import measure, morphology
from scipy import ndimage  # Correct import for distance transforms
import os
import argparse
import tifffile
import pandas as pd

def calculate_ermicc(mito_label, er_mask, pixel_size_nm=1, max_dist_nm=20):
    """
    mito_label: 2D array, unique label per mitochondrion.
    er_mask: 2D binary array (1 for ER, 0 for bg).
    pixel_size_nm: conversion factor from pixels to nanometers.
    max_dist_nm: threshold distance for contact.
    """
    # Convert unit to pixels
    max_dist_px = max_dist_nm / pixel_size_nm
    
    # Distance map of ER: distance from each pixel to the nearest ER pixel
    er_dist_map = ndimage.distance_transform_edt(1 - er_mask)
    
    # Identify contact mask (mitochondria within threshold of ER)
    contact_mask = (er_dist_map <= max_dist_px) & (mito_label > 0)
    
    # Create labeled masks to preserve the mito labels
    labeled_contact_mask = np.where(er_dist_map <= max_dist_px, mito_label, 0)
    global_interface_mask = np.zeros(mito_label.shape, dtype=mito_label.dtype)

    # Get all region properties to compute longest distance (Feret diameter)
    props = measure.regionprops(mito_label)
    # Create a dictionary for quick lookup by label ID
    prop_dict = {p.label: p for p in props}
    feret_diameters = {prop.label: prop.feret_diameter_max for prop in props}
    
    results = []
    # Get all unique mitochondria labels
    mito_ids = np.unique(mito_label[mito_label > 0])
    
    for m_id in mito_ids:
        # Mask for current mitochondrion
        mito_mask = (mito_label == m_id)
        
        # Perimeter of this mitochondrion
        perimeter_px = measure.perimeter(mito_mask)
        perimeter_nm = perimeter_px * pixel_size_nm
        
        # Interface: part of the mitochondria perimeter within distance threshold
        # Find edges of current mito
        mito_boundary = mito_mask ^ morphology.binary_erosion(mito_mask)
        interface_mask = mito_boundary & contact_mask
        
        # Accumulate to global interface mask preserving the mito label
        global_interface_mask[interface_mask] = m_id
        
        interface_len_px = np.sum(interface_mask)
        interface_len_nm = interface_len_px * pixel_size_nm
        
        # Average distance of the contact interface
        if interface_len_px > 0:
            avg_dist_px = np.mean(er_dist_map[interface_mask])
            avg_dist_nm = avg_dist_px * pixel_size_nm
            
            # ERMICC = interface_length / (mito_perimeter * distance)
            if perimeter_nm * avg_dist_nm > 0:
                ermicc = interface_len_nm / (perimeter_nm * avg_dist_nm)
            else:
                ermicc = 0
        else:
            avg_dist_nm = 0
            ermicc = 0

        prop = prop_dict.get(m_id)
        # major_axis_length is in pixels; convert to nm
        major_axis_length_nm = prop.major_axis_length * pixel_size_nm if prop else 0
        
        # Longest distance (Feret diameter) within mitochondrion
        longest_dist_px = feret_diameters.get(m_id, 0.0)
        longest_dist_nm = longest_dist_px * pixel_size_nm
        
            
        results.append({
            "mito_id": m_id,
            "major_axis_length_nm": major_axis_length_nm,
            "longest_distance_nm": longest_dist_nm,
            "interface_length_nm": interface_len_nm,
            "mito_perimeter_nm": perimeter_nm,
            "avg_distance_nm": avg_dist_nm,
            "ermicc": ermicc
        })
        
    return results, labeled_contact_mask, global_interface_mask


def main():
    parser = argparse.ArgumentParser(description="Measure global and per-mito contact.")
    parser.add_argument("--mito_folder", required=True)
    parser.add_argument("--er_folder", required=True)
    parser.add_argument("--pixel_size", type=float, required=True)
    parser.add_argument("--min_contact_distance", type=float, required=True)
    parser.add_argument("--output_folder", default="ERMICC_measure")

    args = parser.parse_args()
    os.makedirs(args.output_folder, exist_ok=True)
    
    summary_results = []
    all_mito_results = []
    
    for filename in os.listdir(args.mito_folder):
        if filename.lower().endswith(('.tif', '.tiff')):
            mito_path = os.path.join(args.mito_folder, filename)
            er_path = os.path.join(args.er_folder, filename)
            
            if os.path.exists(er_path):
                print(f"Processing {filename}...")

                mito_labeled = tifffile.imread(mito_path)
                er_mask = tifffile.imread(er_path) > 0

                # Check if dimension of both images are identical
                if mito_labeled.shape != er_mask.shape:
                    print(f"Skipping {filename}: dimension mismatch")
                    continue

                results, contact_mask, interface_mask = calculate_ermicc(mito_labeled, er_mask, pixel_size_nm=args.pixel_size, max_dist_nm=args.min_contact_distance)
                
                # Add filename to each mitochondrion result and extend the global list
                for r in results:
                    r["filename"] = filename
                all_mito_results.extend(results)

                summary_results.append({
                    "filename": filename,
                    "mito_count": len(results),
                    "global_ermicc": np.mean([r["ermicc"] for r in results]) if len(results) > 0 else 0
                })

                # Save multi-channel QC stack (Mito, ER, contact mask, interface mask)
                # Note: channels 1, 3, and 4 preserve the original mitochondrion label IDs
                qc_stack = np.stack([
                    mito_labeled.astype(np.uint16),
                    er_mask.astype(np.uint16) * 255,
                    contact_mask.astype(np.uint16),
                    interface_mask.astype(np.uint16)
                ], axis=0)

                tifffile.imwrite(
                    os.path.join(args.output_folder, f"QC_{filename}"),
                    qc_stack,
                    imagej=True,
                    metadata={'axes': 'CYX'}
                )

    # Reorder and save summary results
    df_summary = pd.DataFrame(summary_results)
    cols_summary = ['filename'] + [c for c in df_summary.columns if c != 'filename']
    df_summary[cols_summary].to_csv(os.path.join(args.output_folder, "summary_results_ermicc.csv"), index=False)

    # Reorder and save per-mito results
    df_mito = pd.DataFrame(all_mito_results)
    cols_mito = ['filename'] + [c for c in df_mito.columns if c != 'filename']
    df_mito[cols_mito].to_csv(os.path.join(args.output_folder, "per_mito_results_ermicc.csv"), index=False)

    print("Analysis complete. Files saved: summary_results_ermicc.csv and per_mito_results_ermicc.csv")

if __name__ == "__main__":
    main()