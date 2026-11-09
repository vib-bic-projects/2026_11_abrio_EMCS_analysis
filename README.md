# 2026 11 abrio EMCS analysis

## Conda environment

## Run the script
```bash
python measure_ermicc.py --mito_folder <mito_folder> --er_folder <er_folder> --pixel_size 1.4286 --min_contact_distance 20.0 --output_folder <output_folder>
```
### Parameters

| Parameter | Required | Default | Description |
| --------- | -------- | ------- | ----------- |
| `--mito_folder` | Yes | – | Folder containing the mitochondria **label** images (`.tif` / `.tiff`), with one unique integer ID per mitochondrion (0 = background). |
| `--er_folder` | Yes | – | Folder containing the ER masks (`.tif` / `.tiff`). Any pixel value > 0 is treated as ER. Each file must have the **same filename and dimensions** as its counterpart in `--mito_folder`; images without a matching ER file, or with a size mismatch, are skipped. |
| `--pixel_size` | Yes | – | Physical pixel size in nm/pixel (e.g. `1.4286`). Used to convert all lengths and distances from pixels to nanometers. |
| `--min_contact_distance` | Yes | – | Contact distance threshold $d_{\text{max}}$ in nm (e.g. `20.0`). Mitochondrial boundary pixels located at a distance **≤** this value from the ER are counted as contact. Despite its name, it acts as the *maximum* distance for a contact. |
| `--output_folder` | No | `ERMICC_measure` | Folder where results are written (created if it does not exist). |

### Outputs

Written in `--output_folder`:

- `summary_results_ermicc.csv`: one row per image (`filename`, `mito_count`, `global_ermicc` = mean ERMICC of all mitochondria in the image).
- `per_mito_results_ermicc.csv`: one row per mitochondrion (`filename`, `mito_id`, `major_axis_length_nm`, `longest_distance_nm`, `interface_length_nm`, `mito_perimeter_nm`, `avg_distance_nm`, `ermicc`).
- `QC_<filename>`: 4-channel ImageJ TIFF for quality control: (1) mitochondria labels, (2) ER mask, (3) contact mask (mito pixels within the threshold, with their label IDs), (4) contact interface (boundary pixels in contact, with their label IDs).

## Analysis steps:
### Step-by-Step Measurement & Formulas

The ERMICC (ER-Mitochondria Contact Coefficient) is computed through the following sequential steps for each mitochondrion:

#### Step 1: ER Distance Map
A Euclidean Distance Transform (EDT) is calculated on the background of the ER binary mask $ER$ to compute the distance from each pixel to the nearest ER pixel. The physical distance in nanometers is:
$$D_{\text{ER}}(x, y) = \text{EDT}(1 - ER(x, y)) \times \mathrm{pixel\_size}$$

#### Step 2: Mitochondria Perimeter
For each mitochondrion $i$ with mask $M_i$, the boundary pixels are identified using morphological erosion:
$$\text{Boundary}_i(x, y) = M_i(x, y) \oplus \text{Erode}(M_i(x, y))$$

The physical perimeter length in nanometers is computed using the contours of the mask:
$$P_i = \text{perimeter}(M_i) \times \mathrm{pixel\_size}$$

#### Step 3: Contact Interface Length
The contact interface is defined as the subset of the mitochondrion's boundary pixels that reside within the contact distance threshold ($d_{\text{max}}$) of the ER:
$$\text{Interface}_i(x, y) = \text{Boundary}_i(x, y) \cap [D_{\text{ER}}(x, y) \leq d_{\text{max}}]$$

The physical interface length in nanometers is the sum of these pixels scaled by the pixel size:
$$I_i = \left( \sum_{x, y} \text{Interface}_i(x, y) \right) \times \mathrm{pixel\_size}$$

#### Step 4: Average Interface Distance
If the interface length is greater than zero ($I_i > 0$), we calculate the average distance to the ER along the contact interface:
$$\bar{d}_i = \frac{1}{N_{\text{interface}}} \sum_{(x,y) \in \text{Interface}_i} D_{\text{ER}}(x, y)$$
*(where $N_{\text{interface}}$ is the number of pixels in the interface)*

#### Step 5: ERMICC Calculation
The ER-Mitochondria Contact Coefficient (ERMICC) is calculated as:
$$\text{ERMICC}_i = \frac{I_i}{P_i \times \bar{d}_i}$$

If there is no contact ($I_i = 0$), or if the denominator is $0$, the coefficient is defined as:
$$\text{ERMICC}_i = 0$$

#### Step 6 Major axis length and Feret diameter (longest distance)

| Feature | `axis_major_length` | `feret_diameter_max` |
| ------- | ------------------- | -------------------- |
| Definition | The major axis length of an ellipse that has the same normalized second central moments as the region. | The longest Euclidean distance between any two points on the region's convex hull contour. |
| Measurement Method | Calculated via ellipse fitting. | Calculated directly on the convex hull boundary points. |
| Shape Dependency | Heavily influenced by the overall spatial distribution of pixels (moments). | Captures the absolute extreme points of the object, unaffected by inner pixel distribution. |
| Outlier Sensitivity | Robust to minor pixel noise or small boundary irregularities. | Highly sensitive to small spikes, edges, or sharp protrusions on the object boundary. |
