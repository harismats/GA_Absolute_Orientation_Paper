# Code for Paper "A Geometric Algebra Solution to the Absolute Orientation Problem"

This repository contains the MATLAB code used in the paper *"A Geometric Algebra Solution to the Absolute Orientation Problem"*. The code demonstrates several approaches to the absolute orientation problem, including:

1. A traditional SVD-based method.
2. A Characteristic Multivector (CM) method implemented using the Clifford Multivector Toolbox.
3. MATLAB’s built-in Procrustes method.

The code computes the optimal rigid-body transformation (rotation and translation) between two 3D point clouds and compares the performance of these methods.

## Requirements

- **MATLAB** (version R2020a or later recommended)
- **Clifford Multivector Toolbox for MATLAB**  
  (Ensure that the toolbox is added to the MATLAB path via `addpath`.)
- Sample point cloud data (e.g., `apple.ply`) – available in the repository or downloadable from the sources mentioned in the paper.

## File Structure

- **`full_implementation.m`**  
  Main script that:
  - Loads a point cloud.
  - Applies random rotations and translations.
  - Registers the transformed point cloud using the three methods.
  - Computes error metrics (RMSE, Hausdorff distance, Angular error, etc.).
  - Generates and saves figures comparing the methods.

- **`rigidtform.m`**  
  Function to compute the rigid transformation between two point clouds.  
  *Methods available:*  
    1. SVD-based method  
    2. Characteristic Multivector (CM) method  
    3. MATLAB’s Procrustes method

- **`matrix_F_optimized.m`**  
  Computes the F matrix (cross-correlation) from the point clouds using the Clifford multivector toolbox.

- **`matrix_G_optimized.m`**  
  Computes the G matrix (auto-correlation) from the point clouds using the Clifford multivector toolbox.

- **`CM_FG_matrices.m`**  
  Recovers the rotor (rotation operator) from the F and G matrices using Clifford Algebra operations.

- **`hausdorff_distance.m`**  
  Computes the Hausdorff distance between two sets of points.

- **`mean_absolute_error.m`**  
  Computes the Mean Absolute Error between two datasets.

## How to Run

1. **Set Up MATLAB Environment:**
   - Open MATLAB and change the current folder to the repository folder.
   - Ensure that the Clifford Multivector Toolbox and all function files are in your MATLAB path. (This is typically done via the `addpath` commands at the top of `full_implementation.m`.)

2. **Run the Main Script:**
   - Execute the file `full_implementation.m`.
   - The script will load the sample point cloud, apply random transformations, perform registration using the three methods, compute error metrics, and generate several comparison figures.
   - Figures will be saved automatically to the specified output folder.

3. **Review Figures and Output:**
   - Check the generated plots (with multi-line titles) in your output folder.
   - Use the figures and error metrics for further analysis or to replicate the results from the paper.

## Customization

- **Point Cloud Data:**  
  You can change the point cloud by modifying the file path in the main script.
- **Number of Test Cases:**  
  Adjust the variable `no_test_cases` in `full_implementation.m` to change the number of random transformation tests.
- **Output Folder:**  
  Update the variable `output_folder` with the desired path for saving figures.

## Contact

For any questions or issues with the code, please contact:  
**[Haris Matsantonis]**  
**[cm2128@cam.ac.uk]**

---

This README provides a detailed yet concise guide for users to understand and replicate the code from the paper.
