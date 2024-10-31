# BforStrain
These codes compute strain rates on a trianglular 2D mesh using a geodetic velocity field. The method is based on the paper by K. Johnson (2024) available at https://doi.org/10.1029/2023JB027472 . Body-forces in a thin elastic plate are used to compute a velocity and strain rate field computed at the centroids of the triangles in the mesh.  The inversion solves for the distribution of body forces (at the triangle nodes) that best-fits the geodetic velocity field.  You can specify a range of smoothing values (called beta) that minimize the magnitude of the body forces. There is also an optional second step that further minimizes strain rates below a specified threshold.

Paramters can be specified in params.py. The data used in Johnson (2024) from Zeng (2022) https://doi.org/10.1785/0220220180 is currently loaded, along with creeping fault traces from the 2023 US National Seismic Hazard Model. 

To run the demo code, use the code below to download the repository and initialize a new environment

```
git clone https://github.com/jhdorsett/BforStrain
cd BforStrain
conda env create -f environment.yml
conda activate BforStrain
```

Then the code can be executed with the data provided using ``run.ipynb``.

Conceptual project structure:
```
BforStrain
   | bForStrain.py
   |   |-- Mesh
   |   |   |-- create_gps_obs
   |   |   |-- create_creeping_segs
   |   |   |-- create_creeping_segs
   |   |   |-- crop
   |   |   |-- save
   |   |   |-- plotting_functions
   |   |   |-- coordinate_transforms
   |   |-- Inversion
   |   |   |-- bodyforce_greens_functions
   |   |   |-- creeping_greens_functions
   |   |   |-- prepare_inversion
   |   |   |-- invert
   |   |   |-- generate_uncertainty
   |   |   |-- post_process_results
   | disloc3d.py
   |   |-- dc3d
   |   |-- dc3d_wrapper
   | mesh2d.py
   |   |-- smooth2d
   |   |-- helper_functions
   | tools.py
   |   |-- helper_functions
   | data/
   |   |-- geodetic veolcities
   |   |-- creeping faults
   |   |-- coastline data (plotting)
```
