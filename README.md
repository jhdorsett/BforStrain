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

Then you can execute the code using the data provided using ``run.ipynb``.
