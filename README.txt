These Matlab codes compute strain rates on a triangle mesh using a geodetic velocity field. The method is based on the manuscript K.Johnson (in prep). Body-forces in a thin elastic plate are used to compute a velocity and strain rate field computed at the centroids of the triangles in the mesh.  The inversion solves for the distribution of body forces (at the triangle nodes) that best-fits the geodetic velocity field and slip rates across creeping faults (if included).  You can specify a range of smoothing values (called beta) that minimize the magnitude of the body forces. 

An example fault creep rate file and geodetic data set is included for testing purposes. 

You need to run the scripts in the following order:
1. setup_mesh.m
2. build_bodyforce_Greens.m
3. invert_bodyforce.m
4. plot_inversion_results.m

Note:  outputs from 1 and 2 and 3 are saved in mat files, so you do not need to rerun all steps every time

1. setup_mesh.m creates the triangular mesh using mesh2d. You need to specify a number of parameters in the Input Section at the top of the script.  Here you specify the origin of the Cartesian coordinate system, load the GPS data file, the creeping fault data file (if you have one), and the nominal node spacing and mesh domain boundaries. You can also experiment with different levels of mesh refinement by change the value of 'refine_mesh'.

2. build_bodyforce_Greens.m computes the Greens Functions which give strains and velocities at centroids of triangles due to unit body force vectors. The script also computes fault creep Greens Functions if you have specified creeping fault segments in the setup_mesh.m script.  The only input you need to specify in this script is the name of the mat file generated in step 1. 

3.  invert_bodyforce.m computes the inversion of your data for the distribution of body forces (and fault creep if that applies).  You can specify a single value or a range of values for smoothing parameter, beta.  The inversion loops over all beta values and computes 'num' realizations of the strain rate and velocity fields for each beta.  You can specify the beta values, num, and the relative weight placed on fitting the creep rate data for creeping faults (if this applies).  

4. plot_inversion_results.m plots the result of the inversion. There are no inputs required to run this script as long as all the parameters and results are currently loaded into the Matlab workspace.

