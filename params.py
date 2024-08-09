
### Alpides configuration

gps_file = 'data/GNSS_velocities_alpides.txt'    
creeping_file = []
origin = [32.5, 35] 
lon_range = [-15,142]
lat_range = [5, 60]
nom_node_spacing = 400.

## West US configuration

gps_file = 'data/Zeng_vels_augmented_PA_JDF.txt'
creeping_file = 'data/creeping_faults.txt'  
origin = [-120,34] 
lon_range = [-127, -96]
lat_range = [26, 54]
nom_node_spacing = 250.

# mesh decimation factor (speeds up runtime, useful for debugging large scenarios)
# set = 1 for no decimation
decimate = 4

# creeping patch length
patchL = 15.

#option to refine mesh in vicinity of GPS data
# 0 = no refinement, 1-4 provides some
refine_mesh = 0



refine = 1
nu = 0.25
Gshear = 1

betas = [40]


uncertainty = False
# number of realizations of strain rate for each beta value
num = 50
# relative weight on fitting creep rate data (creeping faults)
Wcreep = 1
# optional two-step minimization of strain rates below threshold value
# set twostep = True or False
twostep = False
# relative weight (gamma) on minimizing strain rates below strain_threshold
# (micro-strain per year)
gamma = 400
# 
strain_threshold = 9e-3  # micro-strain per year
minimize_xbound = [-600, 1000]
minimize_ybound = [-800, 1000]