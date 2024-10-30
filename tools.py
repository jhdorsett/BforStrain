import numpy as np
from mesh2d import smooth2, tricon2
from scipy.spatial import Voronoi, Delaunay
import params
import pandas as pd
import geopandas as gpd
import disloc3d

def make_enu_crs(olon, olat, unit = None):
    # uses well known text (wkt) format to generate a CRS compatible with geodataframes
    # equivalent to a local ENU projection

    if unit == "kilometer": unit = "kilometre"
    if unit == "km": unit = "kilometre"
    if unit == "m": unit = "meter"
    
    if unit is None:
        unit = "kilometre"
        print("Assuming units in kilometers as default")
        print("Explicitly define unit as 'kilometer', 'km', 'meter', or 'm' to supress this warning.")

    elif unit not in ["kilometre", "meter"]:
        raise ValueError("Define unit as 'kilometer', 'km', 'meter', or 'm'.")
    
    # Set the unit value for the WKT
    unit_value = "1000" if unit == "kilometre" else "1"

    wkt_crs = f"""
    PROJCS["Local ENU",
        GEOGCS["WGS 84",
            DATUM["WGS_1984",
                SPHEROID["WGS 84",6378137,298.257223563]],
            PRIMEM["Greenwich",0],
            UNIT["degree",0.0174532925199433]],
        PROJECTION["Orthographic"],
        PARAMETER["latitude_of_origin",{olat}],
        PARAMETER["central_meridian",{olon}],
        PARAMETER["false_easting",0],
        PARAMETER["false_northing",0],
        UNIT["{unit}",{unit_value}]]
    """
    return wkt_crs

def llh2local(lons,lats,origin=params.origin):
    """
    Converts from longitude and latitude to local coordinates (x, y) in kilometers.
    """
    
    # WGS84 Ellipsoid Constants
    a = 6378137.0  # Semi-major axis in meters
    e = 0.0818191908426  # Eccentricity

    # Convert degrees to radians
    lons =  np.radians(lons)
    lats =  np.radians(lats)
    origin = np.radians(origin)
    
    # Calculate differences in latitudes and longitudes
    dlat = lats - origin[1]
    dlon = lons - origin[0]
    
    # Ensure longitude differences are within [-pi, pi]
    dlon = (dlon + np.pi) % (2 * np.pi) - np.pi
    
    # Compute average latitude and its sine/cosine
    lat_avg = (lats + origin[1]) / 2
    lat_avg_sin = np.sin(lat_avg)
    lat_avg_cos = np.cos(lat_avg)
    
    # Meridian radius of curvature (M)
    M = a * (1 - e**2) / (np.sqrt(1 - e**2 * lat_avg_sin**2))**3
    
    # Prime vertical radius of curvature (N)
    N = a / np.sqrt(1 - e**2 * lat_avg_sin**2)
    
    # Compute local coordinates in meters
    x = dlon * N * lat_avg_cos  # Easting
    y = dlat * M                # Northing
    
    # Convert to kilometers
    x_km = x / 1000
    y_km = y / 1000
    
    return np.array((x_km, y_km)).T

def local2llh(xs, ys, origin=params.origin):
    """
    Converts from local coordinates to longitude and latitude.
    """
    # Convert origin to radians
    origin = np.radians(origin)
    
    a = 6378137.0  # Semi-major axis in meters
    e = 0.0818191908426  # Eccentricity

    # Unpack local coordinates and convert to meters
    x = xs * 1000  # Easting in meters
    y = ys * 1000  # Northing in meters
    
    # Initial guess for latitude: use origin latitude
    lat = origin[1]
    tol = 1e-10
    max_iter = 10
    iter_count = 0
    diff = np.inf
    
    # Iteratively solve for latitude
    while diff > tol and iter_count < max_iter:
        lat_old = lat
        lat_avg = (lat + origin[1]) / 2
        lat_avg_sin = np.sin(lat_avg)
        
        # Meridian radius of curvature (M)
        M = a * (1 - e**2) / (np.sqrt(1 - e**2 * lat_avg_sin**2))**3
        
        # Update latitude using the y local coordinate
        lat = origin[1] + y / M
        diff = np.max(np.abs(lat - lat_old))
        iter_count += 1
    
    # Compute longitude using updated latitude
    lat_avg = (lat + origin[1]) / 2
    N = a / np.sqrt(1 - e**2 * np.sin(lat_avg)**2)
    lon = origin[0] + x / (N * np.cos(lat_avg))
    
    # Convert longitude to [-180, 180] degrees
    lon = (lon + np.pi) % (2 * np.pi) - np.pi
    
    # Convert back to degrees
    lon = np.degrees(lon)
    lat = np.degrees(lat)
    
    return np.array((lon, lat)).T

### Read in borders...
# naturalearthdata.com data downloaded from https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_admin_0_countries.geojson
borders = gpd.read_file('data/ne_10m_admin_0_countries.geojson')
###############
# Coordinate transforms defined using pyproj, called via wrapper functions
local_crs = make_enu_crs(olon=params.origin[0],olat=params.origin[1],unit='km')
borders_xy = borders.to_crs(local_crs)

def plot_borders(ax,lonlat=False,width=0.5):
     
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    if lonlat:
        borders.boundary.plot(ax=ax, linewidth=width, color='black')
    else:
        borders_xy.boundary.plot(ax=ax, linewidth=width, color='black')
        
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

class build_G_pointforce:
    """
        for a mesh object, builds body force green functions following sandwell and wessel
    
        mesh object requirements:
    
        nodes         2 by n numpy array  (created from geodetic data file+regular grid)
        tri_centroids 2 by n numpy array  (computed automatically in make_tri_centroids)
        
        from inputs:
        
        nu
    """
    def __init__(self,mesh):
        # for legibility in math below...

        print("")
        print("Computing Body Force Green's Functions")

        nu = params.nu

        x_diff = mesh.tri_centroids[:, 0][:, np.newaxis] - mesh.nodes[:, 0]
        y_diff = mesh.tri_centroids[:, 1][:, np.newaxis] - mesh.nodes[:, 1]
        r = np.sqrt(x_diff**2 + y_diff**2)
    
        Ue_x = (3 - nu) * np.log(r) + (1 + nu) * (y_diff**2) / r**2
        Un_y = (3 - nu) * np.log(r) + (1 + nu) * (x_diff**2) / r**2
        Ue_y = -(1 + nu) * x_diff * y_diff / r**2
        Un_x = Ue_y

        Exx_x = (3 - nu) * x_diff / r**2 - 2 * (1 + nu) * y_diff**2 * x_diff / r**4
        Eyy_y = (3 - nu) * y_diff / r**2 - 2 * (1 + nu) * x_diff**2 * y_diff / r**4

        Exx_y = -(1 + nu) * (-2 * x_diff**2 * y_diff / r**4 + y_diff / r**2)
        Eyy_x = -(1 + nu) * (-2 * y_diff**2 * x_diff / r**4 + x_diff / r**2)

        dq_dy = (3 - nu) * y_diff / r**2 + 2 * (1 + nu) * (-y_diff**3 / r**4 + y_diff / r**2)
        dp_dx = (3 - nu) * x_diff / r**2 + 2 * (1 + nu) * (-x_diff**3 / r**4 + x_diff / r**2)

        Exy_x = 0.5 * (Exx_y + dq_dy)
        Exy_y = 0.5 * (dp_dx + Eyy_x)
        
        self.Omega_x = 0.5 * (dq_dy - Exx_y)
        self.Omega_y = 0.5 * (Eyy_x - dp_dx)    

        self.GVe = np.concatenate([Ue_x, Ue_y], axis=1)
        self.GVn = np.concatenate([Un_x, Un_y], axis=1)

        self.GExx = np.concatenate([Exx_x, Exx_y], axis=1)
        self.GExy = np.concatenate([Exy_x, Exy_y], axis=1)
        self.GEyy = np.concatenate([Eyy_x, Eyy_y], axis=1)

class build_G_creep:
    """
        for a mesh object, discretizes given fault patches down to specified length, 
        builds creeping greens functions using Okada formulation (disloc3d)

        mesh object requirements where :
    
        nodes:         2 by n numpy array  (created from geodetic data file+regular grid)
        tri_centroids: 2 by n numpy array  (computed automatically in make_tri_centroids)
        pm:            7 by m numpy array  (fault patch model for m faults read in)
        Patch_id:      1 by m numpy array  
    """
    def __init__(self,mesh):
        
        def piecewise_Gcreep(faultnums,A1,A2):
                # helper function for 3 different geometric creep scenarios (fault is 1 patch, 2 patches, or 3+ patches)
            # combines G1 and G2 greens functions depending on fine mesh system
            # A1 = uniform slip on patch
            # A2 = slip tapering from 0 - 1 linearly 
            cnt = 0
            A = np.zeros(A1.shape)
            for k in range(len(faultnums)):
                # n provides the amount of patches (lines in the matrix) that are part of one continuous fault
                n = np.sum(np.floor(mesh.Patch_id) == faultnums[k])
                
                # assign slip based on configuration defined for piecewise_G above
                if n == 1:
                    A[:, cnt] = A1[:, cnt]

                # possibly can be removed and generalized but only 2 samples so I don't want to mess with it too much rn
                elif n == 2:
                    A[:, cnt] = A1[:, cnt] - A2[:, cnt]
                    A[:, cnt + 1] = A2[:, cnt] + A1[:, cnt + 1]
                
                else:
                    A[:, cnt] = A1[:, cnt] - A2[:, cnt]

                    for j in range(cnt + 1, cnt + n-1):
                        A[:, j] = A2[:, j - 1] + A1[:, j] - A2[:, j]

                    A[:, cnt+n-1] = A2[:, cnt+ n - 2] + A1[:, cnt + n-1] 
                
                cnt = cnt + n

            return(A)
    
        
        npatches = mesh.pm.shape[0]
        faultnums = np.unique(np.floor(mesh.Patch_id))
        nobs = len(mesh.ind)
        
        xobs = np.vstack((mesh.tri_centroids[:, :2].T, np.zeros((1, mesh.tri_centroids.shape[0]))))

        G1east_creep  = np.zeros((nobs, npatches))
        G1north_creep = np.zeros((nobs, npatches))
        G1Exx_creep   = np.zeros((nobs, npatches))
        G1Exy_creep   = np.zeros((nobs, npatches))
        G1Eyy_creep   = np.zeros((nobs, npatches))

        G2east_creep  = np.zeros((nobs, npatches))
        G2north_creep = np.zeros((nobs, npatches))
        G2Exx_creep   = np.zeros((nobs, npatches))
        G2Exy_creep   = np.zeros((nobs, npatches))
        G2Eyy_creep   = np.zeros((nobs, npatches))

        print("Computing Creepings Green's functions calculation")
        past = 0
        for k in range(npatches):
            # Divide into small segments
            nhe = int(np.ceil(mesh.pm[k, 0] / params.refine))
            pf = patchfault(mesh.pm[k, :], nhe, 1)
            
            # run disoc3d on each fine grained fault patch and increment over each coarse piece
            # k is looping over large patches, sum up displacement for each large patch in small patch (J)
            # G2 are additional GFs used below e.g. piecewise_G that apply a taper to slip on each patch
            
            for j in range(nhe):
                taper = j / nhe
                # create disloc3d fault plane with 1 unit of right lateral along strike slip
                m1 = np.hstack([pf[j, :], -1, 0, 0])                
                
                # dont store stress value since we don't use it                
                U1, D, _ = disloc3d.dc3d_wrapper(pf[j,:],[-1,0,0],xobs.T)
                # update G1 matrices
                G1east_creep[:, k] += U1[:, 0]
                G1north_creep[:, k] += U1[:, 1]
                
                # commented out lines for indexing if return_2d=False or unspecified
                G1Exy_creep[:, k] += 0.5 * (D[:, 1] + D[:, 3])
                G1Eyy_creep[:, k] += D[:, 4]
                G1Exx_creep[:, k] += D[:, 0]

                # update G2 matrices 
                G2east_creep[:, k]  += taper * U1[:, 0]
                G2north_creep[:, k] += taper * U1[:, 1]

                # commented out lines for indexing if return_2d=False or unspecified
                G2Exy_creep[:, k] += taper * 0.5 * (D[:, 1] + D[:, 3])
                G2Eyy_creep[:, k] += taper * D[:, 4]
                G2Exx_creep[:, k] += taper * D[:, 0]

                del U1, D

            pct_done = round(k/npatches*100)
            if round(pct_done,-1) > past:
                print(" ",round(pct_done,-1),"% completed")
                past = round(pct_done,-1)

        self.Geast=piecewise_Gcreep(faultnums,G1east_creep,G2east_creep)
        self.Gnorth=piecewise_Gcreep(faultnums,G1north_creep,G2north_creep)
        self.GExx=piecewise_Gcreep(faultnums,G1Exx_creep,G2Exx_creep)
        self.GExy=piecewise_Gcreep(faultnums,G1Exy_creep,G2Exy_creep)
        self.GEyy=piecewise_Gcreep(faultnums,G1Eyy_creep,G2Eyy_creep)
     
def build_mesh(mesh):
    
    def make_patches_creep(mesh):

        def discretize_patches(mesh):
            SegEnds = mesh.SegEnds
            creeping_faults = mesh.creeping_faults

            # Initialize variables
            PatchEnds = np.empty((0, 4), dtype=float)
            PatchCreepRates = np.empty((0,1), dtype=float)
            Patch_id = np.empty((0,1), dtype=int)

            # discretize 
            for k in range(SegEnds.shape[0]):
                patchlength = np.sqrt((SegEnds[k, 2] - SegEnds[k, 0]) ** 2 + (SegEnds[k, 3] - SegEnds[k, 1]) ** 2)
                numpatch = int(np.ceil(patchlength / params.patchL))
                xs = np.linspace(SegEnds[k, 0], SegEnds[k, 2], numpatch + 1)
                ys = np.linspace(SegEnds[k, 1], SegEnds[k, 3], numpatch + 1)
                PatchEnds = np.vstack((PatchEnds, np.column_stack((xs[:-1], ys[:-1], xs[1:], ys[1:]))))
                PatchCreepRates = np.vstack((PatchCreepRates, np.full((numpatch, 1), creeping_faults[k, 5])))
                Patch_id = np.vstack((Patch_id, np.full((numpatch,1), creeping_faults[k, 5])))

            mesh.PatchCreepRates = PatchCreepRates
            mesh.PatchEnds = PatchEnds
            mesh.Patch_id = Patch_id

        creeping_faults = mesh.creeping_faults

        #Convert lon,lat to x,y
        x1 = llh2local(creeping_faults[:,1],creeping_faults[:,2]) 
        x2 = llh2local(creeping_faults[:,3],creeping_faults[:,4])

        # crop to bounding box, since each segment contains 2 points 
        _, keeppatch1 = mesh.bounds.crop(x1,True)
        _, keeppatch2 = mesh.bounds.crop(x2,True)
        keeppatch = keeppatch1 | keeppatch2
        
        SegEnds = np.column_stack((x1,x2))
        mesh.SegEnds = SegEnds[keeppatch]
        mesh.creeping_faults = creeping_faults[keeppatch]
        
        # discretizes each creeping fault segment down to the length specified in inputs file
        # creates PatchEnds, PatchCreepRates, and Patch_id
        discretize_patches(mesh)

        segends1 = mesh.PatchEnds[:, :2]
        segends2 = mesh.PatchEnds[:, 2:]

        num_rows = segends1.shape[0] + segends2.shape[0]
        node_creep = np.zeros((num_rows, 2))

        node_creep[::2] = segends1
        node_creep[1::2] = segends2

        # Create 'edge_creep' using NumPy
        edge_creep = np.column_stack((np.arange(0, num_rows, 2), np.arange(1, num_rows, 2)))

        # Convert the NumPy array to a pandas DataFrame because the numpy unique function automatically sorts
        # we need to preserve the order
        df = pd.DataFrame(np.round(node_creep, 4), columns=['x', 'y'])
        ic, unique_vals = pd.factorize(df.apply(tuple, axis=1))

        mesh.node_creep = np.array(unique_vals.tolist())

        new_edge = edge_creep.copy()
        for k in range(len(ic)):
            new_edge[edge_creep == (k)] = ic[k]

        mesh.edge_creep = new_edge
        
    def refine(mesh):
        # Compute the Voronoi diagram
        vor = Voronoi(mesh.xy_gps)
        # Extract the vertices of the Voronoi diagram
        vx = vor.vertices[:, 0]
        vy = vor.vertices[:, 1]
        _, unique_indices = np.unique(vx, return_index=True) 
        vx = vx[unique_indices]
        vy = vy[unique_indices]
    
        nodes_refine = np.column_stack((vx, vy))

        # crop it using our box
        nodes_refine = mesh.bounds.crop(nodes_refine)

        if params.refine_mesh == 4:
            mesh.nodes = np.vstack((mesh.nodes, nodes_refine))
        elif params.refine_mesh == 3:
            mesh.nodes = np.vstack((mesh.nodes, nodes_refine[::2]))
        elif params.refine_mesh == 2:
            mesh.nodes = np.vstack((mesh.nodes, nodes_refine[::3]))
        elif params.refine_mesh == 1:
            mesh.nodes = np.vstack((mesh.nodes, nodes_refine[::5]))

    # using bounding box, create meshgrid
    mesh.nodes=mesh.bounds.make_grid()


    # crop geodetic data down to bounding box extent
    mesh.xy_gps, ind_gps = mesh.bounds.crop(mesh.xy_gps,True)
    mesh.vel.Ve=mesh.vel.Ve[ind_gps]
    mesh.vel.Vn=mesh.vel.Vn[ind_gps]
    mesh.vel.Sige=mesh.vel.Sige[ind_gps]
    mesh.vel.Sign=mesh.vel.Sign[ind_gps]

    # combine bounding box and geodetic observation points for initial nodes
    mesh.nodes = np.vstack((mesh.nodes, mesh.xy_gps))

    # if desired, refine mesh by computing voronoi diagram 
    # and upscaling mesh with these points
    if bool(params.refine_mesh) == True: refine(mesh)
        
    # generate initial guess for delaunay triangles. maybe switch out with tridel2?
    triDel=Delaunay(mesh.nodes, qhull_options = 'Qt Qbb Qc') #qhull settings consistent with matlab
    mesh.tri=triDel.simplices

    # save fixed edges of mesh
    edge = tricon2(mesh.tri)[0]
    edge_bnd = edge[:,3] < 1
    mesh.edges = edge[edge_bnd,0:2]

    # if we have a set of creeping fault traces that we want to incorporate into the mesh:
    if mesh.creeping:
        make_patches_creep(mesh)        
        # update nodes and list of fixed edges 
        mesh.edges = np.vstack((mesh.edges, mesh.edge_creep+mesh.nodes.shape[0]))
        mesh.nodes = np.vstack((mesh.nodes, mesh.node_creep))

    else:
        mesh.SegEnds = []
        mesh.PatchEnds = []
        mesh.PatchCreepRates = []
        mesh.Patch_id = []
        mesh.edge_creep = []
        mesh.node_creep = []

    # now that we have generated an initial nodeset, fixed edge set, and triangulation,
    # smooth it all
    mesh = smooth2(mesh) 
    mesh.elts=mesh.nodes[mesh.tri]

def patchfault(m,i,j):
    dip = m[3] * np.pi / 180
    strike = -m[4] * np.pi / 180
    sin_dip = np.sin(dip)
    cos_dip = np.cos(dip)
    iw = m[0] / i
    jw = m[1] / j
    is_ = np.arange(1, i + 1) 
    js = np.arange(1, j + 1)

    n = i * j
    #c1 = -m[1] * cos_dip
    #c2 = 0.5 * (m[0] + iw)
    #c3 = m[2] - j * jw * sin_dip

    # Calculate midpoints, depths of each patch
    p = np.outer(cos_dip * (jw * js - m[1]), np.ones(i))
    q = np.outer(np.ones(j), (iw * is_) - 0.5 * (m[0] + iw))
    r = np.outer(m[2] - jw * sin_dip * (j - js), np.ones(i))
    mp = np.column_stack((p.flatten(), q.flatten(), r.flatten()))

    # Adjust midpoints for strike
    R = np.array([[np.cos(strike), -np.sin(strike), 0],
                  [np.sin(strike), np.cos(strike), 0],
                  [0, 0, 1]])
    mp = np.dot(mp, R.T)

    # Adjust midpoints for offset from origin
    mp[:, 0] += m[5]
    mp[:, 1] += m[6]

    # Form patch-models
    pm = np.zeros((n, 7))
    pm[:, 0] = np.ones(n) * iw
    pm[:, 1] = np.ones(n) * jw
    pm[:, 2] = mp[:, 2]
    pm[:, 3:5] = np.ones(n).reshape(-1, 1) * m[3:5]
    pm[:, 5:7] = mp[:, 0:2]

    return(pm)

def make_triangular_patch_stuff(tri, p):
    strikevec_faces = []
    strike_faces = []
    dipvec_faces = []
    dip_faces = []
    centroids_faces = []
    normal_faces = []
    area_faces = []

    for j in range(tri.shape[0]):
        temp1 = [p]
        temp2 = [tri[j, :]]

        vec1 = p[temp2[0][0], :] - p[temp2[0][1], :]
        vec2 = p[temp2[0][2], :] - p[temp2[0][1], :]
        cross_face = np.cross(vec1, vec2)
        veclength = np.linalg.norm(cross_face)
        normal = cross_face / veclength
        strikevec = np.array([1, -normal[0] / normal[1], 0])
        strikevec = strikevec / np.linalg.norm(strikevec)
        dipvec = np.cross(normal, strikevec)

        if dipvec[2] > 0:
            dipvec = -dipvec
        if normal[2] < 0:
            normal = -normal
        strikevec = np.cross(normal, dipvec)

        normal_faces.append(normal)
        strikevec_faces.append(strikevec)
        dipvec_faces.append(dipvec)
        strike_faces.append(90 - np.arctan2(strikevec[1], strikevec[0]) * 180 / np.pi)
        dip_faces.append(np.abs(np.arctan(dipvec[2] / np.sqrt(dipvec[0] ** 2 + dipvec[1] ** 2)) * 180 / np.pi))
        centroids_faces.append([np.mean(temp1[0][temp2[0], 0]), np.mean(temp1[0][temp2[0], 1]), np.mean(temp1[0][temp2[0], 2])])
        area_faces.append(0.5 * np.abs(np.linalg.norm(np.cross(vec1, vec2))))

    patch_stuff = {
        'strikevec_faces': strikevec_faces,
        'strike_faces': strike_faces,
        'dipvec_faces': dipvec_faces,
        'dip_faces': dip_faces,
        'centroids_faces': centroids_faces,
        'normal_faces': normal_faces,
        'area_faces': area_faces
    }
    return patch_stuff
