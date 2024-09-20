import numpy as np
from mesh2d import smooth2, tricon2
from scipy.spatial import Voronoi, Delaunay
import params
import pandas as pd
import pyproj
import geopandas as gpd

###############
# Coordinate transforms defined using pyproj, called via wrapper functions
local_crs = pyproj.CRS.from_proj4(f"+proj=tmerc +lat_0={params.origin[1]} +lon_0={params.origin[0]} +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")
llh_to_local = pyproj.Transformer.from_crs(pyproj.CRS("EPSG:4326"), local_crs, always_xy=True)
local_to_llh = pyproj.Transformer.from_crs(local_crs, pyproj.CRS("EPSG:4326"), always_xy=True)

def llh2local(lons,lats): 
     return(np.array(llh_to_local.transform(lons,lats)).T)
def local2llh(xs,ys): 
     return(np.array(local_to_llh.transform(xs,ys)).T)

### Read in borders...
# naturalearthdata.com data downloaded from https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_admin_0_countries.geojson
borders = gpd.read_file('data/ne_10m_admin_0_countries.geojson')
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
                U1, D = disloc3D_wrapper(m1, xobs, return_stress=False,return_2d=True)
                
                # update G1 matrices
                G1east_creep[:, k] += U1[0, :]
                G1north_creep[:, k] += U1[1, :]
                
                # commented out lines for indexing if return_2d=False or unspecified
                #G1Exy_creep[:, k] += 0.5 * (D[1, :] + D[3, :])
                #G1Eyy_creep[:, k] += D[4, :]
                G1Exx_creep[:, k] += D[0, :]
                G1Exy_creep[:, k] += 0.5 * (D[1, :] + D[2, :])
                G1Eyy_creep[:, k] += D[3, :]

                # update G2 matrices 
                G2east_creep[:, k]  += taper * U1[0, :]
                G2north_creep[:, k] += taper * U1[1, :]

                # commented out lines for indexing if return_2d=False or unspecified
                #G2Exy_creep[:, k] += taper * 0.5 * (D[1, :] + D[3, :])
                #G2Eyy_creep[:, k] += taper * D[4, :]
                G2Exx_creep[:, k] += taper * D[0, :]
                G2Exy_creep[:, k] += taper * 0.5 * (D[1, :] + D[2, :])
                G2Eyy_creep[:, k] += taper * D[3, :]

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

def checkInputs(disloc,coordinates):
        if disloc.c < disloc.W * np.sin(disloc.delta) and np.array_equal(coordinates[2, :], -np.abs(coordinates[2, :])):
                print('warning: physically impossible')
                return(False)
        elif disloc.c >= disloc.W * np.sin(disloc.delta) and np.array_equal(coordinates[2, :], -np.abs(coordinates[2, :])):
                return(True)
        else:
                print('warning: All z should be negative.')
                return(False)
     
def disloc3D(disloc,coordinates):
    X, Y, z = coordinates[:3, :]
    x,y=disloc.coordTrans(X,Y) #rorate into local system
    L = disloc.L
    W = disloc.W
    c = disloc.c
    delta = disloc.delta #dip in rad
    angle_Str = disloc.strike # clockwise is positive
    [Xc, Yc] = disloc.centroids #rupture centroids
    [slip_str, slip_dip, tensile] = disloc.slip

    nu=disloc.nu
    Gshear=disloc.Gshear
    youngs, alpha = disloc.elastic() # method in elastic class to compute these from nu and Gshear defined in class

    # integrating
    d = c - z
    p = y * np.cos(delta) + d * np.sin(delta)

    xi = np.array([x, x, x - L, x - L])
    eta = np.array([p, p - W, p, p - W])

    q = np.ones((4,x.size)) * y * np.sin(delta) - np.ones((4,x.size)) * d * np.cos(delta)

    R = np.sqrt(xi**2 + eta**2 + q**2)

    y_ = eta * np.cos(delta) + q * np.sin(delta)
    d_ = eta * np.sin(delta) - q * np.cos(delta)
    c_ = d_ + np.ones((4,x.size)) * z

    # For displacement
    X11 = 1 / (R * (R + xi))
    X32 = (2 * R + xi) / (R**3 * (R + xi)**2)
    X53 = (8 * R**2 + 9 * R * xi + 3 * xi**2) / (R**5 * (R + xi)**3)

    Y11 = 1 / (R * (R + eta))
    Y32 = (2 * R + eta) / (R**3 * (R + eta)**2)
    Y53 = (8 * R**2 + 9 * R * eta + 3 * eta**2) / (R**5 * (R + eta)**3)

    h = q * np.cos(delta) - np.ones((4,x.size)) * z
    Z32 = np.sin(delta) / R**3 - h * Y32
    Z53 = 3 * np.sin(delta) / R**5 - h * Y53

    Y0 = Y11 - xi**2 * Y32
    Z0 = Z32 - xi**2 * Z53

    # Selecting a right root for theta
    qsign = np.sign(q)
    theta = np.arctan2(xi * eta, np.abs(q) * R)
    theta = qsign * theta

    X = np.sqrt(xi**2 + q**2)

    if np.abs(np.cos(delta)) < 0.000001:
        I3 = 1/2 * (eta / (R + d_) + y_ * q / ((R + d_)**2) - np.log(R + eta))
        I4 = 1/2 * (xi * y_ / ((R + d_)**2))
    else:
        I3 = 1/np.cos(delta) * y_ / (R + d_) - 1 / np.cos(delta)**2 * (
                np.log(R + eta) - np.sin(delta) * np.log(R + d_))
        I4 = np.sin(delta) / np.cos(delta) * xi / (R + d_) + 2 / (np.cos(delta)**2) * np.arctan2(
                eta * (X + q * np.cos(delta)) + X * (R + X) * np.sin(delta),
                xi * (R + X) * np.cos(delta))

    I1 = -(xi / (R + d_)) * np.cos(delta) - I4 * np.sin(delta)
    I2 = np.log(R + d_) + I3 * np.sin(delta)

    D11 = 1 / (R * (R + d_))

    if np.abs(np.cos(delta)) < 0.000001:
        K1 = (xi * q) / (R + d_) * D11
        K3 = np.sin(delta) / (R + d_) * (xi**2 * D11 - 1)
    else:
        K1 = xi / np.cos(delta) * (D11 - Y11 * np.sin(delta))
        K3 = 1 / np.cos(delta) * (q * Y11 - y_ * D11)

    K2 = 1 / R + K3 * np.sin(delta)
    K4 = xi * Y11 * np.cos(delta) - K1 * np.sin(delta)

    J5 = -(d_ + y_**2 / (R + d_)) * D11
    J2 = xi * y_ / (R + d_) * D11

    if np.abs(np.cos(delta)) < 0.000001:
        J6 = -y_ / (R + d_)**2 * (xi**2 * D11 - 1/2)
        J3 = -xi / (R + d_)**2 * (q**2 * D11 - 1/2)
    else:
        J6 = 1 / np.cos(delta) * (K3 - J5 * np.sin(delta))
        J3 = 1 / np.cos(delta) * (K1 - J2 * np.sin(delta))

    J1 = J5 * np.cos(delta) - J6 * np.sin(delta)
    J4 = -xi * Y11 - J2 * np.cos(delta) + J3 * np.sin(delta)

    # ki
    E = np.sin(delta) / R - y_ * q / R**3
    F = d_ / R**3 + xi**2 * Y32 * np.sin(delta)
    G = 2 * X11 * np.sin(delta) - y_ * q * X32
    H = d_ * q * X32 + xi * q * Y32 * np.sin(delta)
    P = np.cos(delta) / R**3 + q * Y32 * np.sin(delta)
    Q = 3 * c_ * d_ / R**5 - (np.ones((4, 1)) * z * Y32 + Z32 + Z0) * np.sin(delta)

    # li
    E_ = np.cos(delta) / R + d_ * q / R**3
    F_ = y_ / R**3 + xi**2 * Y32 * np.cos(delta)
    G_ = 2 * X11 * np.cos(delta) + d_ * q * X32
    H_ = y_ * q * X32 + xi * q * Y32 * np.cos(delta)
    P_ = np.sin(delta) / R**3 - q * Y32 * np.cos(delta)
    Q_ = (3 * c_ * y_) / R**5 + q * Y32 - (np.ones((4, 1)) * z * Y32 + Z32 + Z0) * np.cos(delta)

    # strike slip 
    if slip_str != 0:
        # displacement
        #uA
        Su1A = theta/2 + alpha / 2 * xi * q * Y11
        Su2A = alpha / 2 * q / R
        Su3A = (1 - alpha) / 2 * np.log(R + eta) - alpha / 2 * q**2 * Y11
        #uB
        Su1B = -xi * q * Y11 - theta - (1 - alpha) / alpha * I1 * np.sin(delta)
        Su2B = -q / R + (1 - alpha) / alpha * y_ / (R + d_) * np.sin(delta)
        Su3B = q**2 * Y11 - (1 - alpha) / alpha * I2 * np.sin(delta)
        #uC
        Su1C = (1 - alpha) * xi * Y11 * np.cos(delta) - alpha * xi * q * Z32
        Su2C = (1 - alpha) * (np.cos(delta) / R + 2 * q * Y11 * np.sin(delta)) - alpha * c_ * q / R**3
        Su3C = (1 - alpha) * q * Y11 * np.cos(delta) - alpha * (c_ * eta / R**3 - np.ones((4, 1)) * z * Y11 + xi**2 * Z32)        
        # displacement gradient
        #jA
        Sj1A = -(1 - alpha) / 2 * q * Y11 - alpha / 2 * (xi**2) * q * Y32
        Sj2A = - alpha / 2 * xi * q / R**3
        Sj3A = (1 - alpha) / 2 * xi * Y11 + alpha / 2 * xi * q**2 * Y32
        #jB
        Sj1B = xi**2 * q * Y32 - (1 - alpha) / alpha * J1 * np.sin(delta)
        Sj2B = xi * q / R**3 - (1 - alpha) / alpha * J2 * np.sin(delta)
        Sj3B = -xi * q**2 * Y32 - (1 - alpha) / alpha * J3 * np.sin(delta)
        #jC
        Sj1C = (1 - alpha) * Y0 * np.cos(delta) - alpha * q * Z0
        Sj2C = -(1 - alpha) * xi * (np.cos(delta) / R**3 + 2 * q * Y32 * np.sin(delta)) + alpha * (
                        3 * c_ * xi * q) / R**5
        Sj3C = -(1 - alpha) * xi * q * Y32 * np.cos(delta) + alpha * xi * (
                        (3 * c_ * eta) / R**5 - np.ones((4, 1)) * z * Y32 - Z32 - Z0)
        #kA
        Sk1A = (1 - alpha) / 2 * xi * Y11 * np.sin(delta) + d_ / 2 * X11 + alpha / 2 * xi * F
        Sk2A = alpha / 2 * E
        Sk3A = (1 - alpha) / 2 * (
                        np.cos(delta) / R + q * Y11 * np.sin(delta)) - alpha / 2 * q * F
        #kB
        Sk1B = -xi * F - d_ * X11 + (1 - alpha) / alpha * (xi * Y11 + J4) * np.sin(delta)
        Sk2B = -E + (1 - alpha) / alpha * (1 / R + J5) * np.sin(delta)
        Sk3B = q * F - (1 - alpha) / alpha * (q * Y11 - J6) * np.sin(delta)
        #kC
        Sk1C = -(1 - alpha) * xi * P * np.cos(delta) - alpha * xi * Q
        Sk2C = 2 * (1 - alpha) * (d_ / R**3 - Y0 * np.sin(delta)) * np.sin(delta) - y_ / R**3 * np.cos(
                delta) - alpha * ((c_ + d_) / R**3 * np.sin(delta) - eta / R**3 - 3 * c_ * y_ * q / R**5)
        Sk3C = -(1 - alpha) * q / R**3 + (y_ / R**3 - Y0 * np.cos(delta)) * np.sin(
                delta) + alpha * ((c_ + d_) / R**3 * np.cos(delta) + 3 * c_ * d_ * q / R**5 - (
                        Y0 * np.cos(delta) + q * Z0) * np.sin(delta))
        #lA
        Sl1A = (1 - alpha) / 2 * xi * Y11 * np.cos(delta) + y_ / 2 * X11 + alpha / 2 * xi * F_
        Sl2A = alpha / 2 * E_
        Sl3A = -(1 - alpha) / 2 * (
                        np.sin(delta) / R - q * Y11 * np.cos(delta)) - alpha / 2 * q * F_
        #lB
        Sl1B = -xi * F_ - y_ * X11 + (1 - alpha) / alpha * K1 * np.sin(delta)
        Sl2B = -E_ + (1 - alpha) / alpha * y_ * D11 * np.sin(delta)
        Sl3B = q * F_ + (1 - alpha) / alpha * K2 * np.sin(delta)
        #lC
        Sl1C = (1 - alpha) * xi * P_ * np.cos(delta) - alpha * xi * Q_
        Sl2C = 2 * (1 - alpha) * (y_ / R**3 - Y0 * np.cos(delta)) * np.sin(delta) + d_ / R**3 * np.cos(
                delta) - alpha * ((c_ + d_) / R**3 * np.cos(delta) + 3 * c_ * d_ * q / R**5)
        Sl3C = (y_ / R**3 - Y0 * np.cos(delta)) * np.cos(delta) - alpha * (
                        (c_ + d_) / R**3 * np.sin(delta) - 3 * c_ * y_ * q / R**5 - Y0 * np.sin(
                delta)**2 + q * Z0 * np.cos(delta))
        # displacement
        # u1A_
        Su1A_ = theta / 2 + alpha / 2 * xi * q * Y11
        Su2A_ = alpha / 2 * q / R
        Su3A_ = (1 - alpha) / 2 * np.log(R + eta) - alpha / 2 * q**2 * Y11
        # displacement gradient
        # jA_
        Sj1A_ = -(1 - alpha) / 2 * q * Y11 - alpha / 2 * xi**2 * q * Y32
        Sj2A_ = - alpha / 2 * xi * q / R**3
        Sj3A_ = (1 - alpha) / 2 * xi * Y11 + alpha / 2 * xi * q**2 * Y32
        # kA_
        Sk1A_ = (1 - alpha) / 2 * xi * Y11 * np.sin(delta) + d_ / 2 * X11 + alpha / 2 * xi * F
        Sk2A_ = alpha / 2 * E
        Sk3A_ = (1 - alpha) / 2 * (np.cos(delta) / R + q * Y11 * np.sin(delta)) - alpha / 2 * q * F
        # lA_
        Sl1A_ = (1 - alpha) / 2 * xi * Y11 * np.cos(delta) + y_ / 2 * X11 + alpha / 2 * xi * F_
        Sl2A_ = alpha / 2 * E_
        Sl3A_ = - (1 - alpha) / 2 * (np.sin(delta) / R - q * Y11 * np.cos(delta)) - alpha / 2 * q * F_

        # displacement
        Sux = 1 / (2 * np.pi) * slip_str * (Su1A - Su1A_ + Su1B + np.ones((4, 1)) * z * Su1C)
        Suy = 1 / (2 * np.pi) * slip_str * (
                (Su2A - Su2A_ + Su2B + np.ones((4, 1)) * z * Su2C) * np.cos(delta) -
                (Su3A - Su3A_ + Su3B + np.ones((4, 1)) * z * Su3C) * np.sin(delta))
        Suz = 1 / (2 * np.pi) * slip_str * (
                (Su2A - Su2A_ + Su2B - np.ones((4, 1)) * z * Su2C) * np.sin(delta) +
                (Su3A - Su3A_ + Su3B - np.ones((4, 1)) * z * Su3C) * np.cos(delta))

        # displacement gradients
        Sduxdx = 1 / (2 * np.pi) * slip_str * (Sj1A - Sj1A_ + Sj1B + np.ones((4, 1)) * z * Sj1C)
        Sduydx = 1 / (2 * np.pi) * slip_str * (
                (Sj2A - Sj2A_ + Sj2B + np.ones((4, 1)) * z * Sj2C) * np.cos(delta) -
                (Sj3A - Sj3A_ + Sj3B + np.ones((4, 1)) * z * Sj3C) * np.sin(delta))
        Sduzdx = 1 / (2 * np.pi) * slip_str * (
                (Sj2A - Sj2A_ + Sj2B - np.ones((4, 1)) * z * Sj2C) * np.sin(delta) +
                (Sj3A - Sj3A_ + Sj3B - np.ones((4, 1)) * z * Sj3C) * np.cos(delta))

        Sduxdy = 1 / (2 * np.pi) * slip_str * (Sk1A - Sk1A_ + Sk1B + np.ones((4, 1)) * z * Sk1C)
        Sduydy = 1 / (2 * np.pi) * slip_str * (
                (Sk2A - Sk2A_ + Sk2B + np.ones((4, 1)) * z * Sk2C) * np.cos(delta) -
                (Sk3A - Sk3A_ + Sk3B + np.ones((4, 1)) * z * Sk3C) * np.sin(delta))
        Sduzdy = 1 / (2 * np.pi) * slip_str * (
                (Sk2A - Sk2A_ + Sk2B - np.ones((4, 1)) * z * Sk2C) * np.sin(delta) +
                (Sk3A - Sk3A_ + Sk3B - np.ones((4, 1)) * z * Sk3C) * np.cos(delta))

        Sduxdz = 1 / (2 * np.pi) * slip_str * (Sl1A + Sl1A_ + Sl1B + Su1C + np.ones((4, 1)) * z * Sl1C)
        Sduydz = 1 / (2 * np.pi) * slip_str * (
                (Sl2A + Sl2A_ + Sl2B + Su2C + np.ones((4, 1)) * z * Sl2C) * np.cos(delta) -
                (Sl3A + Sl3A_ + Sl3B + Su3C + np.ones((4, 1)) * z * Sl3C) * np.sin(delta))
        Sduzdz = 1 / (2 * np.pi) * slip_str * (
                (Sl2A + Sl2A_ + Sl2B - Su2C - np.ones((4, 1)) * z * Sl2C) * np.sin(delta) +
                (Sl3A + Sl3A_ + Sl3B - Su3C - np.ones((4, 1)) * z * Sl3C) * np.cos(delta))
    else:
        Sux, Suy, Suz = 0, 0, 0
        Sduxdx, Sduydx, Sduzdx = 0, 0, 0
        Sduxdy, Sduydy, Sduzdy = 0, 0, 0
        Sduxdz, Sduydz, Sduzdz = 0, 0, 0

    # dip-slip
    if slip_dip != 0:
        # displacement
        #uA
        Du1A = alpha / 2 * q / R
        Du2A = theta / 2 + alpha / 2 * eta * q * X11
        Du3A = (1 - alpha) / 2 * np.log(R + xi) - alpha / 2 * q**2 * X11
        #uB
        Du1B = -q / R + (1 - alpha) / alpha * I3 * np.sin(delta) * np.cos(delta)
        Du2B = -eta * q * X11 - theta - (1 - alpha) / alpha * xi / (R + d_) * np.sin(delta) * np.cos(delta)
        Du3B = q**2 * X11 + (1 - alpha) / alpha * I4 * np.sin(delta) * np.cos(delta)
        #uC
        Du1C = (1 - alpha) * np.cos(delta) / R - q * Y11 * np.sin(delta) - alpha * c_ * q / R**3
        Du2C = (1 - alpha) * y_ * X11 - alpha * c_ * eta * q * X32
        Du3C = -d_ * X11 - xi * Y11 * np.sin(delta) - alpha * c_ * (X11 - q**2 * X32)
        # displacement gradient
        #jA
        Dj1A = -alpha / 2 * xi * q / R**3
        Dj2A = -q / 2 * Y11 - alpha / 2 * eta * q / R**3
        Dj3A = (1 - alpha) / 2 * 1 / R + alpha / 2 * q**2 / R**3
        #jB
        Dj1B = xi * q / R**3 + (1 - alpha) / alpha * J4 * np.sin(delta) * np.cos(delta)
        Dj2B = eta * q / R**3 + q * Y11 + (1 - alpha) / alpha * J5 * np.sin(delta) * np.cos(delta)
        Dj3B = -q**2 / R**3 + (1 - alpha) / alpha * J6 * np.sin(delta) * np.cos(delta)
        #jC
        Dj1C = -(1 - alpha) * xi / R**3 * np.cos(delta) + xi * q * Y32 * np.sin(delta) + alpha * (
                        3 * c_ * xi * q / R**5)
        Dj2C = -(1 - alpha) * y_ / R**3 + alpha * 3 * c_ * eta * q / R**5
        Dj3C = d_ / R**3 - Y0 * np.sin(delta) + alpha * c_ / R**3 * (1 - 3 * q**2 / R**2)
        #kA
        Dk1A = alpha / 2 * E
        Dk2A = (1 - alpha) / 2 * d_ * X11 + xi / 2 * Y11 * np.sin(delta) + alpha / 2 * eta * G
        Dk3A = (1 - alpha) / 2 * y_ * X11 - alpha / 2 * q * G
        #kB
        Dk1B = -E + (1 - alpha) / alpha * J1 * np.sin(delta) * np.cos(delta)
        Dk2B = -eta * G - xi * Y11 * np.sin(delta) + (1 - alpha) / alpha * J2 * np.sin(delta) * np.cos(delta)
        Dk3B = q * G + (1 - alpha) / alpha * J3 * np.sin(delta) * np.cos(delta)
        #kC
        Dk1C = -(1 - alpha) * eta / R**3 + Y0 * np.sin(delta)**2 - alpha * (
                        (c_ + d_) / R**3 * np.sin(delta) - 3 * c_ * y_ * q / R**5)
        Dk2C = (1 - alpha) * (X11 - y_**2 * X32) - alpha * c_ * (
                        (d_ + 2 * q * np.cos(delta)) * X32 - y_ * eta * q * X53)
        Dk3C = xi * P * np.sin(delta) + y_ * d_ * X32 + alpha * c_ * (
                        (y_ + 2 * q * np.sin(delta)) * X32 - y_ * q**2 * X53)
        #lA
        Dl1A = alpha / 2 * E_
        Dl2A = (1 - alpha) / 2 * y_ * X11 + xi / 2 * Y11 * np.cos(delta) + alpha / 2 * eta * G_
        Dl3A = -(1 - alpha) / 2 * d_ * X11 - alpha / 2 * q * G_
        #lB
        Dl1B = -E_ - (1 - alpha) / alpha * K3 * np.sin(delta) * np.cos(delta)
        Dl2B = -eta * G_ - xi * Y11 * np.cos(delta) - (1 - alpha) / alpha * xi * D11 * np.sin(delta) * np.cos(delta)
        Dl3B = q * G_ - (1 - alpha) / alpha * K4 * np.sin(delta) * np.cos(delta)
        #lB
        Dl1C = -q / R**3 + Y0 * np.sin(delta) * np.cos(delta) - alpha * (
                        (c_ + d_) / R**3 * np.cos(delta) + 3 * c_ * d_ * q / R**5)
        Dl2C = (1 - alpha) * y_ * d_ * X32 - alpha * c_ * (
                        (y_ - 2 * q * np.sin(delta)) * X32 + d_ * eta * q * X53)
        Dl3C = -xi * P_ * np.sin(delta) + X11 - d_**2 * X32 - alpha * c_ * (
                        (d_ - 2 * q * np.cos(delta)) * X32 - d_ * q**2 * X53)
        # displacement
        # u1A_
        Du1A_ = alpha / 2 * q / R
        Du2A_ = theta / 2 + alpha / 2 * eta * q * X11
        Du3A_ = (1 - alpha) / 2 * np.log(R + xi) - alpha / 2 * q**2 * X11
        # displacement gradient
        # jA_
        Dj1A_ = - alpha / 2 * xi * q / R**3
        Dj2A_ = - q / 2 * Y11 - alpha / 2 * eta * q / R**3
        Dj3A_ = (1 - alpha) / 2 * 1 / R + alpha / 2 * q**2 / R**3
        # kA_
        Dk1A_ = alpha / 2 * E
        Dk2A_ = (1 - alpha) / 2 * d_ * X11 + xi / 2 * Y11 * np.sin(delta) + alpha / 2 * eta * G
        Dk3A_ = (1 - alpha) / 2 * y_ * X11 - alpha / 2 * q * G
        # lA_
        Dl1A_ = alpha / 2 * E_
        Dl2A_ = (1 - alpha) / 2 * y_ * X11 + xi / 2 * Y11 * np.cos(delta) + alpha / 2 * eta * G_
        Dl3A_ = - (1 - alpha) / 2 * d_ * X11 - alpha / 2 * q * G_


        # displacement
        Dux = 1 / (2 * np.pi) * slip_dip * (Du1A - Du1A_ + Du1B + np.ones((4, 1)) * z * Du1C)
        Duy = 1 / (2 * np.pi) * slip_dip * (
                (Du2A - Du2A_ + Du2B + np.ones((4, 1)) * z * Du2C) * np.cos(delta) -
                (Du3A - Du3A_ + Du3B + np.ones((4, 1)) * z * Du3C) * np.sin(delta))
        Duz = 1 / (2 * np.pi) * slip_dip * (
                (Du2A - Du2A_ + Du2B - np.ones((4, 1)) * z * Du2C) * np.sin(delta) +
                (Du3A - Du3A_ + Du3B - np.ones((4, 1)) * z * Du3C) * np.cos(delta))

        # displacement gradients
        Dduxdx = 1 / (2 * np.pi) * slip_dip * (Dj1A - Dj1A_ + Dj1B + np.ones((4, 1)) * z * Dj1C)
        Dduydx = 1 / (2 * np.pi) * slip_dip * (
                (Dj2A - Dj2A_ + Dj2B + np.ones((4, 1)) * z * Dj2C) * np.cos(delta) -
                (Dj3A - Dj3A_ + Dj3B + np.ones((4, 1)) * z * Dj3C) * np.sin(delta))
        Dduzdx = 1 / (2 * np.pi) * slip_dip * (
                (Dj2A - Dj2A_ + Dj2B - np.ones((4, 1)) * z * Dj2C) * np.sin(delta) +
                (Dj3A - Dj3A_ + Dj3B - np.ones((4, 1)) * z * Dj3C) * np.cos(delta))
        
        Dduxdy = 1 / (2 * np.pi) * slip_dip * (Dk1A - Dk1A_ + Dk1B + np.ones((4, 1)) * z * Dk1C)
        Dduydy = 1 / (2 * np.pi) * slip_dip * (
                (Dk2A - Dk2A_ + Dk2B + np.ones((4, 1)) * z * Dk2C) * np.cos(delta) -
                (Dk3A - Dk3A_ + Dk3B + np.ones((4, 1)) * z * Dk3C) * np.sin(delta))
        Dduzdy = 1 / (2 * np.pi) * slip_dip * (
                (Dk2A - Dk2A_ + Dk2B - np.ones((4, 1)) * z * Dk2C) * np.sin(delta) +
                (Dk3A - Dk3A_ + Dk3B - np.ones((4, 1)) * z * Dk3C) * np.cos(delta))

        Dduxdz = 1 / (2 * np.pi) * slip_dip * (Dl1A + Dl1A_ + Dl1B + Du1C + np.ones((4, 1)) * z * Dl1C)
        Dduydz = 1 / (2 * np.pi) * slip_dip * (
                (Dl2A + Dl2A_ + Dl2B + Du2C + np.ones((4, 1)) * z * Dl2C) * np.cos(delta) -
                (Dl3A + Dl3A_ + Dl3B + Du3C + np.ones((4, 1)) * z * Dl3C) * np.sin(delta))
        Dduzdz = 1 / (2 * np.pi) * slip_dip * (
                (Dl2A + Dl2A_ + Dl2B - Du2C - np.ones((4, 1)) * z * Dl2C) * np.sin(delta) +
                (Dl3A + Dl3A_ + Dl3B - Du3C - np.ones((4, 1)) * z * Dl3C) * np.cos(delta))
        
    else:
        Dux, Duy, Duz = 0, 0, 0
        Dduxdx, Dduydx, Dduzdx = 0, 0, 0
        Dduxdy, Dduydy, Dduzdy = 0, 0, 0
        Dduxdz, Dduydz, Dduzdz = 0, 0, 0


    if tensile != 0:
        # displacement
        # uA
        Tu1A = -(1 - alpha) / 2 * np.log(R + eta) - alpha / 2 * q**2 * Y11
        Tu2A = -(1 - alpha) / 2 * np.log(R + xi) - alpha / 2 * q**2 * X11
        Tu3A = theta / 2 - alpha / 2 * q * (eta * X11 + xi * Y11)
        # uB
        Tu1B = q**2 * Y11 - (1 - alpha) / alpha * I3 * np.sin(delta)**2
        Tu2B = q**2 * X11 + (1 - alpha) / alpha * xi / (R + d_) * np.sin(delta)**2
        Tu3B = q * (eta * X11 + xi * Y11) - theta - (1 - alpha) / alpha * I4 * np.sin(delta)**2
        # uC
        Tu1C = -(1 - alpha) * (np.sin(delta) / R + q * Y11 * np.cos(delta)) - alpha * (
                np.ones((4, 1)) * z * Y11 - q**2 * Z32)
        Tu2C = (1 - alpha) * 2 * xi * Y11 * np.sin(delta) + d_ * X11 - alpha * c_ * (
                X11 - q**2 * X32)
        Tu3C = (1 - alpha) * (y_ * X11 + xi * Y11 * np.cos(delta)) + alpha * q * (
                c_ * eta * X32 + xi * Z32)
        # displacement gradient
        # jA
        Tj1A = - (1 - alpha) / 2 * xi * Y11 + alpha / 2 * xi * q**2 * Y32
        Tj2A = - (1 - alpha) / 2 * 1 / R + alpha / 2 * q**2 / R**3
        Tj3A = - (1 - alpha) / 2 * q * Y11 - alpha / 2 * q**3 * Y32
        # jB
        Tj1B = -xi * q**2 * Y32 - (1 - alpha) / alpha * J4 * np.sin(delta)**2
        Tj2B = -q**2 / R**3 - (1 - alpha) / alpha * J5 * np.sin(delta)**2
        Tj3B = q**3 * Y32 - (1 - alpha) / alpha * J6 * np.sin(delta)**2
        # jC
        Tj1C = (1 - alpha) * xi / R**3 * np.sin(delta) + xi * q * Y32 * np.cos(delta) + alpha * xi * (
                3 * c_ * eta / R**5 - 2 * Z32 - Z0)
        Tj2C = (1 - alpha) * 2 * Y0 * np.sin(delta) - d_ / R**3 + alpha * c_ / R**3 * (
                1 - 3 * q**2 / R**2)
        Tj3C = -(1 - alpha) * (y_ / R**3 - Y0 * np.cos(delta)) - alpha * (
                3 * c_ * eta * q / R**5 - q * Z0)
        # kA
        Tk1A = -(1 - alpha) / 2 * (np.cos(delta) / R + q * Y11 * np.sin(delta)) - alpha / 2 * q * F
        Tk2A = -(1 - alpha) / 2 * y_ * X11 - alpha / 2 * q * G
        Tk3A = (1 - alpha) / 2 * (d_ * X11 + xi * Y11 * np.sin(delta)) + alpha / 2 * q * H
        # kB
        Tk1B = q * F - (1 - alpha) / alpha * J1 * np.sin(delta)**2
        Tk2B = q * G - (1 - alpha) / alpha * J2 * np.sin(delta)**2
        Tk3B = -q * H - (1 - alpha) / alpha * J3 * np.sin(delta)**2
        # kC
        Tk1C = (1 - alpha) * (q / R**3 + Y0 * np.sin(delta) * np.cos(delta)) + alpha * (
                np.ones((4, 1)) * z / R**3 * np.cos(delta) + 3 * c_ * d_ * q / R**5 - q * Z0 * np.sin(delta))
        Tk2C = -(1 - alpha) * 2 * xi * P * np.sin(delta) - y_ * d_ * X32 + alpha * c_ * (
                (y_ + 2 * q * np.sin(delta)) * X32 - y_ * q**2 * X53)
        Tk3C = -(1 - alpha) * (xi * P * np.cos(delta) - X11 + y_**2 * X32) + alpha * c_ * (
                (d_ + 2 * q * np.cos(delta)) * X32 - y_ * eta * q * X53) + alpha * xi * Q
        # lA
        Tl1A = (1 - alpha) / 2 * (np.sin(delta) / R - q * Y11 * np.cos(delta)) - alpha / 2 * q * F_
        Tl2A = (1 - alpha) / 2 * d_ * X11 - alpha / 2 * q * G_
        Tl3A = (1 - alpha) / 2 * (y_ * X11 + xi * Y11 * np.cos(delta)) + alpha / 2 * q * H_
        # lB
        Tl1B = q * F_ + (1 - alpha) / alpha * K3 * np.sin(delta)**2
        Tl2B = q * G_ + (1 - alpha) / alpha * xi * D11 * np.sin(delta)**2
        Tl3B = -q * H_ + (1 - alpha) / alpha * K4 * np.sin(delta)**2
        # lC
        Tl1C = -eta / R**3 + Y0 * np.cos(delta)**2 - alpha * (
                np.ones((4, 1)) * z / R**3 * np.sin(delta) - 3 * c_ * y_ * q / R**5 - Y0 * np.sin(delta)**2 + q * Z0 * np.cos(delta))
        Tl2C = (1 - alpha) * 2 * xi * P_ * np.sin(delta) - X11 + d_**2 * X32 - alpha * c_ * (
                (d_ - 2 * q * np.cos(delta)) * X32 - d_ * q**2 * X53)
        Tl3C = (1 - alpha) * (
                xi * P_ * np.cos(delta) + y_ * d_ * X32) + alpha * c_ * (
                (y_ - 2 * q * np.sin(delta)) * X32 + d_ * eta * q * X53) + alpha * xi * Q_
        # displacement
        # uA_
        Tu1A_ = - (1 - alpha) / 2 * np.log(R + eta) - alpha / 2 * q**2 * Y11
        Tu2A_ = - (1 - alpha) / 2 * np.log(R + xi) - alpha / 2 * q**2 * X11
        Tu3A_ = theta / 2 - alpha / 2 * q * (eta * X11 + xi * Y11)
        # displacement gradient
        # jA_
        Tj1A_ = - (1 - alpha) / 2 * xi * Y11 + alpha / 2 * xi * q**2 * Y32
        Tj2A_ = - (1 - alpha) / 2 * 1 / R + alpha / 2 * q**2 / R**3
        Tj3A_ = - (1 - alpha) / 2 * q * Y11 - alpha / 2 * q**3 * Y32
        # kA_
        Tk1A_ = - (1 - alpha) / 2 * (np.cos(delta) / R + q * Y11 * np.sin(delta)) - alpha / 2 * q * F
        Tk2A_ = - (1 - alpha) / 2 * y_ * X11 - alpha / 2 * q * G
        Tk3A_ = (1 - alpha) / 2 * (d_ * X11 + xi * Y11 * np.sin(delta)) + alpha / 2 * q * H
        # lA_
        Tl1A_ = (1 - alpha) / 2 * (np.sin(delta) / R - q * Y11 * np.cos(delta)) - alpha / 2 * q * F_
        Tl2A_ = (1 - alpha) / 2 * d_ * X11 - alpha / 2 * q * G_
        Tl3A_ = (1 - alpha) / 2 * (y_ * X11 + xi * Y11 * np.cos(delta)) + alpha / 2 * q * H_
    
        # Tensile
        # Displacement
        Tux = (1 / (2 * np.pi)) * tensile * (Tu1A - Tu1A_ + Tu1B + np.ones((4, 1)) * z * Tu1C)
        Tuy = (1 / (2 * np.pi)) * tensile * ((Tu2A - Tu2A_ + Tu2B + np.ones((4, 1)) * z * Tu2C) * np.cos(delta) -
                                                (Tu3A - Tu3A_ + Tu3B + np.ones((4, 1)) * z * Tu3C) * np.sin(delta))
        Tuz = (1 / (2 * np.pi)) * tensile * ((Tu2A - Tu2A_ + Tu2B - np.ones((4, 1)) * z * Tu2C) * np.sin(delta) +
                                                (Tu3A - Tu3A_ + Tu3B - np.ones((4, 1)) * z * Tu3C) * np.cos(delta))

        # Displacement gradients
        Tduxdx = (1 / (2 * np.pi)) * tensile * (Tj1A - Tj1A_ + Tj1B + np.ones((4, 1)) * z * Tj1C)
        Tduydx = (1 / (2 * np.pi)) * tensile * ((Tj2A - Tj2A_ + Tj2B + np.ones((4, 1)) * z * Tj2C) * np.cos(delta) -
                                                (Tj3A - Tj3A_ + Tj3B + np.ones((4, 1)) * z * Tj3C) * np.sin(delta))
        Tduzdx = (1 / (2 * np.pi)) * tensile * ((Tj2A - Tj2A_ + Tj2B - np.ones((4, 1)) * z * Tj2C) * np.sin(delta) +
                                                (Tj3A - Tj3A_ + Tj3B - np.ones((4, 1)) * z * Tj3C) * np.cos(delta))
        Tduxdy = (1 / (2 * np.pi)) * tensile * (Tk1A - Tk1A_ + Tk1B + np.ones((4, 1)) * z * Tk1C)
        Tduydy = (1 / (2 * np.pi)) * tensile * ((Tk2A - Tk2A_ + Tk2B + np.ones((4, 1)) * z * Tk2C) * np.cos(delta) -
                                                (Tk3A - Tk3A_ + Tk3B + np.ones((4, 1)) * z * Tk3C) * np.sin(delta))
        Tduzdy = (1 / (2 * np.pi)) * tensile * ((Tk2A - Tk2A_ + Tk2B - np.ones((4, 1)) * z * Tk2C) * np.sin(delta) +
                                                (Tk3A - Tk3A_ + Tk3B - np.ones((4, 1)) * z * Tk3C) * np.cos(delta))
        Tduxdz = (1 / (2 * np.pi)) * tensile * (Tl1A + Tl1A_ + Tl1B + Tu1C + np.ones((4, 1)) * z * Tl1C)
        Tduydz = (1 / (2 * np.pi)) * tensile * ((Tl2A + Tl2A_ + Tl2B + Tu2C + np.ones((4, 1)) * z * Tl2C) * np.cos(delta) -
                                                (Tl3A + Tl3A_ + Tl3B + Tu3C + np.ones((4, 1)) * z * Tl3C) * np.sin(delta))
        Tduzdz = (1 / (2 * np.pi)) * tensile * ((Tl2A + Tl2A_ + Tl2B - Tu2C - np.ones((4, 1)) * z * Tl2C) * np.sin(delta) +
                                                (Tl3A + Tl3A_ + Tl3B - Tu3C - np.ones((4, 1)) * z * Tl3C) * np.cos(delta))
    else:
        Tux, Tuy, Tuz = 0, 0, 0
        Tduxdx, Tduydx, Tduzdx = 0, 0, 0
        Tduxdy, Tduydy, Tduzdy = 0, 0, 0
        Tduxdz, Tduydz, Tduzdz = 0, 0, 0

    factor=np.ones((xi.shape))
    factor[1,:]=factor[1,:]*-1
    factor[2,:]=factor[2,:]*-1

    G1 = np.sum(factor * (Sux + Dux + Tux),axis=0)
    G2 = np.sum(factor * (Suy + Duy + Tuy),axis=0)
    G3 = np.sum(factor * (Suz + Duz + Tuz),axis=0)

    Dg11 = np.sum(factor * (Sduxdx + Dduxdx + Tduxdx),axis=0)
    Dg12 = np.sum(factor * (Sduxdy + Dduxdy + Tduxdy),axis=0)
    Dg13 = np.sum(factor * (Sduxdz + Dduxdz + Tduxdz),axis=0)

    Dg21 = np.sum(factor * (Sduydx + Dduydx + Tduydx),axis=0)
    Dg22 = np.sum(factor * (Sduydy + Dduydy + Tduydy),axis=0)
    Dg23 = np.sum(factor * (Sduydz + Dduydz + Tduydz),axis=0)

    Dg31 = np.sum(factor * (Sduzdx + Dduzdx + Tduzdx),axis=0)
    Dg32 = np.sum(factor * (Sduzdy + Dduzdy + Tduzdy),axis=0)
    Dg33 = np.sum(factor * (Sduzdz + Dduzdz + Tduzdz),axis=0)

    # Coordinate transformation
    Gx = np.cos(angle_Str) * (-G2) - np.sin(angle_Str) * G1
    Gy = np.sin(angle_Str) * (-G2) + np.cos(angle_Str) * G1
    Gz = G3

    displacement = np.array([Gx, Gy, Gz])

    Dg11_ = Dg22
    Dg12_ = -Dg21
    Dg13_ = -Dg23

    Dg21_ = -Dg12
    Dg22_ = Dg11
    Dg23_ = Dg13

    Dg31_ = -Dg32
    Dg32_ = Dg31
    Dg33_ = Dg33

    # Coordinate transformation
    Dgxx = (np.cos(angle_Str) * Dg11_ - np.sin(angle_Str) * Dg21_) * np.cos(angle_Str) + \
    (np.cos(angle_Str) * Dg12_ - np.sin(angle_Str) * Dg22_) * (-np.sin(angle_Str))
    Dgyx = (np.sin(angle_Str) * Dg11_ + np.cos(angle_Str) * Dg21_) * np.cos(angle_Str) - \
    (np.sin(angle_Str) * Dg12_ + np.cos(angle_Str) * Dg22_) * np.sin(angle_Str)
    Dgzx = Dg31_ * np.cos(angle_Str) - Dg32_ * np.sin(angle_Str)

    Dgxy = (np.cos(angle_Str) * Dg11_ - np.sin(angle_Str) * Dg21_) * np.sin(angle_Str) + \
    (np.cos(angle_Str) * Dg12_ - np.sin(angle_Str) * Dg22_) * np.cos(angle_Str)
    Dgyy = (np.sin(angle_Str) * Dg11_ + np.cos(angle_Str) * Dg21_) * np.sin(angle_Str) + \
    (np.sin(angle_Str) * Dg12_ + np.cos(angle_Str) * Dg22_) * np.cos(angle_Str)
    Dgzy = np.sin(angle_Str) * Dg31_ + np.cos(angle_Str) * Dg32_

    Dgxz = np.cos(angle_Str) * Dg13_ - np.sin(angle_Str) * Dg23_
    Dgyz = np.sin(angle_Str) * Dg13_ + np.cos(angle_Str) * Dg23_
    Dgzz = Dg33_

    gradient = np.array([[Dgxx, Dgxy, Dgxz],
                    [Dgyx, Dgyy, Dgyz],
                    [Dgzx, Dgzy, Dgzz]])
    
    gradient = gradient.reshape(9, x.size)
    # Strain components
    Ex = Dgxx
    Ey = Dgyy
    Ez = Dgzz
    Exy = 0.5 * (Dgyx + Dgxy)
    Eyz = 0.5 * (Dgyz + Dgzy)
    Ezx = 0.5 * (Dgzx + Dgxz)

    # Stress components
    Sx = youngs / ((1 + nu) * (1 - 2 * nu)) * (Ex + nu * (Ey + Ez - Ex))
    Sy = youngs / ((1 + nu) * (1 - 2 * nu)) * (Ey + nu * (Ex + Ez - Ey))
    Sz = youngs / ((1 + nu) * (1 - 2 * nu)) * (Ez + nu * (Ey + Ex - Ez))
    Sxy = 2 * Gshear * Exy
    Syz = 2 * Gshear * Eyz
    Szx = 2 * Gshear * Ezx

    # Stress vector
    Stress = np.array([Sx, Sxy, Szx, Sy, Syz, Sz])

    return(displacement,gradient,Stress)

def disloc3D_wrapper(m,coordinates,return_stress=False,return_2d=False):

        class makeFault:

                nu = params.nu
                Gshear = params.Gshear

                def __init__(self, m):
                        self.L = m[0]
                        self.W = m[1]
                        self.c = m[2]
                        self.delta = np.deg2rad(m[3])
                        self.strike = -np.deg2rad(m[4]) # clockwise is positive
                        self.centroids = np.array([m[5],m[6]])
                        self.slip = m[7], m[8], m[9]
                        
                        self.rotate = np.array([[-np.sin(self.strike), np.cos(self.strike)],
                                                [-np.cos(self.strike), -np.sin(self.strike)]])
                def coordTrans(self,X,Y):
                        coordinates_vector = np.array([X-self.centroids[0],Y-self.centroids[1]])
                        x, y = np.dot(self.rotate, coordinates_vector)
                        x = x + 0.5 * self.L
                        return(x,y)
                
                @classmethod
                def elastic(cls):
                        if cls.nu == 0.5: cls.nu = 0.4999
                        youngs = 2 * cls.Gshear * (1 + cls.nu)
                        lambda_ = cls.nu * youngs / ((1 + cls.nu) * (1 - 2 * cls.nu))
                        alpha = (lambda_ + cls.Gshear) / (lambda_ + 2 * cls.Gshear)
                        return(youngs,alpha)
                
        disloc=makeFault(np.squeeze(m))

        if checkInputs(disloc,coordinates) == True: 
              displacement,gradient,Stress=disloc3D(disloc,coordinates)
        else: 
              displacement, gradient, Stress = np.nan


        if return_2d:
            # in this wrapper function, ensure efficient memory use when called by only returning neccesary components
            # indicies for gradient 2D correspond to xx, xy, yx, yy
            # indicies for stress correspond to xx, xy, yy (since xy = yx)
            displacement_2D=displacement[[0,1],:]
            gradient_2D=gradient[[0,1,3,4],:]
            stress_2D=Stress[[0,1,3],:]
            if return_stress:
                 return(displacement_2D,gradient_2D,stress_2D)
            else:
                 return(displacement_2D,gradient_2D)
        else:        
            if return_stress:
                return(displacement,gradient,Stress)
            else:
                return(displacement,gradient)

# cholesky factorization much faster than numpy's standard multivariate_normal
# which uses svd for legacy reasons (?)
# still much slower than matlab's mvnrnd, but improved
# currently not using.
def chol_mvnrnd(mean, cov, num):
    result = np.empty((num,mean.size))
    for i in range(num):
        result[i,:] = mean + np.linalg.cholesky(cov) @ np.random.standard_normal(mean.size)
    return(result)