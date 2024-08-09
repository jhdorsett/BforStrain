import numpy as np
import tools
import params
import scipy
import pickle
import pyvista as pv
import matplotlib.path as MPLpoly
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from datetime import datetime
import time

class bounding_box:

    def __init__(self,range=False):
        # by default, will create bounding box in local CRS using 
        # the lon/lat ranges specified in input file
        # specify range as range[:,0]=minx,maxx
        #                  range[:,1]=miny,maxy
        if not range:
            [self.minx, self.miny] = tools.llh2local(params.lon_range[0],params.lat_range[0])
            [self.maxx, self.maxy] = tools.llh2local(params.lon_range[1],params.lat_range[1])

        else:
            [self.minx, self.miny] = range[0,:]
            [self.maxx, self.maxy] = range[1,:]

        if (self.maxx < self.minx) | (self.maxy < self.miny):
            print("Warning: Check coordinate ranges for mesh domain in params.py file.")
            print("Attempting to fix.")
            if self.maxx < self.minx: self.maxx, self.minx = self.minx, self.maxx
            if self.maxy < self.miny: self.maxy, self.miny = self.miny, self.maxy

        self.nx=round((self.maxx-self.minx)/params.nom_node_spacing)+1
        self.ny=round((self.maxy-self.miny)/params.nom_node_spacing)+1

    def make_grid(self):
        X,Y=np.meshgrid(np.linspace(self.minx,self.maxx,self.nx),np.linspace(self.miny,self.maxy,self.ny))
        nodes = np.column_stack((X.ravel(), Y.ravel()))
        return(nodes)
    
    def crop(self,vals,index=False):
        within_bounds = (
            (vals[:, 0] > self.minx) & 
            (vals[:, 0] < self.maxx) & 
            (vals[:, 1] > self.miny) & 
            (vals[:, 1] < self.maxy)
        )
        subset = vals[within_bounds]
        if index:
            return(subset,within_bounds)
        else:
            return(subset)
                
class Velocity:
    def __init__(self, Ve,Vn,Sige,Sign):
        self.Ve = Ve
        self.Vn = Vn
        self.Sige = Sige
        self.Sign = Sign

# since pickle files just return objects, these are the same function

def load(filename):
    if filename[-4:] != (".pkl"):
        print("Mesh must be a .pkl file. Use Mesh.save_mesh, not Mesh.save_mesh_vtk.")

    pkl_file = open(filename,'rb')
    return(pickle.load(pkl_file))

class Mesh:
    def __init__(self):

        data = np.loadtxt(params.gps_file,delimiter=',')
        data = data[::params.decimate,:]
        self.lon = data[:,0]
        self.lat = data[:,1]

        ind = self.lon > 180 
        if sum(ind) > 0:
            self.lon[ind] = self.lon[ind]-360
            print("Longitude data appears to span 0-360, applying correction.")

        self.xy_gps = tools.llh2local(self.lon,self.lat) 

        self.vel = Velocity(data[:, 2],data[:, 3],data[:, 4],data[:, 5])
        # set minimum velocity sigma
        self.vel.Sige[self.vel.Sige < 0.5] = 0.5
        self.vel.Sign[self.vel.Sign < 0.5] = 0.5

        if bool(params.creeping_file):
            self.creeping = True
            self.creeping_faults = np.loadtxt(params.creeping_file,delimiter=',')
        else:
            self.creeping = False

        print("Building mesh from coordinates and extent given in params.py")
        self.bounds = bounding_box()

    def construct(self):
        tools.build_mesh(self)

        if self.creeping:
            self.getPM()

        self.assign_centroids()

    def assign_centroids(self):
        
        elts=self.nodes[self.tri]
        self.tri_centroids = np.mean(elts, axis=1)
        data_tri=np.zeros(self.xy_gps.shape[0])

        for k, tri in enumerate(elts):
            path=MPLpoly.Path(tri)
            inside=path.contains_points(self.xy_gps)
            data_tri[inside]=k

        Ve_centroids = np.full(self.tri.shape[0], np.nan)
        Vn_centroids = np.full(self.tri.shape[0], np.nan)
        Sige_centroids = np.full(self.tri.shape[0], np.nan)
        Sign_centroids = np.full(self.tri.shape[0], np.nan)

        for k in range(self.tri.shape[0]):
            ind = (data_tri == k )  
            if np.sum(ind) > 0:
                W_sige = 1.0 / self.vel.Sige[ind]
                W_sign = 1.0 / self.vel.Sign[ind]
                
                Ve_centroids[k] = np.sum(W_sige * self.vel.Ve[ind]) / np.sum(W_sige)
                Vn_centroids[k] = np.sum(W_sign * self.vel.Vn[ind]) / np.sum(W_sign)

                Sige_centroids[k] = np.sqrt(np.sum(self.vel.Sige[ind]**2))
                Sign_centroids[k] = np.sqrt(np.sum(self.vel.Sign[ind]**2))

        self.vel_centroids = Velocity(Ve_centroids,Vn_centroids,Sige_centroids,Sign_centroids)
        self.ind = np.squeeze(~np.isnan(Ve_centroids))        

    def getPM(self):
        # Calculate lengths and angles
        delta_x = self.PatchEnds[:, 2] - self.PatchEnds[:, 0]
        delta_y = self.PatchEnds[:, 3] - self.PatchEnds[:, 1]
        lengths = np.sqrt(delta_x**2 + delta_y**2)
        angles = np.arctan2(delta_y, delta_x)
        strike = 90 - angles * 180 / np.pi

        # Calculate centers
        center_x = (self.PatchEnds[:, 0] + self.PatchEnds[:, 2]) / 2
        center_y = (self.PatchEnds[:, 1] + self.PatchEnds[:, 3]) / 2

        # Create pm array
        self.pm = np.stack((
            lengths,
            np.full(lengths.shape, 1e6),
            np.full(lengths.shape, 1e6),
            np.full(lengths.shape, 90),
            strike,
            center_x,
            center_y
        ), axis=-1)

    def plot_mesh(self,score=False):

        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        colormap=plt.get_cmap('viridis')
        if isinstance(score, bool):
            colors = colormap(self.score)
        else:
            min_score = score.min()
            max_score = score.max()
            normalized_score = (score - min_score) / (max_score - min_score)
            colors=colormap(normalized_score)

        # Loop through each triangle and plot
        for tri,color in zip(self.tri,colors):
            polygon = Polygon(self.nodes[tri], closed=True, edgecolor='k', facecolor=color[:3], alpha=color[0])
            ax.add_patch(polygon)

        ax.set_xlabel('km')
        x_min, x_max = self.nodes[:, 0].min(), self.nodes[:, 0].max()
        y_min, y_max = self.nodes[:, 1].min(), self.nodes[:, 1].max()
        
        ax.set_xlim(x_min , x_max )
        ax.set_ylim(y_min , y_max )
        ax.set_title('Smooth2 score metric')
        plt.show()

    def plot_vels(self):

        fig1, ax1=plt.subplots()
        V_centroids = np.sqrt(self.vel_centroids.Ve**2+self.vel_centroids.Vn**2)
        ax1.tripcolor(self.nodes[:,0], self.nodes[:,1],self.tri,linewidth=0.25,facecolors=V_centroids,edgecolors='k')
        ax1.quiver(self.xy_gps[:,0],self.xy_gps[:,1],self.vel.Ve,self.vel.Vn, color='red')
        fig1.gca().set_aspect('equal')
        ax1.set_title('velocity data')
        # Make the plot interactive
        
        # Keep the script running to allow viewing of the plots
        plt.show()
    
    def plot_results(self,X,Y):
        XY = np.sqrt(X**2+Y**2)
        #plt.tripcolor(nodes[:,0], nodes[:,1],triDel.simplices.copy(),linewidth=0.25,facecolors=areas)
        fig=plt.tripcolor(self.nodes[:,0], self.nodes[:,1],self.tri,linewidth=0.25,facecolors=XY,edgecolors='k')
        plt.figsize=(8,8)
        plt.quiver(self.xy_gps[:,0],self.xy_gps[:,1],self.vel.Ve,self.vel.Vn, color='red')
        #plt.plot(PatchEnds[:, [0, 2]].T, PatchEnds[:, [1, 3]].T, color='blue')  # Plot all line segments at once
        plt.gca().set_aspect('equal')
        plt.title('modeled velocity')
        # Make the plot interactive
        plt.colorbar()
        plt.show()

    def save(self,filename=False):
        import pickle
        # if file unspecified, generate default name
        if not(filename):
            timestamp = datetime.now().strftime("%d-%m-%Y_%H:%M")
            filename = f"{timestamp}_mesh.pkl"
        output = open(filename,'wb')
        pickle.dump(self,output)
        output.close()

    def save_vtk(self,results=False,filename=False,save_data=False):

        # if file unspecified, generate default name
        if not(filename):
            timestamp = datetime.now().strftime("%d-%m-%Y_%H:%M")
            filename = f"{timestamp}.vtk"

        if save_data:
            # first create file with observed velocities. this is an unstructured point cloud
            gps_loc = np.hstack((self.xy_gps,np.zeros((self.xy_gps.shape[0],1))))
            gps_obs = pv.PolyData(gps_loc)

            gps_obs.point_data['Observed Velocity']  = np.column_stack((self.vel.Ve,self.vel.Vn))
            gps_obs.point_data['Velocity Uncertainty']  = np.column_stack((self.vel.Sige,self.vel.Sign))

            gps_obs.save(f"{filename[:-4]}_data.vtk",binary = False)

        # format nodes and triangles for pyvista
        mesh_points = np.hstack((self.nodes,np.zeros((self.nodes.shape[0],1))))
        cell_types = np.full(len(self.tri), pv.CellType.TRIANGLE, dtype=np.int8)
        tri_idx=np.full((self.tri.shape[0],1),3)
        connectivity = np.hstack((tri_idx,self.tri))
        connectivity = connectivity.ravel()

        # Create a PyVista mesh object and populate with velocity observations
        grid = pv.UnstructuredGrid(connectivity, cell_types, mesh_points)

        grid.cell_data['Observed Velocity']  = np.column_stack((self.vel_centroids.Ve,self.vel_centroids.Vn))
        grid.cell_data['Velocity Uncertainty']  = np.column_stack((self.vel_centroids.Sige,self.vel_centroids.Sign))
        
        grid.lines = np.hstack([[2, *edge] for edge in self.edges])

        if results:
        
            if len(params.betas) > 1:
                for beta in params.betas:

                    grid.cell_data['Modeled Velocity']  = np.column_stack((results.ve[beta],results.vn[beta]))
                    grid.cell_data['Modeled Strain']  = np.column_stack((results.Exx[beta],results.Exy[beta],results.Eyy[beta]))

                    # save a separate file for each beta
                    filename=f"{filename[:-4]}_beta_{beta}.vtk"

                    grid.save(filename, binary = False)

                else:
                    grid.cell_data['Modeled Velocity']  = np.column_stack((results.ve[0],results.vn[0]))
                    grid.cell_data['Modeled Strain']  = np.column_stack((results.Exx[0],results.Exy[0],results.Eyy[0]))

        if (not save_data) & (not results):
            print("Warning: you are saving a mesh with no attributes.")

        grid.save(filename, binary = False)


class inversion_results:

    def __init__(self):
        # build realizations
        self.Exx = []
        self.Exy = []
        self.Eyy = []

        self.Ve = []
        self.Vn = []

        self.mhats=[]
        self.dhats=[]
        self.chi2=[]

        if params.uncertainty:
            self.Cov_Ve = []
            self.Cov_Vn = []
            self.Cov_exx = []
            self.Cov_eyy = []
            self.Cov_exy = []

    def save(self,filename=False):
        import pickle
        # if file unspecified, generate default name
        if not(filename):
            timestamp = datetime.now().strftime("%d-%m-%Y_%H:%M")
            filename = f"{timestamp}_results.pkl"
        output = open(filename,'wb')
        pickle.dump(self,output)
        output.close()


def Inversion(mesh):
    # Generates Ge_x, Ge_y, Gn_x, Gn_y
    results=inversion_results()

    d = np.concatenate([mesh.vel_centroids.Ve[mesh.ind], mesh.vel_centroids.Vn[mesh.ind]])
    sig = np.concatenate([mesh.vel_centroids.Sige[mesh.ind], mesh.vel_centroids.Sign[mesh.ind]])

    Gpoint = tools.build_G_pointforce(mesh)

    if mesh.creeping:
        Gcreep = tools.build_G_creep(mesh)

        GVe = np.concatenate([Gpoint.GVe,  Gcreep.Geast], axis=1)
        GVn = np.concatenate([Gpoint.GVn, Gcreep.Gnorth], axis=1)
        GExx = np.concatenate([Gpoint.GExx, Gcreep.GExx], axis=1)
        GExy = np.concatenate([Gpoint.GExy, Gcreep.GExy], axis=1)
        GEyy = np.concatenate([Gpoint.GEyy, Gcreep.GEyy], axis=1)

    else:
        GVe = Gpoint.GVe
        GVn = Gpoint.GVn
        GExx = Gpoint.GExx
        GExy = Gpoint.GExy
        GEyy = Gpoint.GEyy


    G = np.concatenate([GVe[mesh.ind,:],GVn[mesh.ind,:]])

    numobs = mesh.nodes.shape[0]*2
    
    L = np.eye((numobs))
    
    # if creeping...
    if mesh.creeping:
        numfaults = mesh.PatchCreepRates.shape[0]
        L = np.hstack((L,np.zeros((numobs, numfaults))))
        Lcreep = np.hstack((np.zeros((numfaults,numobs)),
                            np.eye(numfaults)))


    for beta in params.betas:

        GL = np.vstack((G/sig[:, np.newaxis], beta * L))
        d0=np.concatenate((d/sig  ,np.zeros((L.shape[0], ))))

        if mesh.creeping:
            GL=np.vstack((GL, params.Wcreep * Lcreep))
            d0=np.concatenate((d0,np.squeeze(params.Wcreep * mesh.PatchCreepRates)))


        print("Performing inversion")

        tic = time.time()
    
        mhat0 = scipy.linalg.lstsq(GL, d0, lapack_driver='gelsy')[0]

        T = time.time() - tic
        print(f"Inversion time: {T:.1f} seconds")

        if params.twostep:

            Exxs = GExx @ mhat0
            Exys = GExy @ mhat0
            Eyys = GEyy @ mhat0
            max_shear = np.sqrt((Exxs-Eyys)**2 + Exys**2)

            ind_boundaries = (
                (mesh.tri_centroids[:, 0] > params.minimize_xbound[0]) & 
                (mesh.tri_centroids[:, 0] < params.minimize_xbound[1]) & 
                (mesh.tri_centroids[:, 1] > params.minimize_ybound[0]) & 
                (mesh.tri_centroids[:, 1] < params.minimize_ybound[1])
            )

            ind_threshold = (max_shear < params.strain_threshold) & ind_boundaries
            pct_cells=sum(ind_threshold)/ind_threshold.shape[0]*100
            print(f"Minimizing strain rates in: {pct_cells:.1f}% of cells")
           
            L_iter = np.vstack((GExx[ind_threshold, :], GExy[ind_threshold, :], GEyy[ind_threshold, :]))
            GL_iter = np.vstack((GL, params.gamma * L_iter))
            d0_iter = np.concatenate((d0, np.zeros(3 * np.sum(ind_threshold))))


            tic = time.time()


            mhat = scipy.sparse.linalg.lsqr(GL_iter, d0_iter, x0=mhat0, iter_lim=1500)[0]
            T = time.time() - tic
            print(f"Completed second inversion step in: {T:.1f} seconds")

        else:
            mhat = mhat0


        dhat = G @ mhat
        
        # Propagate errors to strain rate
        if params.uncertainty:
            if params.twostep:
                if mesh.creeping:
                    Ginv = (G.T @ G + beta**2 * L.T @ L + params.Wcreep**2 * Lcreep.T @ Lcreep + params.gamma * L_iter.T @ L_iter)
                else:
                    Ginv = (G.T @ G + beta**2 * L.T @ L + params.gamma * L_iter.T @ L_iter)
            else:
                if mesh.creeping:
                    Ginv = (G.T @ G + beta**2 * L.T @ L + params.Wcreep**2 * Lcreep.T @ Lcreep)
                else:
                    Ginv = (G.T @ G + beta**2 * L.T @ L)

            tic = time.time()
            Gsharp = scipy.linalg.lstsq(Ginv,G.T,lapack_driver='gelsy')[0]
            T = time.time() - tic
            print(f"Completed uncertainty inversion in: {T:.1f} seconds")

            # Error propagation for strain rates
            Cov_bf = Gsharp @ Gsharp.T  # Note, data covariance is excluded because Gw is weighted
            Cov_exx = GExx @ Cov_bf @ GExx.T
            Cov_exy = GExy @ Cov_bf @ GExy.T
            Cov_eyy = GEyy @ Cov_bf @ GEyy.T

            # Error propagation for velocities
            Cov_Ve = GVe @ Cov_bf @ GVe.T
            Cov_Vn = GVn @ Cov_bf @ GVn.T

            results.Cov_exx.append(Cov_exx)
            results.Cov_eyy.append(Cov_eyy)
            results.Cov_exy.append(Cov_exy)
            results.Cov_Ve.append(Cov_Ve)
            results.Cov_Vn.append(Cov_Vn)

        # computed strain rates and velocities
        Exxs = GExx @ mhat
        Exys = GExy @ mhat
        Eyys = GEyy @ mhat
            
        Ves = GVe @ mhat
        Vns = GVn @ mhat

        chi2=sum((d/sig - dhat/sig)**2)/(d.shape[0])

        # for each results vector, generate a set of n (specified as params.num) noisey data points
        # noise generated using covariance matrices computed above
        if params.uncertainty:
            print("Generating uncertainty matrices. Warning: this is very slow.")
            J = np.zeros(Cov_exx.shape[0])  
            # Generate Exx realizations
            Exxs = (Exxs + np.random.multivariate_normal(J, Cov_exx, params.num)).T
            Exys = (Exys + np.random.multivariate_normal(J, Cov_exy, params.num)).T
            Eyys = (Eyys + np.random.multivariate_normal(J, Cov_eyy, params.num)).T
            Ves  = (Ves  + np.random.multivariate_normal(J, Cov_Ve,  params.num)).T
            Vns  = (Vns  + np.random.multivariate_normal(J, Cov_Vn,  params.num)).T

        results.mhats.append(mhat)
        results.dhats.append(dhat)
        results.chi2.append(chi2)

        results.Exx.append(Exxs)
        results.Exy.append(Exys)
        results.Eyy.append(Eyys)
        results.Ve.append(Ves)
        results.Vn.append(Vns)
        print(f"Chi-2 for beta = {beta}: {chi2:.2f}")
        
    return(results)



