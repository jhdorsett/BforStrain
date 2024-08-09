# This is a partial port of D. Engwirda's MATLAB Mesh2D package (2017) 
# Python refactorization by J. Dorsett (2024)
# Based on Persson and Strang's "A simple mesh generator in MATLAB"
# Ported all dependencies for smooth2 function. 

import numpy as np
import triangle as tr
from ismember import ismember
import time
import scipy.sparse as sp
import matplotlib.path as MPLpoly

# class version of makeopt function
class makeopt:
    def __init__(self):
        self.iter=32
        self.disp=4
        self.vtol=1.0e-2
        self.dbug=False

def evalhfn(vert, edge, EMAT):
    # No HFUN - HVRT is mean edge-length at vertices
    evec = vert[edge[:, 1], :] - vert[edge[:, 0], :]
    elen = np.sqrt(np.sum(evec**2, axis=1))

    # updated for sparse EMAT
    hvrt = (EMAT @ elen) / np.maximum(np.array(EMAT.sum(axis=1)).flatten(), np.finfo(float).eps)

    free = np.ones((vert.shape[0],), dtype=bool)
    free[edge[:, 0]] = False
    free[edge[:, 1]] = False

    hvrt[free] = np.inf

    return hvrt

def tricon2(tt,cc=None):
    nt=tt.shape[0]

    ee = np.zeros((nt * 3, 2),dtype='int')
    ee[:nt, :] = tt[:, [0, 1]]
    ee[nt:2*nt, :] = tt[:, [1, 2]]
    ee[2*nt:, :] = tt[:, [2, 0]]

    ee_sorted = np.sort(ee, axis=1)
    # convert int to double inside of unique function because this is faster
    ed, iv, jv = np.unique(ee_sorted[:, 0] * (2**31) + ee_sorted[:, 1], return_index=True, return_inverse=True)
    ee = ee_sorted[iv]


    tt = np.hstack((tt, np.zeros((nt, 3), dtype=int)))

    # Update the extended columns with jv values
    tt[:, 3] = jv[:nt]
    tt[:, 4] = jv[nt:2*nt]
    tt[:, 5] = jv[2*nt:3*nt]



    ne = ee.shape[0]  # Number of unique edges
    # Extend ee array with three additional columns for storing triangle indices
    ee_ext = np.hstack((ee, np.zeros((ne, 3), dtype=int)))

    # Initialize ep array with 3s
    ep = np.full(ne, 3, dtype=int)

    for ti in range(nt):
        ei = tt[ti, 3]
        ej = tt[ti, 4]
        ek = tt[ti, 5]
        
        ee_ext[ei, ep[ei]-1] = ti + 1
        ee_ext[ej, ep[ej]-1] = ti + 1
        ee_ext[ek, ep[ek]-1] = ti + 1

        ep[ei] += 1
        ep[ej] += 1
        ep[ek] += 1

    if cc is not None:
        # Check membership and get indices
        cc_sort=np.sort(cc)
        cc_sort_dec=cc_sort.astype(float)[:, 0] * (2**31) + cc_sort.astype(float)[:, 1]
        #ip = np.isin(ed, cc_sort)
        # Initialize arrays
        LX = np.isin(ed, cc_sort_dec)
        ip = np.zeros_like(ed, dtype=int)

        # Map indices where LX is True to their corresponding values in X
        indices = np.where(LX)[0]
        #ip[indices] = np.arange(1, np.sum(LX) + 1)
        ip[indices] = ismember(ed,cc_sort_dec)[1]+1
        ee_ext[:, 4] = ip

    return(ee_ext,tt)

# fixed for zero indexed
def triarea(pp, tt):
    # Extract vertices for each triangle
    verts = pp[tt]  

    # Calculate edge vectors
    ev12 = verts[:, 1, :] - verts[:, 0, :]
    ev13 = verts[:, 2, :] - verts[:, 0, :]

    # Determine dimensionality of vertices
    dim = pp.shape[1]

    # Calculate area based on dimension
    if dim == 2:
        area = 0.5 * (ev12[:, 0] * ev13[:, 1] - ev12[:, 1] * ev13[:, 0])
    elif dim == 3:
        avec = np.cross(ev12, ev13)
        area = 0.5 * np.linalg.norm(avec, axis=1)
    else:
        raise ValueError('Unsupported dimension.')

    return area

# fixed for zero indexed
def triscr2(pp, tt):
    # Compute signed area-length ratios for triangles
    scal = 4.0 * np.sqrt(3.0) / 3.0

    area = triarea(pp, tt)  # Compute triangle areas using triarea function

    # Compute squared edge lengths
    lrms = (np.sum((pp[tt[:, 1], :] - pp[tt[:, 0], :]) ** 2, axis=1) +
            np.sum((pp[tt[:, 2], :] - pp[tt[:, 1], :]) ** 2, axis=1) +
            np.sum((pp[tt[:, 2], :] - pp[tt[:, 0], :]) ** 2, axis=1))

    lrms = (lrms / 3.0) ** 1.00

    # Compute area-length ratios
    tscr = scal * area / lrms

    return tscr

def setset2(a, b):
    a_tuples = [tuple(row) for row in a]
    b_tuples = [tuple(row) for row in b]

    # Check membership
    b_set = set(b_tuples)
    is_member = np.array([row in b_set for row in a_tuples])

    return(is_member)

def constrained_delaunay(vertices,segments,flags=False): 

    # Convert vertices and segments to triangle's input format
    vertices = np.array(vertices, dtype='float64')
    segments = np.array(segments, dtype='int32')

    # sets default Delaunay behavior
    # https://www.cs.cmu.edu/~quake/triangle.switch.html
    if not flags:
        flags='pq0'

    t = tr.triangulate({'vertices': vertices,'segments': segments}, flags)
    nodes=t['vertices']
    constraints=t['segments']
    tria=t['triangles']
    return nodes,constraints,tria
  

def deltri2(vert,conn,node,PSLG,part):
    vert,conn,tria = constrained_delaunay(vert,conn)

    # calc inside status...
    tnum = np.zeros(len(tria))
     
    tmid = vert[tria[:,0],:]+vert[tria[:,1],:]+vert[tria[:,2],:]  
    tmid = tmid / 3.


    vmin = np.min(node, axis=0)
    vmax = np.max(node, axis=0)
    
    box=MPLpoly.Path(np.array([[vmin[0], vmin[1]], 
                                [vmax[0], vmin[1]], 
                                [vmax[0], vmax[1]], 
                                [vmin[0], vmax[1]]]))
    
    stat=box.contains_points(tmid)
    tnum[stat]=1
    tria = tria[tnum>0,:]
    tnum = tnum[tnum>0]


    area=triarea(vert,tria)
    tria[area<0,:]=tria[area<0][:,[0,2,1]]

    return(vert,conn,tria,tnum)

def relax_verts(vert,conn, edge,IMAT,JMAT,EMAT,iter,free):
    for isub in range(np.maximum(2,np.minimum(8,iter))):
        hvrt = evalhfn(vert,edge,EMAT)
        hmid=0.5*(hvrt[edge[:,0]]+hvrt[edge[:,1]])
        evec = vert[edge[:,1],:] - vert[edge[:,0],:]
        elen=np.sqrt(np.sum(evec**2,1))

        scal = 1-elen/hmid
        scal = np.minimum(1,scal)
        scal = np.maximum(-1,scal)

        scalscal=np.vstack((scal,scal)).T

        ipos=vert[edge[:,0],:] -.67*scalscal*evec
        jpos=vert[edge[:,1],:] +.67*scalscal*evec

        scal=np.maximum(np.abs(scal)**1,2**-52**.75)
        scalscal=np.vstack((scal,scal)).T

        vnew=IMAT@(scalscal*ipos)+JMAT@(scalscal*jpos)
        vsum=np.maximum((EMAT@scal),2**-52**.75)

        vsumvsum=np.vstack((vsum,vsum)).T

        vnew=vnew/vsumvsum

        #assert fixed pooints
        vnew[conn.reshape(-1),:]=vert[conn.reshape(-1),:]
        vnew[free,:]=vert[free,:]
        vert = vnew
    return(vnew)

def undo_verts(vert,tria,vold,original_score,iter):
    new_score = np.ones((tria.shape[0]))
    btri = np.ones((tria.shape[0]), dtype=bool)

    umax = 8
    for undo in range(1, umax+1):
        new_score[btri]=triscr2(vert,tria[btri,0:3])
        # Determine if triangles need "unwinding"
        smin = 0.70
        smax = 0.90
        sdel = 0.025

        stol = smin + (iter+1) * sdel
        stol = min(smax, stol)
        # this may need checking when iter > 1
        
        btri=((new_score <= stol) & (new_score < original_score))
        if not np.any(btri):
            break
        
        ivrt = np.unique(tria[btri, :3])
        bvrt = np.zeros((vert.shape[0],), dtype=bool)
        bvrt[ivrt] = True
            
        if undo != umax:
            bnew = 0.75 ** undo
            bold = 1.0 - bnew
        else:
            bnew = 0.0
            bold = 1.0 - bnew
            
        vert[bvrt, :] = bold * vold[bvrt, :] + bnew * vert[bvrt, :]
        
        btri = np.any(bvrt[tria[:, :3]], axis=1)

    return(vert,new_score)

def density_metric(vert,conn,edge,hvrt,free,vdeg):

    h_mid = (hvrt[edge[:, 0]] + hvrt[edge[:, 1]]) * 0.5 
    # calculate locations of edge midpoints
    edge_vec = vert[edge[:, 1], :]-vert[edge[:, 0], :] 
    edge_length = np.sqrt(np.sum(edge_vec ** 2, axis=1))
    # length ratio over spacing function. tells us if edges are too long or short
    scal = edge_length/h_mid
    # keep vector contains index of vertices that we want to keep and not remove
    keep = np.zeros(vert.shape[0], dtype=bool)
    # verticies connected to many edges
    keep[vdeg > 4] = True
    # vertices that our connecting boundary array
    keep[np.unique(conn)] = True
    # keep free surfaces (vertices that don't connect to edges)
    keep[free] = True

    # vbnd are verts that are bound because they are used in connecting boundary
    vbnd = np.zeros(vert.shape[0], dtype=bool)
    vbnd[conn[:, 0]] = True
    vbnd[conn[:, 1]] = True

    # boolean array for vertices where we are going to NOT allow density changes
    edge_bnd = vbnd[edge[:, 0]] | vbnd[edge[:, 1]]

    # Density control
    lmax = 5.0 / 4.0
    lmin = 1.0 / lmax

    # using length ratio created above, arrays less and more 
    # for where our edges are too long and short
    less = scal <= lmin
    more = scal >= lmax

    less[edge_bnd] = False
    more[edge_bnd] = False

    less_index = np.where(less)[0]
    for less_position in range(len(less_index)):
        edge_position = less_index[less_position]
        vert_A = edge[edge_position, 0]
        vert_B = edge[edge_position, 1]
        # If still disjoint
        if keep[vert_A] and keep[vert_B]:
            keep[vert_A] = False
            keep[vert_B] = False
        else:
            less[edge_position] = False

    edge_bnd[less] = keep[edge[less, 0]] & keep[edge[less, 1]]
    more[edge_bnd] = False

    return(less,more,keep)

def reindex_verts(vert,edge,less,keep,edge_mid,tria):
    # use -1 as my "empty" values, since 0 is a real index
    redo = np.zeros(vert.shape[0], dtype=int)-1

    # Identify keep and less indices
    itop = np.sum(keep)
    iend = np.sum(less)

    redo[keep] = np.arange(itop)
    redo[edge[less, 0]] = np.arange(itop , itop + iend)
    redo[edge[less, 1]] = np.arange(itop, itop + iend)


    # Create new vertices and triangles
    vnew = np.vstack((vert[keep, :], edge_mid[less, :]))
    tnew = redo[tria[:, :3]]
    ttmp=np.sort(tnew,axis=1)
    okay = np.all(np.diff(ttmp, axis=1) != 0, axis=1)
    okay = okay & (ttmp[:, 0] >= 0)
    tnew=tnew[okay,:]

    return(vnew,tnew,okay,redo)

def make_MATS(vert,edge):
    
    nvert = vert.shape[0]
    nedg = edge.shape[0]

    # Initialize lists to store data for the sparse matrices
    IMAT_data = []
    IMAT_row = []
    IMAT_col = []

    JMAT_data = []
    JMAT_row = []
    JMAT_col = []

    # Fill the data for IMAT and JMAT
    for i in range(nedg):
        IMAT_row.append(edge[i, 0])
        IMAT_col.append(i)
        IMAT_data.append(1)

        JMAT_row.append(edge[i, 1])
        JMAT_col.append(i)
        JMAT_data.append(1)

    IMAT = sp.csr_matrix((IMAT_data, (IMAT_row, IMAT_col)), shape=(nvert, nedg))
    JMAT = sp.csr_matrix((JMAT_data, (JMAT_row, JMAT_col)), shape=(nvert, nedg))
    return(IMAT, JMAT)

    
def smooth2(mesh): 

    vert = mesh.nodes
    conn = mesh.edges
    tria = mesh.tri

    print()
    print(" Smooth triangulation...")
    print()
    print(" -----------------------------------------------------------------------")
    print("      |ITER.|          |MOVE(X)|          |DTRI(X)|          |T. ITER|  ")
    print(" -----------------------------------------------------------------------")

    # jump into smooth2!
    tnum=np.ones(tria.shape[0])
    node=vert
    opts=makeopt()
    pmax=max(tnum).astype(int)
    # currently, don't preallocate part

    # it's dynamic but only has size 1 for smooth2
    part = {}
    PSLG = conn
    for ppos in range(1, pmax + 1):
        tsel = (tnum == ppos)
        tcur = tria[tsel, :]
        ecur, tcur = tricon2(tcur)
        ebnd = (ecur[:, 3] == 0)
        same = setset2(PSLG, ecur[ebnd, :2])
        part[ppos-1] = np.where(same)[0]

    vmin = np.min(vert, axis=0)
    vmax = np.max(vert, axis=0)

    # Compute the adjustment values
    vdel = vmax - 1.0 * vmin
    vmin = vmin - 0.5 * vdel
    vmax = vmax + 0.5 * vdel

    # Create the bounding box vertices
    vbox = np.array([
        [vmin[0], vmin[1]],
        [vmax[0], vmin[1]],
        [vmax[0], vmax[1]],
        [vmin[0], vmax[1]]
    ])
    vert = np.vstack((vert,vbox))

    t = 0
    for iter in range(opts.iter):
        t_0 = time.time()
    
        edge,tria = tricon2(tria,conn)
        
        original_score = triscr2(vert,tria[:,0:3])

        IMAT, JMAT = make_MATS(vert,edge)

        EMAT = IMAT + JMAT

        vdeg = np.array(EMAT.sum(axis=1)).flatten()

        free = (vdeg == 0)

        # calculate locations of edge midpoints
        edge_beg = vert[edge[:, 0], :] 
        edge_end = vert[edge[:, 1], :] 
        #edge_vec = edge_beg-edge_end
        #edge_length = np.sqrt(np.sum(edge_vec ** 2, axis=1))

        # save old vertices for later...
        old_vert=vert
        # perform the relaxation step, moving verts using distance function and fixing at IMAT, JMAT, EMAT
        vert = relax_verts(vert, conn, edge, IMAT, JMAT, EMAT, iter,free)
        # if the quality metric has gotten substantially worse after relax_verts, 
        # replace vert with old_vert
        vert, new_score = undo_verts(vert,tria,old_vert,original_score,iter)

        # update the score for this larger iteration, overwrite original_score
        original_score = new_score

        # save vdel for later...
        vdel = np.sum((vert-old_vert)**2,axis=1)
        #central point of each edge
        edge_mid = (vert[edge[:, 0], :] + vert[edge[:, 1], :])*0.5

        # compute h function at updated vertices and store the values at midpoints
        hvrt = evalhfn(vert, edge, EMAT)

        # compute areas in the mesh that should be down/upsampled based on geometry + h function
        less, more, keep = density_metric(vert,conn,edge,hvrt,free,vdeg)

        # in those areas, go ahead and dump/upscale vertices and update triangles
        vnew,tnew,okay,redo = reindex_verts(vert,edge,less,keep,edge_mid,tria)

        # do a quick check to see if they're bad...
        new_score = triscr2(vnew,tnew)

        stol = 0.8
        tbad = (new_score < stol) & (new_score < original_score[okay])
        # where they are bad, store the array
        vbad = np.zeros(vnew.shape[0], dtype=bool)
        vbad[tnew[tbad, :]] = True

        # Filter edge merge
        lidx = np.where(less)[0]
        # ZERO INDEXED FOR PYTHON
        ebad = vbad[redo[edge[lidx, 0]]] | vbad[redo[edge[lidx, 1]]]
        less[lidx[ebad]] = False
        keep[edge[lidx[ebad], :2]] = True

        # Reindex vert/conn
        redo = np.zeros(vert.shape[0], dtype=int)

        itop = np.sum(keep)
        iend = np.sum(less)

        redo[keep] = np.arange(itop)
        redo[edge[less, 0]] = np.arange(itop , itop + iend)
        redo[edge[less, 1]] = np.arange(itop, itop + iend)

        vert = np.vstack((vert[keep, :], edge_mid[less, :], edge_mid[more, :]))
        conn = redo[conn[:, :2]]
        vert, conn, tria, tnum = deltri2(vert,conn,node,PSLG,part)


        vdel= vdel/ (hvrt**2)
        move = vdel > opts.vtol**2

        ntri=len(tria)
        nmov=sum(move)

        t = t + time.time()-t_0

        if np.mod(iter+1,opts.disp) == 0:
            #print(iter+1,nmov,ntri,t/opts.disp)
            print(f"{iter+1:11d} {nmov:18d} {ntri:18d} {(round(t/opts.disp,3)):20f}")
            t = 0
        
        if nmov == 0:
            break
        
    keep = np.zeros(len(vert), dtype=bool)
    keep[tria] = True
    keep[conn] = True

    redo = np.zeros(len(vert), dtype=int)
    redo[keep] = np.arange(np.sum(keep))

    # Update 'conn' and 'tria' using 'redo'
    conn = redo[conn]
    tria = redo[tria]

    # Update 'vert' using 'keep'
    vert = vert[keep, :]

    mesh.nodes = vert
    mesh.edges = conn
    mesh.tri = tria
    mesh.tnum = tnum
    mesh.score = triscr2(vert,tria)
    #return(vert,conn,tria,tnum)
    return(mesh)
