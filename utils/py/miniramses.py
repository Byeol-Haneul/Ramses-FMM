import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.io import FortranFile
from astropy.io import ascii
import os

import time

class Cool:
    """
    This is the class for RAMSES cooling table.
    """
    def __init__(self,n1,n2):
        """
        This function initialize the cooling table. 
        Args:
            n1: number of points for the gas density axis
            n2: number of points for the gas temperature axis
        """
        self.n1 = n1
        self.n2 = n2
        self.nH = np.zeros([n1])
        self.T2 = np.zeros([n2])
        self.cool = np.zeros([n1,n2])
        self.heat = np.zeros([n1,n2])
        self.spec = np.zeros([n1,n2,6])
        self.xion = np.zeros([n1,n2])

def clean(dat,n1,n2):
    dat = np.array(dat)
    dat = dat.reshape(n2, n1)
    return dat

def clean_spec(dat,n1,n2):
    dat = np.array(dat)
    dat = dat.reshape(6, n2, n1)
    return dat

def rd_cool(filename):
    """This function reads a RAMSES cooling table file (unformatted Fortran binary) 
    and store it in a cooling object.

    Args:
        filename: the complete path (including the name) of the cooling table file.

    Returns:
        A cooling table (Cool) object.

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filename, 'r') as f:
        n1, n2 = f.read_ints('i')
        c = Cool(n1,n2)
        nH = f.read_reals('f8')
        T2 = f.read_reals('f8')
        cool = f.read_reals('f8')
        heat = f.read_reals('f8')
        cool_com = f.read_reals('f8')
        heat_com = f.read_reals('f8')
        metal = f.read_reals('f8')
        cool_prime = f.read_reals('f8')
        heat_prime = f.read_reals('f8')
        cool_com_prime = f.read_reals('f8')
        heat_com_prime = f.read_reals('f8')
        metal_prime = f.read_reals('f8')
        mu = f.read_reals('f8')
        n_spec = f.read_reals('f8')
        c.nH = nH
        c.T2 = T2
        c.cool = clean(cool,n1,n2)
        c.heat = clean(heat,n1,n2)
        c.metal = clean(metal,n1,n2)
        c.spec = clean_spec(n_spec,n1,n2)
        c.xion = c.spec[0]
        for i in range(0,n2):
            c.xion[i,:] = c.spec[0,i,:] - c.nH
        return c

class Map:
    """This class defines a map object.
    """
    def __init__(self,nx,ny):
        """This function initalize a map object.
        
        Args:
            nx: number of pixels in the x direction
            ny: number of pixels in the y direction
        """
        self.nx = nx
        self.ny = ny
        self.data = np.zeros([nx,ny])

def rd_map(filename):
    """This function reads a RAMSES map file (unformatted Fortran binary)
    as produced by the RAMSES utilities amr2map or part2map and store it in a map object.

    Args:
        filename: the complete path (including the name) of the map file.

    Returns:
        A map (class Map) object.

    Example:
        import miniramses as ram
        map = ram.rd_map("dens.map")
        plt.imshow(map.data,origin="lower")

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filename, 'r') as f:
        t, dx, dy, dz = f.read_reals('f8')
        nx, ny = f.read_ints('i')
        dat = f.read_reals('f4')
    
    dat = np.array(dat)
    dat = dat.reshape(ny, nx)
    m = Map(nx,ny)
    m.data = dat
    m.time = t
    m.nx = nx
    m.ny = ny
    
    return m

class Histo:
    """This class defines a histogram object.
    """
    def __init__(self,nx,ny):
        """This function initalize a histogram object.

        Args:
            nx: number of pixels in the x direction
            ny: number of pixels in the y direction
        """
        self.nx = nx
        self.ny = ny
        self.h = np.zeros([nx,ny])

def rd_histo(filename):
    """This function reads a RAMSES histogram file (unformatted Fortran binary)
    as produced by the RAMSES utilities histo and store it in a Histo object.

    Args:
        filename: the complete path (including the name) of the histo file.

    Returns:
        A histogram (class Histo) object.

    Example:
        import miniramses as ram
        h = ram.rd_histo("histo.dat")
        plt.imshow(h.data,origin="lower")
    
    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filename, 'r') as f:
        nx, ny = f.read_ints('i')
        dat = f.read_reals('f4')
        lxmin, lxmax = f.read_reals('f8')
        lymin, lymax = f.read_reals('f8')

    dat = np.array(dat)
    dat = dat.reshape(ny, nx)
    h = Histo(nx,ny)
    h.data = dat
    h.nx = nx
    h.ny = ny
    h.lxmin = lxmin
    h.lxmax = lxmax
    h.lymin = lymin
    h.lymax = lymax

    return h

class Part:
    def __init__(self,nnp,nndim,star=False,sink=False,peak=False):
        self.np = nnp
        self.ndim = nndim
        self.xp = np.zeros([nndim,nnp])
        self.vp = np.zeros([nndim,nnp])
        self.mp = np.zeros([nnp])
        if(star):
            self.zp = np.zeros([nnp])
            self.tp = np.zeros([nnp])
        if(sink):
            self.fp = np.zeros([nndim,nnp])
            self.tp = np.zeros([nnp])
        if(peak):
            self.pid = np.zeros([nnp],dtype=np.int32)
            
def rd_part(nout,**kwargs):
    """This function reads a RAMSES particle file (unformatted Fortran binary) 
    as produced by the RAMSES code in the snapshot directory output_00* 
    and store it in a variable containing all the particle information (Part object).

    Args:
        nout: the RAMSES snapshot number. For example output_000012 corresponds to nout=12.

    Optional args:

        center: a numpy array containing the coordinates of the center of the sphere restricting
                the region to read in data.

        radius: the radius of the sphere restricting the region to read in data.

    Returns:
        A variable p (class Part) object defined as:
            p.np: number of particles
            p.ndim: number of space dimensions
            p.xp: coordinates of the particles. p.xp[0] gives the x coordinate as a numpy array.
            p.vp: velocities of the particles. p.vp[0] gives the x-component as a numpy array.
            p.mp: array containing the particle masses

    Example:
        import miniramses as ram
        p = ram.rd_part(12,center=[0.5,0.5,0.5],radius=0.1)
        print(np.max(p.xp[0]))
    
    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    
    backup = kwargs.get("backup",False)
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    path = kwargs.get("path","./")
    star = kwargs.get("star",False)
    sink = kwargs.get("sink",False)
    peak = kwargs.get("peak",False)

    car1 = str(nout).zfill(5)
    i = rd_info(nout,path=path,backup=backup)
    ncpu = i.ncpu
    ndim = i.ndim
    levelmin = i.levelmin
    nlevelmax = i.nlevelmax

    #if ( not (center is None)  and not (radius is None) ):
    #    info = rd_info(nout)
    #    cpulist = get_cpu_list(info,**kwargs)
    #    print("Will open only",len(cpulist),"files")
    #else:
    #    cpulist = range(1,ncpu+1)
    cpulist = range(1,ncpu+1)

    prefix="/part."        
    if(star):
        prefix="/star."
    if(sink):
        prefix="/sink."

    npart = 0
    for icpu in cpulist:
        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+prefix+car2
        else:
            filename = path+"/output_"+car1+prefix+car2

        npart2 = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
        npart = npart + npart2

    txt = "Found "+str(npart)+" particles"
    print(txt)
    print("Reading particle data...")

    p = Part(npart,ndim,star,sink,peak)
    p.np = npart
    p.ndim = ndim

    ipart = 0
    for	icpu in	cpulist:
        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+prefix+car2
        else:
            filename = path+"/output_"+car1+prefix+car2

        npart2 = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
        
        for idim in range(0,ndim):
            if(backup):
                offset = 8+idim*npart2*8
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
            else:
                offset = 8+idim*npart2*4
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)

            p.xp[idim,ipart:ipart+npart2] = xp
            
        for idim in range(0,ndim):
            if(backup):
                offset = 8+npart2*8*ndim+idim*npart2*8
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
            else:
                offset = 8+npart2*4*ndim+idim*npart2*4
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)

            p.vp[idim,ipart:ipart+npart2] = xp

        if(backup):
            offset = 8+npart2*8*ndim*2
            xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
        else:
            offset = 8+npart2*4*ndim*2
            xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)

        p.mp[ipart:ipart+npart2] = xp

        if(star):
            if(backup):
                offset = 8+npart2*8*ndim*2+npart2*8
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
            else:
                offset = 8+npart2*4*ndim*2+npart2*4
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)

            p.zp[ipart:ipart+npart2] = xp

            if(backup):
                offset = 8+npart2*8*ndim*2+2*npart2*8
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
            else:
                offset = 8+npart2*4*ndim*2+2*npart2*4
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)

            p.tp[ipart:ipart+npart2] = xp

        if(sink):
            for idim in range(0,ndim):
                if(backup):
                    offset = 8+idim*npart2*8
                    xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
                else:
                    offset = 8+idim*npart2*4
                    xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)

                p.fp[idim,ipart:ipart+npart2] = xp

            if(backup):
                offset = 8+npart2*8*ndim*2+2*npart2*8
                xp = np.fromfile(filename,dtype=np.float64,count=npart2,offset=offset)
            else:
                offset = 8+npart2*4*ndim*2+2*npart2*4
                xp = np.fromfile(filename,dtype=np.float32,count=npart2,offset=offset)

            p.tp[ipart:ipart+npart2] = xp

        ipart = ipart + npart2

    if(peak):
        prefix="/peak_part."        
        if(star):
            prefix="/peak_star."
        if(sink):
            prefix="/peak_sink."
        ipart = 0
        for icpu in cpulist:
            car2 = str(icpu).zfill(5)
            filename = path+"/output_"+car1+prefix+car2
            npart2 = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
            offset = 8
            pid = np.fromfile(filename,dtype=np.int32,count=npart2,offset=offset)
            p.pid[ipart:ipart+npart2] = pid
            ipart = ipart + npart2
    
    if ( not (center is None)  and not (radius is None) ):
        # Filtering particles
        r = np.sqrt((p.xp[0]-center[0])**2+(p.xp[1]-center[1])**2+(p.xp[2]-center[2])**2)
        p.np = np.count_nonzero(r < radius)
        p.mp = p.mp[r < radius]
        p.xp = p.xp[:,r < radius]
        p.vp = p.vp[:,r < radius]
        if(star):
            p.zp = p.zp[r < radius]
            p.tp = p.tp[r < radius]
        if(sink):
            p.tp = p.tp[r < radius]
        if(peak):
            p.pid = p.pid[r < radius]

    return p

class Level:
    def __init__(self,nndim):
        self.level = 0
        self.ngrid = 0
        self.ndim = nndim
        self.xg = np.empty(shape=(nndim,0))
        self.refined = np.empty(shape=(2**nndim,0),dtype=bool)

def rd_amr(nout,**kwargs):

    backup = kwargs.get("backup",False)
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    path = kwargs.get("path","./")

    car1 = str(nout).zfill(5)
    i = rd_info(nout,path=path,backup=backup)
    ncpu = i.ncpu
    ndim = i.ndim
    levelmin = i.levelmin
    nlevelmax = i.nlevelmax

    txt = "ncpu="+str(ncpu)+" ndim="+str(ndim)+" nlevelmax="+str(nlevelmax)
    print(txt)
    print("Time=",i.texp)
    print("Reading grid data...")

    #if ( not (center is None)  and not (radius is None) ):
    #    info = rd_info(nout)
    #    cpulist = get_cpu_list(info,**kwargs)
    #    print("Will open only",len(cpulist),"files")
    #else:
    #    cpulist = range(1,ncpu+1)
    cpulist = range(1,ncpu+1)
    
    amr=[]
    for ilevel in range(0,nlevelmax):
        amr.append(Level(ndim))
        
    amr[0].boxlen = i.boxlen
    
    numbl = np.zeros([nlevelmax,ncpu],dtype=np.int32)
    
    # Reading and computing total AMR grids count
    for icpu in cpulist:

        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)

        if(backup):
            filename = path+"/backup_"+car1+"/amr."+car2
        else:
            filename = path+"/output_"+car1+"/amr."+car2

        skip = 12
        for ilevel in range(levelmin-1,nlevelmax):
            offset = skip+4*(ilevel+1-levelmin)
            numbl[ilevel,icpu-1] = np.fromfile(filename,dtype=np.int32,count=1,offset=offset)[0]
            amr[ilevel].ngrid = amr[ilevel].ngrid + numbl[ilevel,icpu-1]

    # Allocating memory
    for ilevel in range(0,nlevelmax):
        amr[ilevel].xg = np.zeros([ndim,amr[ilevel].ngrid],dtype=float)
        amr[ilevel].refined = np.zeros([2**ndim,amr[ilevel].ngrid],dtype=bool)

    iskip = np.zeros(nlevelmax, dtype=int)
    nvar = ndim+2**ndim

    # Reading and storing data
    for icpu in cpulist:

        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+"/amr."+car2
        else:
            filename = path+"/output_"+car1+"/amr."+car2

        offset = 12 + 4*(nlevelmax+1-levelmin)
        for ilevel in range(levelmin-1,nlevelmax):
            ncache = numbl[ilevel,icpu-1]

            transfer = np.fromfile(filename,dtype=np.int32,count=nvar*ncache,offset=offset)
            transfer = np.reshape(transfer,(ncache,nvar))
            transfer = np.transpose(transfer)

            # Store grid Cartesian index
            for idim in range(0,ndim):
                amr[ilevel].xg[idim,iskip[ilevel]:iskip[ilevel]+ncache] = transfer[idim]
                
            # Store cell refinement map
            for ind in range(0,2**ndim):
                amr[ilevel].refined[ind,iskip[ilevel]:iskip[ilevel]+ncache] = transfer[ndim+ind]
            
            offset = offset + ncache*nvar*4
            iskip[ilevel] = iskip[ilevel] + ncache
            
    return amr

class Hydro:
    def __init__(self,nndim,nnvar):
        self.level = 0
        self.ngrid = 0
        self.ndim = nndim
        self.nvar = nnvar
        self.u = np.empty(shape=(nnvar,2**nndim,0))

def rd_hydro(nout,**kwargs):

    prefix = kwargs.get("prefix","hydro")
    backup = kwargs.get("backup",False)
    center = kwargs.get("center")
    radius = kwargs.get("radius")
    path = kwargs.get("path","./")

    car1 = str(nout).zfill(5)
    i = rd_info(nout,path=path,backup=backup)
    ncpu = i.ncpu
    ndim = i.ndim
    levelmin = i.levelmin
    nlevelmax = i.nlevelmax

    #if ( not (center is None)  and not (radius is None) ):
    #    info = rd_info(nout)
    #    cpulist = get_cpu_list(info,**kwargs)
    #    print("Will open only",len(cpulist),"files")
    #else:
    #    cpulist = range(1,ncpu+1)
    cpulist = range(1,ncpu+1)
    
    # Get number of hydro variables
    car1 = str(nout).zfill(5)
    if(backup):
        filename = path+"/backup_"+car1+"/"+prefix+".00001"
    else:
        filename = path+"/output_"+car1+"/"+prefix+".00001"

    nvar = np.fromfile(filename,dtype=np.int32,count=1,offset=4)[0]
    
    txt = "ncpu="+str(ncpu)+" ndim="+str(ndim)+" nlevelmax="+str(nlevelmax)+" nvar="+str(nvar)
    print(txt)
    print("Reading "+prefix+" data...")

    hydro=[]
    for ilevel in range(0,nlevelmax):
        hydro.append(Hydro(ndim,nvar))
        hydro[ilevel].level = ilevel
        
    numbl = np.zeros([nlevelmax,ncpu],dtype=np.int32)
    
    # Reading and computing total AMR grids count
    for icpu in cpulist:

        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+"/"+prefix+"."+car2
        else:
            filename = path+"/output_"+car1+"/"+prefix+"."+car2

        skip = 16
        for ilevel in range(levelmin-1,nlevelmax):
            offset = skip+4*(ilevel+1-levelmin)
            numbl[ilevel,icpu-1] = np.fromfile(filename,dtype=np.int32,count=1,offset=offset)[0]
            hydro[ilevel].ngrid = hydro[ilevel].ngrid + numbl[ilevel,icpu-1]

    # Allocating memory
    for ilevel in range(0,nlevelmax):
        hydro[ilevel].u = np.zeros([nvar,2**ndim,hydro[ilevel].ngrid],dtype=float)
        hydro[ilevel].nvar = nvar

    iskip = np.zeros(nlevelmax, dtype=int)
    nvartot = nvar*2**ndim
    
    # Reading and storing data
    for icpu in cpulist:

        car2 = str(icpu).zfill(5)
        if(backup):
            filename = path+"/backup_"+car1+"/"+prefix+"."+car2
        else:
            filename = path+"/output_"+car1+"/"+prefix+"."+car2

        offset = 16 + 4*(nlevelmax+1-levelmin)
        
        for ilevel in range(levelmin-1,nlevelmax):
            ncache = numbl[ilevel,icpu-1]

            if(backup):
                transfer = np.fromfile(filename,dtype=np.float64,count=nvartot*ncache,offset=offset)
            else:
                transfer = np.fromfile(filename,dtype=np.float32,count=nvartot*ncache,offset=offset)

            transfer = np.reshape(transfer,(ncache,nvar,2**ndim))            
            transfer = np.transpose(transfer,(1,2,0))

            # Store cell hydro variables
            for ivar in range(0,nvar):
                for ind in range(0,2**ndim):
                    hydro[ilevel].u[ivar,ind,iskip[ilevel]:iskip[ilevel]+ncache] = transfer[ivar,ind]

            if(backup):
                offset = offset + ncache*nvartot*8
            else:
                offset = offset + ncache*nvartot*4

            iskip[ilevel] = iskip[ilevel] + ncache

    return hydro

class Cell:
    def __init__(self,nndim,nnvar):
        self.ncell = 0
        self.ndim = nndim
        self.nvar = nnvar
        self.x = np.empty(shape=(nndim,0))
        self.u = np.empty(shape=(nnvar,0))
        self.dx = np.empty(shape=(0))

def rd_cell(nout,**kwargs):
    """This function reads RAMSES AMR and hydro files (unformatted Fortran binary) 
    as produced by the RAMSES code in the snapshot directory output_00* 
    and store it in a variable containing all the hydro leaf cells information (Cell object).

    Args:
        nout: the RAMSES snapshot number. For example output_000012 corresponds to nout=12.

    Optional args:

        center: a numpy array containing the coordinates of the center of the sphere restricting the region to read in data.

        radius: the radius of the sphere restricting the region to read in data.

    Returns:
        A variable c (class Cell) object defined as:
            c.ncell: number of AMR cells
            c.ndim: number of space dimensions
            c.nvar: number of hydro variables
            c.x: coordinates of the cells. c.x[0] gives the x coordinate as a numpy array.
            c.u: hydro variables in each cell. For example, c.u[0] gives the gas density as a numpy array.
            c.dx: array containing the individual AMR cell sizes.

    Example:
        import miniramses as ram
        c = ram.rd_cell(12,center=[0.5,0.5,0.5],radius=0.1)
        print(np.max(c.dx))

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    
    path = kwargs.get("path","./")
    center = kwargs.get("center")
    radius = kwargs.get("radius")

    a = rd_amr(nout,**kwargs)
    h = rd_hydro(nout,**kwargs)

    nlevelmax = len(a)
    ndim = a[0].ndim
    nvar = h[0].nvar
    boxlen = a[0].boxlen
    
    offset = np.zeros([ndim,2**ndim])
    if (ndim == 1):
        offset[0,:]=[-0.5,0.5]
    if (ndim == 2):
        offset[0,:]=[-0.5,0.5,-0.5,0.5]
        offset[1,:]=[-0.5,-0.5,0.5,0.5]
    if (ndim == 3):
        offset[0,:]=[-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5]
        offset[1,:]=[-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5]
        offset[2,:]=[-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5]

    ncell = 0
    for ilev in range(0,nlevelmax):
        ncell = ncell + np.count_nonzero(a[ilev].refined == False)

    print("Found",ncell,"leaf cells")
    print("Extracting leaf cells...")

    c = Cell(ndim,nvar)
    c.ncell = ncell
    
    for ilev in range(0,nlevelmax):
        dx = 0.5*boxlen/2**ilev
        for ind in range(0,2**ndim):
            nc = np.count_nonzero(a[ilev].refined[ind] == False)
            if (nc > 0):
                xc = np.zeros([ndim,nc])
                for idim in range(0,ndim):
                    xc[idim,:]= (2*a[ilev].xg[idim,np.where(a[ilev].refined[ind] == False)]+1+offset[idim,ind])*dx
                c.x = np.append(c.x,xc,axis=1)
                uc = np.zeros([nvar,nc])
                for ivar in range(0,nvar):
                    uc[ivar,:]= h[ilev].u[ivar,ind,np.where(a[ilev].refined[ind] == False)]
                c.u = np.append(c.u,uc,axis=1)
                dd = np.ones(nc)*dx
                c.dx = np.append(c.dx,dd)

    if ( not (center is None)  and not (radius is None) ):
        # Filtering cells
        r = np.sqrt((c.x[0]-center[0])**2+(c.x[1]-center[1])**2+(c.x[2]-center[2])**2) - dx
        c.ncell = np.count_nonzero(r < radius)
        c.u  = c.u[:,r < radius]
        c.x  = c.x[:,r < radius]
        c.dx = c.dx[r < radius]

    if(ndim==1):
        c.x = c.x[0]
        ind = np.argsort(c.x)
        c.x = c.x[ind]
        c.dx = c.dx[ind]
        for  ivar in range(0,nvar):
            c.u[ivar]=c.u[ivar,ind]
        
    return c

def save_cell(c,filename):

    with open(filename,'wb') as f:
        np.save(f,c.ncell)
        np.save(f,c.ndim)
        np.save(f,c.nvar)
        np.save(f,c.dx)
        np.save(f,c.x)
        np.save(f,c.u)

def load_cell(filename):

    with open(filename,'rb') as f:
        ncell = np.load(f)
        ndim = np.load(f)
        nvar = np.load(f)
        c = Cell(ndim,nvar)
        c.ncell = ncell
        c.ndim = ndim
        c.nvar = nvar
        c.dx = np.append(c.dx,np.load(f))
        c.x  = np.append(c.x, np.load(f),axis=1)
        c.u  = np.append(c.u, np.load(f),axis=1)

    return c

class Snap1d:
    pass

def rd_log(filename,**kwargs):
    """This function reads the standard ouput (aka log file)
    as produced by the RAMSES code for 1D simulations.

    Args:
        filename: the log file name (usually run.log).

    Optional args:

        None

    Returns:
        A variable r (class Run) object defined as:
            c.ncell: number of AMR cells.
            c.lev: level of refinement of the cells.
            c.x: coordinates of the cells.
            c.d: density of the cells.
            c.u: velocity field x-component.
            c.v: velocity field y-component.
            c.w: velocity field z-component.
            c.d: pressure of the cells.
            c.A: magnetic field x-component.
            c.B: magnetic field y-component.
            c.C: magnetic field z-component.

    Example:
        import miniramses as ram
        r = ram.rd_log("run.log")
        plt.plot(r["x"],r["d"]))
    
    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    cmd="grep -n Output "+filename+" > /tmp/out.txt"
    os.system(cmd)
    lines = ascii.read("/tmp/out.txt")
    print("Found "+str(len(lines))+" output(s)")

    r=[]
    for out in range(0,len(lines)):

        i = int(lines["col1"][out][:-1])
        n = int(lines["col3"][out])

        data = ascii.read(filename,header_start=i-3,data_start=i-2,data_end=i+n-2)

        r.append(Snap1d())
        r[out].ncell=n
        r[out].lev=np.array(data["lev"],dtype='int')
        r[out].x=np.array(data["x"])
        r[out].d=np.array(data["d"])
        r[out].p=np.array(data["P"])
        r[out].u=np.array(data["u"])
        if len(data.columns)>5:
            r[out].v=np.array(data["v"])
            r[out].w=np.array(data["w"])
            r[out].A=np.array(data["A"])
            r[out].B=np.array(data["B"])
            r[out].C=np.array(data["C"])

    return r

class Info:
    def __init__(self,nncpu):
        self.bound_key = np.zeros(shape=(nncpu+1),dtype=np.double)
        
def rd_info(nout,**kwargs):
    
    backup = kwargs.get("backup",False)
    path = kwargs.get("path","./")

    car1 = str(nout).zfill(5)
    if(backup):
        filename = path+"/backup_"+car1+"/info.txt"
    else:
        filename = path+"/output_"+car1+"/info.txt"

    info=ascii.read(filename,delimiter="=",format='no_header')

    ncpu=int(info[0][1])

    i = Info(ncpu)

    i.ncpu=ncpu
    i.ndim=int(info[2][1])
    i.levelmin=int(info[3][1])
    i.nlevelmax=int(info[4][1])
    i.boxlen=info[7][1]
    i.time=info[8][1]
    i.texp=info[9][1]
    i.aexp=info[10][1]
    i.gamma=info[16][1]
    i.unit_l=info[17][1]
    i.unit_d=info[18][1]
    i.unit_t=info[19][1]

    return i

def hilbert3d(x,y,z,bit_length):
    
    state_diagram = [ 1, 2, 3, 2, 4, 5, 3, 5,
                      0, 1, 3, 2, 7, 6, 4, 5,
                      2, 6, 0, 7, 8, 8, 0, 7,
                      0, 7, 1, 6, 3, 4, 2, 5,
                      0, 9,10, 9, 1, 1,11,11,
                      0, 3, 7, 4, 1, 2, 6, 5,
                      6, 0, 6,11, 9, 0, 9, 8,
                      2, 3, 1, 0, 5, 4, 6, 7,
                      11,11, 0, 7, 5, 9, 0, 7,
                      4, 3, 5, 2, 7, 0, 6, 1,
                      4, 4, 8, 8, 0, 6,10, 6,
                      6, 5, 1, 2, 7, 4, 0, 3,
                      5, 7, 5, 3, 1, 1,11,11,
                      4, 7, 3, 0, 5, 6, 2, 1,
                      6, 1, 6,10, 9, 4, 9,10,
                      6, 7, 5, 4, 1, 0, 2, 3,
                      10, 3, 1, 1,10, 3, 5, 9,
                      2, 5, 3, 4, 1, 6, 0, 7,
                      4, 4, 8, 8, 2, 7, 2, 3,
                      2, 1, 5, 6, 3, 0, 4, 7,
                      7, 2,11, 2, 7, 5, 8, 5,
                      4, 5, 7, 6, 3, 2, 0, 1,
                      10, 3, 2, 6,10, 3, 4, 4,
                      6, 1, 7, 0, 5, 2, 4, 3]

    state_diagram = np.array(state_diagram)
    state_diagram = state_diagram.reshape((8,2,12),order='F')

    n = len(x)
    order = np.zeros(n,dtype="double")
    x_bit_mask = np.zeros(bit_length  ,dtype="bool")
    y_bit_mask = np.zeros(bit_length  ,dtype="bool")
    z_bit_mask = np.zeros(bit_length  ,dtype="bool")
    i_bit_mask = np.zeros(3*bit_length,dtype=bool)
    
    for ip in  range(0,n):
        
        for i in range(0,bit_length):
            x_bit_mask[i] = x[ip] & (1 << i)
            y_bit_mask[i] = y[ip] & (1 << i)
            z_bit_mask[i] = z[ip] & (1 << i)
            
        for i in range(0,bit_length):
            i_bit_mask[3*i+2] = x_bit_mask[i]
            i_bit_mask[3*i+1] = y_bit_mask[i]
            i_bit_mask[3*i  ] = z_bit_mask[i]
            
        cstate = 0
        for i in range(bit_length-1,-1,-1):
            b2 = 0
            if (i_bit_mask[3*i+2]):
                b2 = 1
            b1 = 0
            if (i_bit_mask[3*i+1]):
                b1 = 1
            b0 = 0
            if (i_bit_mask[3*i  ]):
                b0 = 1
            sdigit = b2*4 + b1*2 + b0
            nstate = state_diagram[sdigit,0,cstate]
            hdigit = state_diagram[sdigit,1,cstate]
            i_bit_mask[3*i+2] = hdigit & (1 << 2)
            i_bit_mask[3*i+1] = hdigit & (1 << 1)
            i_bit_mask[3*i  ] = hdigit & (1 << 0)
            cstate = nstate
            
        order[ip]= 0
        for i in range(0,3*bit_length):
            b0 = 0
            if (i_bit_mask[i]):
                b0 = 1
            order[ip] = order[ip] + float(b0)*2.**i
                
    return order

def hilbert2d(x,y,bit_length):
    
    state_diagram = [ 1, 0, 2, 0, 
                      0, 1, 3, 2, 
                      0, 3, 1, 1, 
                      0, 3, 1, 2, 
                      2, 2, 0, 3, 
                      2, 1, 3, 0, 
                      3, 1, 3, 2, 
                      2, 3, 1, 0 ]
    
    state_diagram = np.array(state_diagram)    
    state_diagram = state_diagram.reshape((4,2,4), order='F')
    
    n = len(x)
    order = np.zeros(n,dtype="double")
    x_bit_mask = np.zeros(bit_length  ,dtype="bool")
    y_bit_mask = np.zeros(bit_length  ,dtype="bool")
    i_bit_mask = np.zeros(2*bit_length,dtype=bool)
    
    for ip in  range(0,n):
        
        for i in range(0,bit_length):
            x_bit_mask[i] = bool(x[ip] & (1 << i))
            y_bit_mask[i] = bool(y[ip] & (1 << i))
            
        for i in range(0,bit_length):
            i_bit_mask[2*i+1] = x_bit_mask[i]
            i_bit_mask[2*i  ] = y_bit_mask[i]
            
        cstate = 0
        for i in range(bit_length-1,-1,-1):
            b1 = 0
            if (i_bit_mask[2*i+1]):
                b1 = 1
            b0 = 0
            if (i_bit_mask[2*i  ]):
                b0 = 1
            sdigit = b1*2 + b0
            nstate = state_diagram[sdigit,0,cstate]
            hdigit = state_diagram[sdigit,1,cstate]
            i_bit_mask[2*i+1] = hdigit & (1 << 1)
            i_bit_mask[2*i  ] = hdigit & (1 << 0)
            cstate = nstate
            
        order[ip]= 0
        for i in range(0,2*bit_length):
            b0 = 0
            if (i_bit_mask[i]):
                b0 = 1
            order[ip] = order[ip] + float(b0)*2.**i
                
    return order

def get_cpu_list(info,**kwargs):

    center = kwargs.get("center")
    radius = kwargs.get("radius")
    center = np.array(center)
    radius = float(radius)
    
    for ilevel in range(0,info.nlevelmax):
        dx = 1/2**ilevel
        if (dx < 2*radius/info.boxlen):
            break

    levelmin = np.max([ilevel,1])
    bit_length = levelmin-1
    nmax = 2**bit_length
    ndim = info.ndim
    ncpu = info.ncpu
    nlevelmax = info.nlevelmax
    dkey = 2**(ndim*(nlevelmax+1-bit_length))
    ibound = [0, 0, 0, 0, 0, 0]
    if(bit_length > 0):
        ibound[0:3] = (center-radius)*nmax/info.boxlen
        ibound[3:6] = (center+radius)*nmax/info.boxlen
        ibound[0:3] = np.array(ibound[0:3]).astype(int)
        ibound[3:6] = np.array(ibound[3:6]).astype(int)
        ndom = 8
        idom = [ibound[0], ibound[3], ibound[0], ibound[3], ibound[0], ibound[3], ibound[0], ibound[3]]
        jdom = [ibound[1], ibound[1], ibound[4], ibound[4], ibound[1], ibound[1], ibound[4], ibound[4]]
        kdom = [ibound[2], ibound[2], ibound[2], ibound[2], ibound[5], ibound[5], ibound[5], ibound[5]]
        order_min = hilbert3d(idom,jdom,kdom,bit_length)
    else:
        ndom = 1
        order_min = np.array([0.])
        
    bounding_min = order_min*dkey
    bounding_max = (order_min+1)*dkey

    cpu_min = np.zeros(ndom, dtype=int)
    cpu_max = np.zeros(ndom, dtype=int)
    for icpu in range(0,ncpu):
        for idom in range(0,ndom):
            if( (info.bound_key[icpu] <= bounding_min[idom]) and (info.bound_key[icpu+1] > bounding_min[idom]) ):
                cpu_min[idom] = icpu+1
            if( (info.bound_key[icpu] < bounding_max[idom]) and (info.bound_key[icpu+1] >= bounding_max[idom]) ):
                cpu_max[idom] = icpu+1


    ncpu_read = 0
    cpu_read = np.zeros(ncpu, dtype=bool)
    cpu_list = []
    for idom in range(0,ndom):
        for icpu in range(cpu_min[idom]-1,cpu_max[idom]):
            if ( not cpu_read[icpu] ):
                cpu_list.append(icpu+1)
                ncpu_read = ncpu_read+1
                cpu_read[icpu] = True

    return cpu_list

def visu(x,y,dx,v,**kwargs):
    '''The simple visualization function visu() make a 2D scatter plot from RAMSES AMR data. 
    
    Args:
    
        x: the x-coordinate of the cells to show on the scatter plot.

        y: the y-coordinate of the cells to show on the scatter plot.
    
        dx: the size of the cells to show on the scatter plot.

        v: the value to show as a color square contained in the cell.
        
    Optional args:

        vmin: minimum value for the input array v to use in the color range

        vmax: maximum value for the input array v to use in the color range 

        log: when set, use the log of the input array v in the color range

        sort: useful only for 3D data. Plot the square symbola in the scatter plot in increasing order of array sort.

    Returns:
    
        Output a scatter plot figure of size 1000 pixels aside.
    
    Exemple:
    
        Example for a 2D or 3D RAMSES dataset using variable c from the object Cell. 
        import miniramses as ram
        c=ram.rd_cell(2)
        ram.visu(c.x[0],c.x[1],c.dx,c.u[0],sort=c.u[0],log=1,vmin=-3,vmax=1)

    Authors: Romain Teyssier (Princeton University, October 2022)
    '''

    xmin=np.min(x-dx/2)
    xmax=np.max(x+dx/2)
    ymin=np.min(y-dx/2)
    ymax=np.max(y+dx/2)
    
    log = kwargs.get("log",None)
    vmin = kwargs.get("vmin",None)
    vmax = kwargs.get("vmax",None)
    sort = kwargs.get("sort",None)
    cmap = kwargs.get("cmap",'viridis')

    if( not (log is None)):
        if vmin==None:
            v = np.log10(abs(v))
        else:
            v = np.log10(abs(v+float(vmin)))            
        
    print("min=",np.min(v)," max=",np.max(v))

    if( not (sort is None)):
        ind = np.argsort(sort)
    else:
        ind = np.arange(0,v.size)

    plt.rcParams['figure.dpi'] = 58
    px = 1/plt.rcParams['figure.dpi'] 
    fig, ax = plt.subplots(figsize=(1000*px,1000*px))
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.scatter(x,y,s=0.0001)
    rescale=np.maximum(xmax-xmin,ymax-ymin)        
    ax.set_aspect("equal")
    plt.scatter(x[ind],y[ind],c=v[ind],s=(dx[ind]*800/rescale)**2,marker="s",vmin=vmin,vmax=vmax,cmap=cmap)
    plt.colorbar(shrink=0.8)

def mk_movie(**kwargs):
    '''The function mk_movie() takes 2D data files containing maps and converts them into a sequence of images, 
    before combining them into a movie. It requires a standard set of python packages and the Linux packages
    ffmpeg and convert (ImageMagick).
    
    Args:
    
        start: starting index of the sequence of numpy array you wish to turn into image frames.
    
        stop: number of arrays you wish to be turned into plots. 
            This will be the variable "snum" for the end product. 
            For now, if you wish to test out the function, 
            you can try out other smaller values to adjust the image for your preferences.
            
        path: path leading to the directory where your files are stored, Default: "."
    
        prefix: starting name of a typical file. Ex: if you have 50 files, called "fig01.npy", "fig02.npy" … "fig50.npy", write in "fig".
    
        fill: This is for the zfill parameter. If your files are standardized into "fig001.npy", "fig002.npy"… "fig100.npy",
            write in 3, for example. If this is not how your files are formatted, write in the number 1.
    
        suffix: suffix at the end of a file: Ex: ".npy", ".map", etc…
    
        cmap: write in what color you wish your array to be displayed in (value for cmap). Options include "Reds", "Blues", and more.
    
        cbar: write "YES" for this parameter if you want your figure to have a colorbar. Write anything else if not.
    
        cbunit: units of the colormapping to be displayed next to the colorbar: Ex: "Concentration [code units]"
              If you do not plan on using a colorbar, write in any script.
    
        tunit: units of time displayed by rd_img. Ex: "seconds", "minutes", "hours", "[code units]"
    
        bunit: units of the box size displayed by rd_img. Ex: "cm", "kpc", "Mpc", "[code units]" …
    
        fname: starting name of each of your images.
    
        mvname: what you want your movie to be called.
        
    Returns:
    
        info: a string stating that the movie was done.
    
    Exemple:

        import miniramses as ram
        info = ram.mk_movie(start=100,stop=2000,path="../movie1",prefix="dens_",fill=5,suffix=".map",cmap="Reds", 
                cbar="YES", cbunit="log Density [H/cc]", tunit="Gyr",
                fname="img", mvname="movie", vmin=-1, vmax=6)
    
    By default, the movie's framerate is 30 frames per second, at a resolution of 420p
    You can edit this function and its parameters according to what fits your model best.
    
    As it runs, the function will print the files it is currently converting.

    Authors: Thomas Decugis and Romain Teyssier (Princeton University, October 2022)
    '''
    start = kwargs.get("start",1)
    stop = kwargs.get("stop",1)
    prefix = kwargs.get("prefix")
    suffix = kwargs.get("suffix")
    fill = kwargs.get("fill",5)
    path = kwargs.get("path",".")
    vmin = kwargs.get("vmin",None)
    vmax = kwargs.get("vmax",None)
    cmap = kwargs.get("cmap","Reds")
    cbar = kwargs.get("cbar",None)
    cbunit = kwargs.get("cbunit",None)
    tunit = kwargs.get("tunit"," ")
    bsize = kwargs.get("bsize",1)
    bunit = kwargs.get("bunit","[code units]")
    fname = kwargs.get("fname","frame")
    mvname = kwargs.get("mvname","movie")
    
    cmd="curl https://tigress-web.princeton.edu/~rt3504/DAT/logo_essai.jpg --output logo_essai.jpg"
    os.system(cmd)
    concom = "convert logo_essai.jpg -resize 280x200 logo_essai.png"
    os.system(concom)

    for snapshot in range(start, stop + 1): 
        ar = path + "/" + str(prefix) + str(snapshot).zfill(fill) + str(suffix)
        print(ar) #prints file that function is working on.

        map =rd_map(ar)
        time = map.time
        array = map.data
        
        if (not (cbar is None)):
            px = 1/plt.rcParams['figure.dpi']
            fig, ax = plt.subplots(figsize=(1000*px,1000*px))
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            print(np.min(array),np.max(array))
            shw = ax.imshow(array, cmap = cmap, vmin=vmin, vmax=vmax, origin="lower", extent=[0,bsize,0,bsize])
            bar = plt.colorbar(shw,shrink=0.8)
            bar.set_label(cbunit, fontsize=18) 
            bar.ax.tick_params(labelsize=18) 
            plt.ylabel(bunit,fontsize=18)

        else:
            plt.imshow(array, cmap = cmap)#if you wish to graph model in a specific way, modify this program

        ax = plt.gca()
        txt = f't = {time:4.2f}' + tunit
        label = ax.set_xlabel(txt, fontsize = 18, color = "black")
        ax.xaxis.set_label_coords(0.1, 0.95)
        ax.tick_params(axis='both', labelsize=18)
        newname = str(fname)+ str(snapshot).zfill(fill) + ".png"
        print(newname)
        plt.savefig(newname) #saves created images as pngs under the name that was given
        if snapshot == start:
            plt.show()
        plt.close(fig)
        com = "convert logo_essai.png -bordercolor white -border 0.1 " + newname + " +swap -geometry +100+850 -composite " + newname
        os.system(com)
    print("Input files converted into frames: done")
    moviecom = "ffmpeg -y -r 30 -f image2 -s 1000x1000 -start_number " +str(start)+" -i " + str(fname) + "%05d.png" + " -vcodec libx264 -crf 25  -pix_fmt yuv420p " + str(mvname) + ".mp4" 
    os.system(moviecom)
    ok = "Movie: done"
    print(ok)
    return ok 

class HaloCat:
   """
   This is the class for RAMSES halo catalogue.
   """
   def __init__(self):
       """
       This function initialize the halo catalogue.
       """
       self.index = np.empty(shape=(0),dtype=int)
       self.patch = np.empty(shape=(0),dtype=int)
       self.npart = np.empty(shape=(0),dtype=int)
       self.x = np.empty(shape=(0))
       self.y = np.empty(shape=(0))
       self.z = np.empty(shape=(0))
       self.u = np.empty(shape=(0))
       self.v = np.empty(shape=(0))
       self.w = np.empty(shape=(0))
       self.mpatch = np.empty(shape=(0))
       self.mass = np.empty(shape=(0))
       self.r200 = np.empty(shape=(0))
       self.rmax = np.empty(shape=(0))
       self.c200 = np.empty(shape=(0))

def rd_halo(nout,**kwargs):
   """
   This function reads and compiles data for position, mass,
   density, and index from the halo catalogue.
   Args:
       nout:output file number
   author: Josiah Taylor
   """
   backup = kwargs.get("backup",False)
   center = kwargs.get("center")
   radius = kwargs.get("radius")
   path = kwargs.get("path","./")

   car1 = str(nout).zfill(5)
   i = rd_info(nout,path=path,backup=backup)
   ncpu = i.ncpu
   ndim = i.ndim

   output = str(nout).zfill(5)
   cat = HaloCat()
   for i in range(0, ncpu):
       name = str(i+1).zfill(5)
       file_name = path+"/output_%s/halo.%s" % (output,name)
       halo_cat = ascii.read(file_name)
       index = halo_cat['index']
       patch = halo_cat['patch']
       npart = halo_cat['npart']
       x = halo_cat['pos_x']
       y = halo_cat['pos_y']
       z = halo_cat['pos_z']
       u = halo_cat['vel_x']
       v = halo_cat['vel_y']
       w = halo_cat['vel_z']
       mpatch = halo_cat['mpatch']
       mass = halo_cat['mass']
       r200 = halo_cat['r200']
       rmax = halo_cat['rmax']
       c200 = halo_cat['c200']
       cat.index = np.append(cat.index,index)
       cat.patch = np.append(cat.patch,patch)
       cat.npart = np.append(cat.npart,npart)
       cat.x = np.append(cat.x,x)
       cat.y = np.append(cat.y,y)
       cat.z = np.append(cat.z,z)
       cat.u = np.append(cat.u,u)
       cat.v = np.append(cat.v,v)
       cat.w = np.append(cat.w,w)
       cat.mpatch = np.append(cat.mpatch,mpatch)
       cat.mass = np.append(cat.mass,mass)
       cat.r200 = np.append(cat.r200,r200)
       cat.rmax = np.append(cat.rmax,rmax)
       cat.c200 = np.append(cat.c200,c200)

   return cat

class ClumpCat:
   """
   This is the class for RAMSES halo catalogue.
   """
   def __init__(self):
       """
       This function initialize the halo catalogue.
       """
       self.index = np.empty(shape=(0),dtype=int)
       self.halo = np.empty(shape=(0),dtype=int)
       self.ncell = np.empty(shape=(0),dtype=int)
       self.x = np.empty(shape=(0))
       self.y = np.empty(shape=(0))
       self.z = np.empty(shape=(0))
       self.u = np.empty(shape=(0))
       self.v = np.empty(shape=(0))
       self.w = np.empty(shape=(0))
       self.mass = np.empty(shape=(0))
       self.dmax = np.empty(shape=(0))
       self.dmin = np.empty(shape=(0))
       self.dsad = np.empty(shape=(0))

def rd_clump(nout,**kwargs):
   """
   This function reads and compiles data for position, mass,
   density, and index from the halo catalogue.
   Args:
       nout:output file number
   author: Josiah Taylor
   """
   backup = kwargs.get("backup",False)
   center = kwargs.get("center")
   radius = kwargs.get("radius")
   path = kwargs.get("path","./")

   car1 = str(nout).zfill(5)
   i = rd_info(nout,path=path,backup=backup)
   ncpu = i.ncpu
   ndim = i.ndim

   output = str(nout).zfill(5)
   cat = ClumpCat()
   for i in range(0, ncpu):
       name = str(i+1).zfill(5)
       file_name = path+"/output_%s/clump.%s" % (output,name)
       halo_cat = ascii.read(file_name)
       index = halo_cat['index']
       halo = halo_cat['halo']
       ncell = halo_cat['ncell']
       x = halo_cat['pos_x']
       y = halo_cat['pos_y']
       z = halo_cat['pos_z']
       u = halo_cat['vel_x']
       v = halo_cat['vel_y']
       w = halo_cat['vel_z']
       mass = halo_cat['mass']
       dmax = halo_cat['rho+']
       dmin = halo_cat['rho-']
       dsad = halo_cat['relevance']
       dsad = dmax/dsad
       cat.index = np.append(cat.index,index)
       cat.halo = np.append(cat.halo,halo)
       cat.ncell = np.append(cat.ncell,ncell)
       cat.x = np.append(cat.x,x)
       cat.y = np.append(cat.y,y)
       cat.z = np.append(cat.z,z)
       cat.u = np.append(cat.u,u)
       cat.v = np.append(cat.v,v)
       cat.w = np.append(cat.w,w)
       cat.mass = np.append(cat.mass,mass)
       cat.dmax = np.append(cat.dmax,dmax)
       cat.dmin = np.append(cat.dmin,dmin)
       cat.dsad = np.append(cat.dsad,dsad)

   if ( not (center is None)  and not (radius is None) ):
       # Filtering clumps
       r = np.sqrt((cat.x-center[0])**2+(cat.y-center[1])**2+(cat.z-center[2])**2)
       cat.index = cat.index[r < radius]
       cat.halo = cat.halo[r < radius]
       cat.ncell = cat.ncell[r < radius]
       cat.x  = cat.x[r < radius]
       cat.y  = cat.y[r < radius]
       cat.z  = cat.z[r < radius]
       cat.u  = cat.u[r < radius]
       cat.v  = cat.v[r < radius]
       cat.w  = cat.w[r < radius]
       cat.mass = cat.mass[r < radius]
       cat.dmax = cat.dmax[r < radius]
       cat.dmin = cat.dmin[r < radius]
       cat.dsad = cat.dsad[r < radius]
        
   return cat

class GraficFile:
    """
    Thid is the empty class for grafic files data
    """

def rd_grafic(filein):
    """This function reads a grafic file (unformatted Fortran binary)
    as produced by the MUSIC code.

    Args:
        filename: the complete path (including the name) of the input grafic file.

    Returns:
        A grafic (class GraficFile) object.

    Example:
        import miniramses as ram
        g = ram.rd_grafic("ic_deltab")
        plt.imshow(g.data[:,:,0],origin="lower")

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(filein, 'r') as f:
        recl = ["i4", "i4", "i4", "f4", "f4", "f4", "f4", "f4", "f4", "f4", "f4"] 
        n1, n2, n3, dx, x1, x2, x3, a, omega_m, omega_l, h0 = f.read_record(*recl)
        n1=int(n1[0])
        n2=int(n2[0])
        n3=int(n3[0])
        print("Reading file "+filein)
        print("Found array of size=",n1,n2,n3)
        dat = np.zeros((n1,n2,n3))
        for k in range(n3):
            plane = f.read_reals('f4')
            dat[:,:,k] = plane.reshape((n1,n2))

    out = GraficFile()
    out.n1=n1
    out.n2=n2
    out.n3=n3
    out.dx=dx[0]
    out.x1=x1[0]
    out.x2=x2[0]
    out.x3=x3[0]
    out.omega_m=omega_m[0]
    out.omega_l=omega_l[0]
    out.h0=h0[0]
    out.data=np.array(dat)

    return out

def wr_grafic(dat,header1,header2,fileout):
    """This function writes a grafic file (unformatted Fortran binary)
    which is the file format produced e.g. by the MUSIC code.

    Args:
        dat: a 3D numpy array of type "f4"

        header1: a 1D numpy array with 3 elements of type "i4". It should contain the 3 dimensions of the input array.

        header2: a 1D mumpy array with 8 elements of type "f4". It should contain dx, xoff1, xoff2, xoff3 and 4 additional constants,

        filename: the complete path (including the name) of the output grafic file.

    Returns:
        Nothing

    Example:
        import miniramses as ram
        dat = np.zeros((512,512,512),dtype="f4")
        dx = 1./512.
        header1 = np.array([512,512,512],dtype="i4")
        header2 = np.array([dx,0,0,0,0,0,0,0],dtype="f4")
        ram.wr_grafic(dat,header1,header2,"ic_d")

    Authors: Romain Teyssier (Princeton University, October 2022)
    """
    with FortranFile(fileout, 'w') as f:
        f.write_record(header1,header2)
        n3 = int(header1[2])
        for k in range(n3):
            plane = dat[:, :, k]
            f.write_record(plane.T)


