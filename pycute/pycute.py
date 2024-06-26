import os
import ctypes

import numpy as np
from numpy import ctypeslib

from . import utils


def wrap_phi(phi):
    mask = phi<0.
    phi[mask] = 2*np.pi+phi[mask]
    return phi


def celestial_to_spherical(pos,degree=True):
    toradians = np.pi/180. if degree else 1.
    ra,dec = pos.T
    theta = (90.-dec)*toradians
    phi = ra*toradians
    phi = wrap_phi(phi)
    return np.asarray([theta,phi]).T


class PyCute(object):

    C_TYPE = ctypes.c_double
    PATH_CUTE = os.path.join(utils.lib_dir,'cute.so')

    def __init__(self):

        self.cute = ctypes.CDLL(self.PATH_CUTE,mode=ctypes.RTLD_LOCAL)
        self.clear()

    def set_verbosity(self,mode='info'):
        self.cute.set_verbosity.argtypes = (ctypes.c_char_p,)
        self.cute.set_verbosity(mode.encode('utf-8'))
        self._verbose = mode

    def set_bin(self,mode,edges=None,size=None,binning='lin'):

        if mode == 'bin':
            assert size is not None
            binning = 'bin'
            edges = np.arange(size+1,dtype=self.C_TYPE)
        elif binning == 'lin':
            if size is None: size = len(edges)-1
            edges = np.linspace(edges[0],edges[-1],size+1,dtype=self.C_TYPE)
        elif binning == 'log':
            if size is None: size = len(edges)-1
            edges = np.logspace(np.log10(edges[0]),np.log10(edges[-1]),size+1,base=10,dtype=self.C_TYPE)
        elif binning == 'custom':
            size = len(edges) - 1
            edges = np.asarray(edges,dtype=self.C_TYPE).flatten()
        else:
            raise ValueError('Wrong binning type/mode.')

        self.cute.set_bin.argtypes = (ctypes.c_char_p,ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(size+1,)),ctypes.c_size_t,ctypes.c_char_p)
        self.cute.set_bin(mode.encode('utf-8'),edges,size,binning.encode('utf-8'))
        self._edges[mode] = edges

        return self._edges[mode].copy()

    def set_pole(self,num=1,ells=[0,2,4,6,8,10,12]):

        nells = len(ells)
        if (np.mod(ells,2)==0).all():
            ells = np.arange(nells)*2
            multitype = 'even'
        elif (np.mod(ells,2)==1).all():
            ells = 1+np.arange(nells)*2
            multitype = 'odd'
        else:
            ells = np.arange(nells)
            multitype = 'all'

        self.cute.set_pole.argtypes = (ctypes.c_size_t,ctypes.c_char_p,ctypes.c_size_t)
        self.cute.set_pole(num,multitype.encode('utf-8'),nells)
        self._ells[num] = ells.tolist()

        return list(self._ells[num])

    def set_los(self,num=1,los='midpoint',n=0):

        vec = np.array([0],dtype=self.C_TYPE)
        if not isinstance(los,str):
            vec = np.array(los,dtype=self.C_TYPE)
            los = 'custom'

        self.cute.set_los.argtypes = (ctypes.c_size_t,ctypes.c_char_p,ctypeslib.ndpointer(dtype=self.C_TYPE),ctypes.c_size_t)
        self.cute.set_los(num,los.encode('utf-8'),vec,n)

        self._los[num] = (los if los != 'custom' else vec,n)

        return tuple(self._los[num])

    def set_n_poles(self,ells,n=1):
        if np.isscalar(ells[0]): ells = [ells]*n
        ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in range(n)]
        return ells

    def set_n_los(self,los,losn=0,n=1):
        if isinstance(los,str) or (np.isscalar(los[0]) and not isinstance(los[0],str)): los = [los]*n
        if np.isscalar(losn): losn = [losn]*n
        los = [self.set_los(ilos+1,los=los[ilos],n=losn[ilos]) for ilos in range(n)]
        return los

    def set_tobin(self,bins,tobin=None,maxtobin=1):
        if tobin is None:
            tobin = [ibin+1 for ibin,bin in enumerate(bins) if bin is not None]
        else:
            if np.isscalar(tobin): tobin = [tobin]
            for tobin_ in tobin: assert (tobin_ >= 1) and (tobin_ <= len(bins)) and (bins[tobin_-1] is not None)
        assert (len(tobin) <= maxtobin)
        if len(tobin) == 1: return tobin[0]
        return tobin

    def set_weight_type(self,weighttype='prod'):
        self.cute.set_weight_type.argtypes = (ctypes.c_char_p,)
        self.cute.set_weight_type(weighttype.encode('utf-8'))
        return weighttype

    def set_2pcf_smu(self,sedges,muedges,position1,weight1,position2=None,weight2=None,sbinning='lin',mubinning='lin',ssize=None,musize=None,los='midpoint',losn=0,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=musize,binning=mubinning)
        self.los = self.set_los(1,los=los,n=losn)
        cross = self.set_catalogues([position1,position2],[weight1,weight2])

        self.run_2pcf_smu(nthreads=nthreads)

    def run_2pcf_smu(self,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._edges['aux'])-1)

        self.s = np.zeros(shape,dtype=self.C_TYPE).flatten()
        self.mu = np.zeros(shape,dtype=self.C_TYPE).flatten()
        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_main_aux.argtypes = (typecounts,typecounts,typecounts,ctypes.c_char_p,ctypes.c_size_t)
        self.cute.run_2pcf_main_aux(self.s,self.mu,self.counts,'s-mu'.encode('utf-8'),nthreads)
        self.s.shape = shape
        self.mu.shape = shape
        self.counts.shape = shape

    def set_2pcf_s(self,sedges,position1,weight1,position2=None,weight2=None,sbinning='lin',ssize=None,weighttype='prod',nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.weighttype = self.set_weight_type(weighttype)
        cross = self.set_catalogues([position1,position2],[weight1,weight2])

        self.run_2pcf_s(nthreads=nthreads)

    def run_2pcf_s(self,nthreads=8):

        shape = (len(self._edges['main'])-1)

        self.s = np.zeros(shape,dtype=self.C_TYPE).flatten()
        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_main.argtypes = (typecounts,typecounts,ctypes.c_char_p,ctypes.c_size_t)
        self.cute.run_2pcf_main(self.s,self.counts,'s-mu'.encode('utf-8'),nthreads)
        self.s.shape = shape
        self.counts.shape = shape

    def set_2pcf_angular(self,thetaedges,position1,weight1,position2=None,weight2=None,thetabinning='lin',thetasize=None,celestial=True,degree=True,nthreads=8):

        #theta edges in radians
        toradians = np.pi/180. if degree else 1.

        if celestial:
            position1 = celestial_to_spherical(position1,degree=degree)
            if position2 is not None: position2 = celestial_to_spherical(position2,degree=degree)
        else:
            position1 = position1.copy()*toradians
            position1[:,1] = wrap_phi(position1[:,1])
            if position2 is not None:
                position2 = position2.copy()*toradians
                position2[:,1] = wrap_phi(position2[:,1])

        self.thetaedges = self.set_bin('main',edges=np.asarray(thetaedges)*toradians,size=thetasize,binning=thetabinning)/toradians
        cross = self.set_catalogues([position1,position2],[weight1,weight2])

        self.run_2pcf_angular(degree=degree,nthreads=nthreads)

    def run_2pcf_angular(self,degree=True,nthreads=8):

        shape = (len(self._edges['main'])-1)

        self.theta = np.zeros(shape,dtype=self.C_TYPE).flatten()
        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_main.argtypes = (typecounts,typecounts,ctypes.c_char_p,ctypes.c_size_t)
        self.cute.run_2pcf_main(self.theta,self.counts,'angular'.encode('utf-8'),nthreads)
        self.theta.shape = shape
        if degree: self.theta /= np.pi/180.
        self.counts.shape = shape

    def set_2pcf_multi(self,sedges,position1,weight1,position2=None,weight2=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',losn=0,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=1)
        self.ells = self.set_pole(1,ells=ells)
        self.los = self.set_los(1,los=los,n=losn)
        cross = self.set_catalogues([position1,position2],[weight1,weight2])

        self.run_2pcf_multi(nthreads=nthreads)

    def run_2pcf_multi(self,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._ells[1]))

        self.s = np.zeros(shape,dtype=self.C_TYPE).flatten()
        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_multi.argtypes = (typecounts,typecounts,ctypes.c_size_t)
        self.cute.run_2pcf_multi(self.s,self.counts,nthreads)
        self.s.shape = shape
        self.counts.shape = shape

    def set_2pcf_scos(self,sedges,muedges,position1,weight1,position2=None,weight2=None,sbinning='lin',mubinning='lin',ssize=None,musize=None,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=musize,binning=mubinning)
        cross = self.set_catalogues([position1,position2],[weight1,weight2])

        self.run_2pcf_scos(nthreads=nthreads)

    def run_2pcf_scos(self,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._edges['aux'])-1)

        self.s = np.zeros(shape,dtype=self.C_TYPE).flatten()
        self.mu = np.zeros(shape,dtype=self.C_TYPE).flatten()
        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_main_aux.argtypes = (typecounts,typecounts,typecounts,ctypes.c_char_p,ctypes.c_size_t)
        self.cute.run_2pcf_main_aux(self.s,self.mu,self.counts,'s-cos'.encode('utf-8'),nthreads)
        self.s.shape = shape
        self.mu.shape = shape
        self.counts.shape = shape

    def set_3pcf_multi(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',losn=0,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=1)
        self.ells = self.set_n_poles(ells,n=2)
        self.los = self.set_n_los(los,losn,n=2)
        if (self.ells[-1] != self.ells[0]) and (position3 is None):
            position3 = position2
            weight3 = weight2

        self.set_catalogues([position1,position2,position3],[weight1,weight2,weight3])
        self.run_3pcf_multi(nthreads=nthreads)

    def run_3pcf_multi(self,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._edges['main'])-1,len(self._ells[1]),len(self._ells[2]))

        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_3pcf_multi.argtypes = (typecounts,ctypes.c_size_t)
        self.cute.run_3pcf_multi(self.counts,nthreads)
        self.counts.shape = shape

    def set_3pcf_multi_double_los(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',losn=0,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=1)
        self.ells = self.set_n_poles(ells,n=2)
        self.los = self.set_n_los(los,losn,n=2)
        if (self.ells[-1] != self.ells[0]) and (position3 is None):
            position3 = position2
            weight3 = weight2

        self.set_catalogues([position1,position2,position3],[weight1,weight2,weight3])
        self.run_3pcf_multi_double_los(nthreads=nthreads)

    def run_3pcf_multi_double_los(self,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._edges['main'])-1,len(self._ells[1]),len(self._ells[2]))

        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_3pcf_multi_double_los.argtypes = (typecounts,ctypes.c_size_t)
        self.cute.run_3pcf_multi_double_los(self.counts,nthreads)
        self.counts.shape = shape

    def set_2pcf_multi_binned(self,sedges,binsize,position1,weight1=None,bin1=None,position2=None,weight2=None,bin2=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',losn=0,tobin=None,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=1)
        self.binedges = self.set_bin('bin',size=binsize)
        self.ells = self.set_pole(1,ells=ells)
        self.los = self.set_los(1,los=los,n=losn)
        self.tobin = self.set_tobin([bin1,bin2],tobin=tobin,maxtobin=1)

        self.set_catalogues([position1,position2],[weight1,weight2],[bin1,bin2])
        self.run_2pcf_multi_binned(tobin=self.tobin,nthreads=nthreads)

    def run_2pcf_multi_binned(self,tobin=2,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self.binedges)-1,len(self._ells[1]))

        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_multi_binned.argtypes = (typecounts,ctypes.c_size_t,ctypes.c_size_t)
        self.cute.run_2pcf_multi_binned(self.counts,tobin,nthreads)
        self.counts.shape = shape

    def set_4pcf_multi_binned(self,sedges,binsize,position1,weight1=None,bin1=None,position2=None,weight2=None,bin2=None,position3=None,weight3=None,bin3=None,position4=None,weight4=None,bin4=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',losn=0,tobin=None,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=1)
        self.binedges = self.set_bin('bin',size=binsize)
        self.ells = self.set_n_poles(ells,n=2)
        self.los = self.set_n_los(los,losn,n=2)
        self.tobin = self.set_tobin([bin1,bin2,bin3,bin4],tobin=tobin,maxtobin=2)

        self.set_catalogues([position1,position2,position3,position4],[weight1,weight2,weight3,weight4],[bin1,bin2,bin3,bin4])
        self.run_4pcf_multi_binned(tobin=self.tobin,nthreads=nthreads)

    def run_4pcf_multi_binned(self,tobin=[2,2],nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._edges['main'])-1,len(self._ells[1]),len(self._ells[2]))

        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
        typetobin = ctypeslib.ndpointer(ctypes.c_size_t,shape=(2))

        self.cute.run_4pcf_multi_binned.argtypes = (typecounts,typetobin,ctypes.c_size_t)
        self.cute.run_4pcf_multi_binned(self.counts,np.array(tobin,dtype=ctypes.c_size_t),nthreads)
        self.counts.shape = shape

    def set_2pcf_multi_radial_legendre(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',losn=0,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=1)
        self.ells = self.set_n_poles(ells,n=2)
        self.los = self.set_los(1,los=los,n=losn)

        self.set_catalogues([position1,position2],[weight1,weight2])
        self.run_2pcf_multi_radial_legendre(nthreads=nthreads)

    def run_2pcf_multi_radial_legendre(self,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._edges['main'])-1,len(self._ells[1]),len(self._ells[2]))

        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_multi_radial_legendre.argtypes = (typecounts,ctypes.c_size_t)
        self.cute.run_2pcf_multi_radial_legendre(self.counts,nthreads)
        self.counts.shape = shape

    def set_2pcf_multi_angular_legendre(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',losn=0,nthreads=8):

        self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
        self.muedges = self.set_bin('aux',edges=muedges,size=1)
        self.ells = self.set_n_poles(ells,n=2)
        self.los = self.set_los(1,los=los,n=losn)

        self.set_catalogues([position1,position2],[weight1,weight2])
        self.run_2pcf_multi_angular_legendre(nthreads=nthreads)

    def run_2pcf_multi_angular_legendre(self,nthreads=8):

        shape = (len(self._edges['main'])-1,len(self._edges['main'])-1,len(self._ells[1]),len(self._ells[2]))

        self.counts = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))

        self.cute.run_2pcf_multi_angular_legendre.argtypes = (typecounts,ctypes.c_size_t)
        self.cute.run_2pcf_multi_angular_legendre(self.counts,nthreads)
        self.counts.shape = shape

    def set_catalogues(self,positions,weights,bins=None):

        self.cute.clear_catalogs()

        def size(tab):
            return tab.shape[-1] if len(tab.shape) > 1 else 1

        sizeposition = size(positions[0])
        sizeweight = size(weights[0]) if weights[0] is not None else 1
        if bins is None: bins = [None]*len(positions)
        for num,(position,weight,bin) in enumerate(zip(positions,weights,bins)):
            if position is not None:
                if size(position) != sizeposition: raise ValueError('You have to work with the same dimensionality')
                if weight is None:
                    weight = np.ones((position.shape[0],sizeweight),dtype=position.dtype)
                else:
                    assert size(weight) == sizeweight
                if bin is None:
                    bin = np.zeros((position.shape[0],1),dtype=np.int8)
                self.set_catalogue(num+1,position,weight,bin)

    def set_catalogue(self,num,position,weight,bin,copy=False):

        #First dimension is number of objects
        n = position.shape[0]

        def size(tab):
            return tab.shape[-1] if len(tab.shape) > 1 else 1

        sizeposition = size(position)
        sizeweight = size(weight)
        assert size(bin) == 1

        def new(tab,dtype=self.C_TYPE):
            if copy: return np.asarray(tab,dtype=dtype).flatten()
            return tab.astype(dtype).flatten()

        self._position[num] = new(position)
        self._weight[num] = new(weight)
        self._bin[num] = new(bin,dtype=ctypes.c_size_t)

        typeposition = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(n*sizeposition))
        typeweight = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(n*sizeweight))
        typebin = ctypeslib.ndpointer(dtype=ctypes.c_size_t,shape=(n))

        self.cute.set_catalog.argtypes = (ctypes.c_size_t,typeposition,typeweight,typebin,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_size_t)
        self.cute.set_catalog(num,self._position[num],self._weight[num],self._bin[num],n,sizeposition,sizeweight)

        self._position[num].shape = (n,sizeposition)
        self._weight[num].shape = (n,sizeweight)
        self._bin[num].shape = (n,1)

    def clear(self):
        self._position = {}
        self._weight = {}
        self._bin = {}
        self._edges = {}
        self._ells = {}
        self._los = {}
        self.cute.clear_catalogs()
        self.cute.clear_bins()
        self.cute.clear_poles()
        self.set_verbosity('quiet')
        self.set_weight_type()
        self.set_verbosity()

    def integrate_legendre(self,ells=None,nthreads=8):

        if ells is not None: self.ells = self.set_pole(ells=ells)

        shape = (len(self._edges['main'])-1,len(self._ells[1]))

        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(self.counts.size))
        self.integral = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typeintegral = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.integral)))

        self.cute.integrate_legendre.argtypes = (typecounts,typeintegral,ctypes.c_size_t)
        self.cute.integrate_legendre(self.counts.flatten(),self.integral,nthreads)

        self.integral.shape = shape

    def integrate_radial_legendre(self,ells=None,nthreads=8):

        if ells is not None: self.ells = self.set_n_poles(ells,n=2)

        shape = (len(self._edges['main'])-1,len(self._edges['main'])-1,len(self._ells[1]),len(self._ells[2]))

        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(self.counts.size))
        self.integral = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typeintegral = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.integral)))

        self.cute.integrate_radial_legendre.argtypes = (typecounts,typeintegral,ctypes.c_size_t)
        self.cute.integrate_radial_legendre(self.counts.flatten(),self.integral,nthreads)

        self.integral.shape = shape

    def integrate_angular_legendre(self,ells=None,nthreads=8):

        if ells is not None: self.ells = self.set_n_poles(ells,n=2)

        shape = (len(self._edges['main'])-1,len(self._edges['main'])-1,len(self._ells[1]),len(self._ells[2]))

        typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(self.counts.size))
        self.integral = np.zeros(shape,dtype=self.C_TYPE).flatten()
        typeintegral = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.integral)))

        self.cute.integrate_angular_legendre.argtypes = (typecounts,typeintegral,ctypes.c_size_t)
        self.cute.integrate_angular_legendre(self.counts.flatten(),self.integral,nthreads)

        self.integral.shape = shape
