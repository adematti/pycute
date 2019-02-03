# coding: utf8

import os
import ctypes
import scipy
from scipy import constants
from numpy import ctypeslib

def wrap_phi(phi):
	mask = phi<0.
	phi[mask] = 2*constants.pi+phi[mask]
	return phi

def celestial_to_spherical(pos):
	ra,dec = pos.T
	theta = (90.-dec)*constants.degree
	phi = ra*constants.degree
	phi = wrap_phi(phi)
	return scipy.asarray([theta,phi]).T

class PyCute(object):

	C_TYPE = ctypes.c_double
	PATH_CUTE = os.path.join(os.path.dirname(os.path.realpath(__file__)),'cute.so')

	def __init__(self):

		self.cute = ctypes.CDLL(self.PATH_CUTE,mode=ctypes.RTLD_GLOBAL)
		self.clear()

	def set_bin(self,mode,edges=None,size=None,binning='lin'):

		if mode == 'bin':
			assert size is not None
			binning = 'bin'
			edges = scipy.arange(size+1,dtype=self.C_TYPE)
		elif binning == 'lin':
			if size is None: size = len(edges)-1
			edges = scipy.linspace(edges[0],edges[-1],size+1,dtype=self.C_TYPE)
		elif binning == 'log':
			if size is None: size = len(edges)-1
			edges = scipy.logspace(scipy.log10(edges[0]),scipy.log10(edges[-1]),size+1,base=10,dtype=self.C_TYPE)
		elif binning == 'custom':
			size = len(edges) - 1
			edges = scipy.asarray(edges,dtype=self.C_TYPE).flatten()
		else:
			raise ValueError('Wrong binning type/mode.')
	
		self.cute.set_bin.argtypes = (ctypes.c_char_p,ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(size+1,)),ctypes.c_size_t,ctypes.c_char_p)
		self.cute.set_bin(mode,edges,size,binning)
	
		return edges
	
	def set_pole(self,num=1,ells=[0,2,4,6,8,10,12]):

		nells = len(ells)
		if (scipy.mod(ells,2)==0).all():
			ells = scipy.arange(nells)*2
			multitype = 'even'
		elif (scipy.mod(ells,2)==1).all():
			ells = 1+scipy.arange(nells)*2
			multitype = 'odd'
		else:
			ells = scipy.arange(nells)
			multitype = 'all'
		
		self.cute.set_pole.argtypes = (ctypes.c_size_t,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.set_pole(num,multitype,nells)

		return ells.tolist()

	def set_2pcf_smu(self,sedges,muedges,position1,weight1,position2=None,weight2=None,sbinning='lin',mubinning='lin',ssize=None,musize=None,los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=musize,binning=mubinning)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_smu(los=los,nthreads=nthreads)

	def run_2pcf_smu(self,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.muedges)-1)
	
		self.s = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.mu = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_main_aux.argtypes = (typecounts,typecounts,typecounts,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_main_aux(self.s,self.mu,self.counts,'s-mu',los,nthreads)
		self.s.shape = shape
		self.mu.shape = shape
		self.counts.shape = shape
	
	def set_2pcf_s(self,sedges,position1,weight1,position2=None,weight2=None,sbinning='lin',ssize=None,nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_s(nthreads=nthreads)

	def run_2pcf_s(self,nthreads=8):
	
		shape = (len(self.sedges)-1)
	
		self.s = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_main.argtypes = (typecounts,typecounts,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_main(self.s,self.counts,'s-mu',nthreads)
		self.s.shape = shape
		self.counts.shape = shape

	def set_2pcf_angular(self,thetaedges,position1,weight1,position2=None,weight2=None,thetabinning='lin',thetasize=None,celestial=True,nthreads=8):

		#theta edges in radians
		toradians = constants.degree if celestial else 1.

		if celestial:
			position1 = celestial_to_spherical(position1)
			if position2 is not None: position2 = celestial_to_spherical(position2)
		else:
			position1[:,1] = wrap_phi(position1[:,1])
			if position2 is not None: position2[:,1] = wrap_phi(position2[:,1])
		
		self.thetaedges = self.set_bin('main',edges=scipy.asarray(thetaedges)*toradians,size=thetasize,binning=thetabinning)/toradians
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
	
		self.run_2pcf_angular(degree=True,nthreads=nthreads)

	def run_2pcf_angular(self,degree=True,nthreads=8):
		
		shape = (len(self.thetaedges)-1)
	
		self.theta = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_main.argtypes = (typecounts,typecounts,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_main(self.theta,self.counts,'angular',nthreads)
		self.theta.shape = shape
		if degree: self.theta /= constants.degree
		self.counts.shape = shape
		
	def set_2pcf_multi(self,sedges,position1,weight1,position2=None,weight2=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=1)
		self.ells = self.set_pole(1,ells=ells)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_multi(los=los,nthreads=nthreads)

	def run_2pcf_multi(self,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.ells))
	
		self.s = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_multi.argtypes = (typecounts,typecounts,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_multi(self.s,self.counts,los,nthreads)
		self.s.shape = shape
		self.counts.shape = shape
		
	def set_2pcf_scos(self,sedges,muedges,position1,weight1,position2=None,weight2=None,sbinning='lin',mubinning='lin',ssize=None,musize=None,nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=musize,binning=mubinning)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_scos(nthreads=nthreads)

	def run_2pcf_scos(self,nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.muedges)-1)
	
		self.s = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.mu = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_main_aux.argtypes = (typecounts,typecounts,typecounts,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_main_aux(self.s,self.mu,self.counts,'s-cos','midpoint',nthreads)
		self.s.shape = shape
		self.mu.shape = shape
		self.counts.shape = shape
			
	def set_3pcf_multi(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=1)
		if scipy.isscalar(ells[0]): ells = [ells,ells]
		self.ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in [0,1]]
		if (self.ells[-1] != self.ells[0]) and (position3 is None):
			position3 = position2
			weight3 = weight2
		
		self.set_catalogues([position1,position2,position3],[weight1,weight2,weight3])
		self.run_3pcf_multi(los=los,nthreads=nthreads)

	def run_3pcf_multi(self,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,len(self.ells[0]),len(self.ells[1]))
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_3pcf_multi.argtypes = (typecounts,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_3pcf_multi(self.counts,los,nthreads)
		self.counts.shape = shape

	def set_3pcf_multi_double_los(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='endpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=1)
		if scipy.isscalar(ells[0]): ells = [ells,ells]
		self.ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in [0,1]]
		if (self.ells[-1] != self.ells[0]) and (position3 is None):
			position3 = position2
			weight3 = weight2
		
		self.set_catalogues([position1,position2,position3],[weight1,weight2,weight3])
		self.run_3pcf_multi_double_los(los=los,nthreads=nthreads)

	def run_3pcf_multi_double_los(self,los='endpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,len(self.ells[0]),len(self.ells[1]))
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_3pcf_multi_double_los.argtypes = (typecounts,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_3pcf_multi_double_los(self.counts,los,nthreads)
		self.counts.shape = shape
		
	def set_2pcf_multi_binned(self,sedges,binsize,position1,weight1,position2,weight2,bin2,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],normalize=False,los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=1)
		self.binedges = self.set_bin('bin',size=binsize)
		self.ells = self.set_pole(1,ells=ells)
		
		self.set_catalogues([position1,position2],[weight1,weight2],[None,bin2])
		self.run_2pcf_multi_binned(normalize=normalize,los=los,nthreads=nthreads)

	def run_2pcf_multi_binned(self,normalize=False,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.binedges)-1,len(self.ells))
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_multi_binned.argtypes = (typecounts,ctypes.c_bool,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_multi_binned(self.counts,normalize,los,nthreads)
		self.counts.shape = shape
		
	def set_4pcf_multi_binned(self,sedges,binsize,position1,weight1,position2,weight2,bin2,position3,weight3,position4,weight4,bin4,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],normalize=False,los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=1)
		self.binedges = self.set_bin('bin',size=binsize)
		if scipy.isscalar(ells[0]): ells = [ells,ells]
		self.ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in [0,1]]
		
		self.set_catalogues([position1,position2,position3,position4],[weight1,weight2,weight3,weight4],[None,bin2,None,bin4])
		self.run_4pcf_multi_binned(normalize=normalize,los=los,nthreads=nthreads)

	def run_4pcf_multi_binned(self,normalize=False,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,len(self.ells[0]),len(self.ells[1]))
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_4pcf_multi_binned.argtypes = (typecounts,ctypes.c_bool,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_4pcf_multi_binned(self.counts,normalize,los,nthreads)
		self.counts.shape = shape

	def set_2pcf_multi_radial_legendre(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=1)
		if scipy.isscalar(ells[0]): ells = [ells,ells]
		self.ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in [0,1]]
		self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_multi_radial_legendre(los=los,nthreads=nthreads)

	def run_2pcf_multi_radial_legendre(self,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,len(self.ells[0]),len(self.ells[1]))
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_multi_radial_legendre.argtypes = (typecounts,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_multi_radial_legendre(self.counts,los,nthreads)
		self.counts.shape = shape

	def set_2pcf_multi_angular_legendre(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',edges=sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',edges=muedges,size=1)
		if scipy.isscalar(ells[0]): ells = [ells,ells]
		self.ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in [0,1]]
		self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_multi_angular_legendre(los=los,nthreads=nthreads)

	def run_2pcf_multi_angular_legendre(self,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,len(self.ells[0]),len(self.ells[1]))
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_multi_angular_legendre.argtypes = (typecounts,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_multi_angular_legendre(self.counts,los,nthreads)
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
				if size(position) != sizeposition: raise ValueError
				if weight is None:
					weight = scipy.ones((position.shape[0],sizeweight),dtype=position.dtype)
				else:
					assert size(weight) == sizeweight
				if bin is None:
					bin = scipy.ones((position.shape[0],1),dtype=scipy.int8)
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
			if copy: return scipy.asarray(tab,dtype=dtype).flatten()
			return tab.astype(dtype).flatten()
		
		self.position[num] = new(position)
		self.weight[num] = new(weight)
		self.bin[num] = new(bin,dtype=ctypes.c_size_t)

		typeposition = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(n*sizeposition))
		typeweight = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(n*sizeweight))
		typebin = ctypeslib.ndpointer(dtype=ctypes.c_size_t,shape=(n))
		
		self.cute.set_catalog.argtypes = (ctypes.c_size_t,typeposition,typeweight,typebin,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_size_t)
		self.cute.set_catalog(num,self.position[num],self.weight[num],self.bin[num],n,sizeposition,sizeweight)
		
		self.position[num].shape = (n,sizeposition)
		self.weight[num].shape = (n,sizeweight)
		self.bin[num].shape = (n,1)

	def clear(self):
		self.position = {}
		self.weight = {}
		self.bin = {}
		self.cute.clear_catalogs()
		self.cute.clear_bins()
		self.cute.clear_poles()
	
	def integrate_legendre(self,ells=None,nthreads=8):
	
		if ells is not None: self.ells = self.set_pole(ells=ells)

		shape = (len(self.sedges)-1,len(self.ells))
		
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(self.counts.size))
		self.integral = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typeintegral = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.integral)))
	
		self.cute.integrate_legendre.argtypes = (typecounts,typeintegral,ctypes.c_size_t)
		self.cute.integrate_legendre(self.counts.flatten(),self.integral,nthreads)
		
		self.integral.shape = shape
	
	def integrate_radial_legendre(self,ells=None,nthreads=8):
	
		if ells is not None:
			if scipy.isscalar(ells[0]): ells = [ells,ells]
			self.ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in [0,1]]

		shape = (len(self.sedges)-1,len(self.sedges)-1,len(self.ells[0]),len(self.ells[1]))
		
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(self.counts.size))
		self.integral = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typeintegral = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.integral)))
	
		self.cute.integrate_radial_legendre.argtypes = (typecounts,typeintegral,ctypes.c_size_t)
		self.cute.integrate_radial_legendre(self.counts.flatten(),self.integral,nthreads)
		
		self.integral.shape = shape
	
	def integrate_angular_legendre(self,ells=None,nthreads=8):
	
		if ells is not None:
			if scipy.isscalar(ells[0]): ells = [ells,ells]
			self.ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in [0,1]]

		shape = (len(self.sedges)-1,len(self.sedges)-1,len(self.ells[0]),len(self.ells[1]))
		
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(self.counts.size))
		self.integral = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typeintegral = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.integral)))
	
		self.cute.integrate_angular_legendre.argtypes = (typecounts,typeintegral,ctypes.c_size_t)
		self.cute.integrate_angular_legendre(self.counts.flatten(),self.integral,nthreads)
		
		self.integral.shape = shape
		
