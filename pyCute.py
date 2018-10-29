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
		self.free()

	def set_bin(self,mode,edges,size=None,binning='lin'):

		if binning == 'lin':
			if size is None: size = len(edges)-1
			edges = scipy.linspace(edges[0],edges[-1],size+1,dtype=self.C_TYPE)
		elif binning == 'log':
			if size is None: size = len(edges)-1
			edges = scipy.logspace(scipy.log10(edges[0]),scipy.log10(edges[-1]),size+1,base=10,dtype=self.C_TYPE)
		else:
			edges = scipy.asarray(edges,dtype=self.C_TYPE)
			binning = 'custom'
	
		self.cute.set_bin.argtypes = (ctypes.c_char_p,ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(size+1,)),ctypes.c_size_t,ctypes.c_char_p)
		self.cute.set_bin(mode,edges,size,binning)
	
		return edges
		
	def set_ells(self,ells=[0,2,4,6,8,10,12]):
		self.nells = len(ells)
		if (scipy.mod(ells,2)==0).all():
			self.ells = scipy.arange(self.nells)*2
			self.multitype = 'even'
		elif (scipy.mod(ells,2)==1).all():
			self.ells = 1+scipy.arange(self.nells)*2
			self.multitype = 'odd'
		else:
			self.ells = scipy.arange(self.nells)
			self.multitype = 'all'

	def set_2pcf_smu(self,sedges,muedges,position1,weight1,position2=None,weight2=None,sbinning='lin',mubinning='lin',ssize=None,musize=None,los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',muedges,size=musize,binning=mubinning)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_smu(los=los,nthreads=nthreads)
		self.free()
	
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

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_s(nthreads=nthreads)
		self.free()
	
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
		
		self.thetaedges = self.set_bin('main',thetaedges*toradians,size=thetasize,binning=thetabinning)/toradians
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
	
		self.run_2pcf_angular(degree=True,nthreads=nthreads)
		self.free()
	
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

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',muedges,size=1)
		self.set_ells(ells=ells)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_multi(los=los,nthreads=nthreads)
		self.free()
	
	def run_2pcf_multi(self,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,self.nells)
	
		self.s = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_multi.argtypes = (typecounts,typecounts,ctypes.c_size_t,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_multi(self.s,self.counts,self.nells,self.multitype,los,nthreads)
		self.s.shape = shape
		self.counts.shape = shape
		
	def set_2pcf_scos(self,sedges,muedges,position1,weight1,position2=None,weight2=None,sbinning='lin',mubinning='lin',ssize=None,musize=None,nthreads=8):

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',muedges,size=musize,binning=mubinning)
		cross = self.set_catalogues([position1,position2],[weight1,weight2])
		
		self.run_2pcf_scos(nthreads=nthreads)
		self.free()
	
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

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',muedges,size=1)
		self.set_ells(ells=ells)
		
		self.set_catalogues([position1,position2,position3],[weight1,weight2,weight3])
		self.run_3pcf_multi(los=los,nthreads=nthreads)
		self.free()
	
	def run_3pcf_multi(self,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,self.nells,self.nells)
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_3pcf_multi.argtypes = (typecounts,ctypes.c_size_t,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_3pcf_multi(self.counts,self.nells,self.multitype,los,nthreads)
		self.counts.shape = shape

	def set_3pcf_multi_double_los(self,sedges,position1,weight1,position2,weight2,position3=None,weight3=None,sbinning='lin',ssize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],los='endpoint',nthreads=8):

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',muedges,size=1)
		self.set_ells(ells=ells)
		
		self.set_catalogues([position1,position2,position3],[weight1,weight2,weight3])
		self.run_3pcf_multi_double_los(los=los,nthreads=nthreads)
		self.free()
	
	def run_3pcf_multi_double_los(self,los='endpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,self.nells,self.nells)
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_3pcf_multi_double_los.argtypes = (typecounts,ctypes.c_size_t,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_3pcf_multi_double_los(self.counts,self.nells,self.multitype,los,nthreads)
		self.counts.shape = shape
		
	def set_2pcf_multi_radial(self,sedges,radialedges,position1,weight1,position2,weight2,sbinning='lin',ssize=None,radialbinning='lin',radialsize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],normalize=False,los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',muedges,size=1)
		self.radialedges = self.set_bin('radial',radialedges,size=radialsize,binning=radialbinning)
		self.set_ells(ells=ells)
		
		self.set_catalogues([position1,position2],[weight1,weight2])
		self.run_2pcf_multi_radial(normalize=normalize,los=los,nthreads=nthreads)
		self.free()
	
	def run_2pcf_multi_radial(self,normalize=False,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.radialedges)-1,self.nells)
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_2pcf_multi_radial.argtypes = (typecounts,ctypes.c_size_t,ctypes.c_char_p,ctypes.c_bool,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_2pcf_multi_radial(self.counts,self.nells,self.multitype,normalize,los,nthreads)
		self.counts.shape = shape
		
	def set_4pcf_multi_radial(self,sedges,radialedges,position1,weight1,position2,weight2,position3,weight3,position4,weight4,sbinning='lin',ssize=None,radialbinning='lin',radialsize=None,ells=[0,2,4,6,8,10,12],muedges=[-1.,1.],normalize=False,los='midpoint',nthreads=8):

		self.sedges = self.set_bin('main',sedges,size=ssize,binning=sbinning)
		self.muedges = self.set_bin('aux',muedges,size=1)
		self.radialedges = self.set_bin('radial',radialedges,size=radialsize,binning=radialbinning)
		self.set_ells(ells=ells)
		
		self.set_catalogues([position1,position2,position3,position4],[weight1,weight2,weight3,weight4])
		self.run_4pcf_multi_radial(normalize=normalize,los=los,nthreads=nthreads)
		self.free()
	
	def run_4pcf_multi_radial(self,normalize=False,los='midpoint',nthreads=8):
	
		shape = (len(self.sedges)-1,len(self.sedges)-1,self.nells,self.nells)
	
		self.counts = scipy.zeros(shape,dtype=self.C_TYPE).flatten()
		typecounts = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(len(self.counts)))
	
		self.cute.run_4pcf_multi_radial.argtypes = (typecounts,ctypes.c_size_t,ctypes.c_char_p,ctypes.c_bool,ctypes.c_char_p,ctypes.c_size_t)
		self.cute.run_4pcf_multi_radial(self.counts,self.nells,self.multitype,normalize,los,nthreads)
		self.counts.shape = shape
	
	def set_catalogues(self,positions,weights):
		
		self.cute.clear_catalogs()
		
		def size(tab):
			return tab.shape[-1] if len(tab.shape) > 1 else 1
		
		sizeposition = size(positions[0])
		sizeweight = size(weights[0]) if weights[0] is not None else 1

		for num,(position,weight) in enumerate(zip(positions,weights)):
			if position is not None:
				if size(position) != sizeposition: raise ValueError
				if weight is None:
					weight = scipy.ones((position.shape[0],sizeweight),dtype=position.dtype)
				else:
					assert size(weight) == sizeweight
				self.set_catalogue(num+1,position,weight)	 
		
	def set_catalogue(self,num,position,weight,copy=False):
	
		#First dimension is number of objects
		n = position.shape[0]
		
		def size(tab):
			return tab.shape[-1] if len(tab.shape) > 1 else 1
		
		sizeposition = size(position)
		sizeweight = size(weight)
		
		def new(tab):
			if copy: return scipy.asarray(tab,dtype=self.C_TYPE).flatten()
			return tab.astype(self.C_TYPE).flatten()
		
		self.position[num] = new(position)
		self.weight[num] = new(weight)

		typeposition = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(n*sizeposition))
		typeweight = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=(n*sizeweight))
		
		self.cute.set_catalog.argtypes = (ctypes.c_size_t,typeposition,typeweight,ctypes.c_size_t,ctypes.c_size_t,ctypes.c_size_t)
		self.cute.set_catalog(num,self.position[num],self.weight[num],n,sizeposition,sizeweight)
		
		self.position[num].shape = (n,sizeposition)
		self.weight[num].shape = (n,sizeweight)

	def free(self):
		self.position = {}
		self.weight = {}
		self.cute.clear_catalogs()

