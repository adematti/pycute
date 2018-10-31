# coding: utf8

import os
import scipy
from scipy import constants
from numpy import testing
from Corrfunc.theory.DDsmu import DDsmu
from Corrfunc.theory.DD import DD
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from kdcount import correlate
from pycute import *
from pycute.pyCute import wrap_phi

nthreads = 8
ns = 10
nspar = 6
nmu = 4
ells = [0,1,2,3,4]
maxsize = 20.
sedges = scipy.linspace(0.001,10.,ns+1)
radialedges = scipy.linspace(-4.,4.,nspar+1)
muedges = scipy.linspace(-1.,1.,nmu+1)
thetaedges = scipy.linspace(0.01,10.,ns+1)
verbose = False

def save_catalogues(nrand=1000,ncat=4):
	rmax = maxsize
	rmin = -rmax
	for icat in range(1,ncat+1):
		cat = []
		for j in range(3): cat.append(scipy.random.uniform(rmin,rmax,nrand))
		cat.append(scipy.random.uniform(0.5,1.,nrand))
		cat = scipy.asarray(cat).T
		scipy.savetxt('cat_{:d}.txt'.format(icat),cat)

def load_catalogues(ncat=2):
	res = []
	for icat in range(1,ncat+1):
		cat = scipy.loadtxt('cat_{:d}.txt'.format(icat),unpack=True)
		res += [cat[:3].T,cat[-1]]
	return res
	
def survey_size(pos):
	diff = pos.max(axis=0)-pos.min(axis=0)
	return scipy.sqrt((diff**2).sum())

def cartesian_to_angular(pos):
	pos = pos.T
	dist = scipy.sqrt((pos**2).sum(axis=0))
	costheta = pos[2]/dist
	theta = 90.-scipy.arccos(costheta)/constants.degree
	phi = wrap_phi(scipy.arctan2(pos[1],pos[0]))/constants.degree
	return scipy.asarray([phi,theta]).T

def reference_s(sedges,position1,weight1,position2=None,weight2=None):
	"""Reference pair counting via corrfunc"""
	if position2 is None:
		X2 = Y2 = Z2 = None
		factor = 0.5
	else:
		X2 = position2[:,0]
		Y2 = position2[:,1]
		Z2 = position2[:,2]
		factor = 1.
	ref = DD(position2 is None,nthreads,sedges,X1=position1[:,0],Y1=position1[:,1],Z1=position1[:,2],weights1=weight1,X2=X2,Y2=Y2,Z2=Z2,weights2=weight2,weight_type='pair_product',periodic=False,verbose=verbose)
	return factor*ref['npairs']*ref['weightavg']
	
def reference_angular(thetaedges,position1,weight1,position2=None,weight2=None):
	"""Reference pair counting via corrfunc"""
	RA2=None if position2 is None else position2[:,0]
	DEC2=None if position2 is None else position2[:,1]
	ref = DDtheta_mocks(position2 is None,nthreads,scipy.asarray(thetaedges),position1[:,0],position1[:,1],weights1=weight1,RA2=RA2,DEC2=DEC2,weights2=weight2,output_thetaavg=False,weight_type='pair_product',verbose=verbose)
	return ref['npairs']*ref['weightavg']
	
def reference_smu(sedges,muedges,position1,weight1,position2=None,weight2=None):
	"""Reference pair counting via kdcount"""
	tree1 = correlate.points(position1,boxsize=None,weights=weight1)
	if position2 is None: tree2 = tree1
	else: tree2 = correlate.points(position2,boxsize=None,weights=weight2)
	bins = correlate.RmuBinning(scipy.asarray(sedges),(len(muedges)-1),observer=(0,0,0),mu_min=muedges[0],mu_max=muedges[-1],absmu=False)
	pc = correlate.paircount(tree2,tree1,bins,np=0,usefast=False,compute_mean_coords=True)
	return pc.sum1
	
def reference_multi(sedges,position1,weight1,position2=None,weight2=None,ells=[0,1,2,3,4]):
	"""Reference pair counting via kdcount"""
	tree1 = correlate.points(position1,boxsize=None,weights=weight1)
	if position2 is None: tree2 = tree1
	else: tree2 = correlate.points(position2,boxsize=None,weights=weight2)
	bins = correlate.MultipoleBinning(scipy.asarray(sedges),ells)
	pc = correlate.paircount(tree2,tree1,bins,np=0,usefast=False,compute_mean_coords=True)
	norm = 2*scipy.asarray(ells)+1
	return pc.sum1.T/norm
	
def test_s():
	position1,weight1,position2,weight2 = load_catalogues()
	pycute = PyCute()
	
	pycute.set_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
	countsref = reference_s(sedges,position1,weight1,position2=position2,weight2=weight2)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
	
	pycute.set_2pcf_s(sedges,position1,weight1,nthreads=nthreads)
	countsref = reference_s(sedges,position1,weight1)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
	
	pycute.set_2pcf_s(sedges,position1,weight1,position2=position2,nthreads=nthreads)
	countsref = reference_s(sedges,position1,weight1,position2=position2)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

def test_angular():
	position1,weight1,position2,weight2 = load_catalogues()
	position1 = cartesian_to_angular(position1)
	position2 = cartesian_to_angular(position2)
	
	pycute = PyCute()
	pycute.set_2pcf_angular(thetaedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads,celestial=True)
	countsref = reference_angular(thetaedges,position1,weight1,position2=position2,weight2=weight2)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
	
	pycute.set_2pcf_angular(thetaedges,position1,weight1,nthreads=nthreads,celestial=True)
	countsref = reference_angular(thetaedges,position1,weight1)
	testing.assert_allclose(countsref,2.*pycute.counts,rtol=1e-7,atol=1e-7)
	
	thetaedges_ = [1.,3.,5.,6.]
	pycute.set_2pcf_angular(thetaedges_,position1,weight1,nthreads=nthreads,celestial=True,thetabinning='custom')
	countsref = reference_angular(thetaedges_,position1,weight1)
	testing.assert_allclose(countsref,2.*pycute.counts,rtol=1e-7,atol=1e-7)
	
	thetaedges_ = scipy.logspace(-2,1,10,base=10)
	pycute.set_2pcf_angular(thetaedges_,position1,weight1,nthreads=nthreads,celestial=True,thetabinning='log')
	countsref = reference_angular(thetaedges_,position1,weight1)
	testing.assert_allclose(countsref,2.*pycute.counts,rtol=1e-7,atol=1e-7)


def test_smu():
	position1,weight1,position2,weight2 = load_catalogues()
	#weight1[:] = 1.; weight2[:] = 1.
	pycute = PyCute()
	
	pycute.set_2pcf_smu(sedges,muedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
	countsref = reference_smu(sedges,muedges,position1,weight1,position2=position2,weight2=weight2)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
	
	pycute.set_2pcf_smu(sedges,muedges,position1,weight1,nthreads=nthreads)
	countsref = reference_smu(sedges,muedges,position1,weight1,position2=position1,weight2=weight1)
	counts = pycute.counts + pycute.counts[:,::-1]
	testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)
	
	sedges_ = [0.,1.,3.,5.,6.]
	pycute.set_2pcf_smu(sedges_,muedges,position1,weight1,nthreads=nthreads,sbinning='custom')
	countsref = reference_smu(sedges_,muedges,position1,weight1,position2=position1,weight2=weight1)
	counts = pycute.counts + pycute.counts[:,::-1]
	testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)
	
	sedges_ = scipy.logspace(-2,1,10,base=10)
	pycute.set_2pcf_smu(sedges_,muedges,position1,weight1,nthreads=nthreads,sbinning='log')
	countsref = reference_smu(sedges_,muedges,position1,weight1,position2=position1,weight2=weight1)
	counts = pycute.counts + pycute.counts[:,::-1]
	testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)
	
	#ref = DDsmu(0,nthreads,pycute.sedges,muedges[-1],nmu,X1=position1[0],Y1=position1[1],Z1=position1[2],weights1=weight1,X2=position2[0],Y2=position2[1],Z2=position2[2],weights2=weight2,weight_type='pair_product',periodic=False,verbose=True,output_savg=True)
	#countsref = (ref['npairs']*ref['weightavg']).reshape((len(sedges)-1,nmu//2))
	
def test_multi():
	position1,weight1,position2,weight2 = load_catalogues()
	pycute = PyCute()
	for ells in [[0,1,2,3,4,5,6],[0,2,4,6,8],[1,3,5,7,9]]:
		pycute.set_2pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,muedges=[-1.,1.])
		countsref = reference_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells)
		testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
		pycute.set_2pcf_multi(sedges,position1,weight1,ells=ells,nthreads=nthreads,muedges=[-1.,1.])
		countsref = reference_multi(sedges,position1,weight1,ells=ells)
		counts = (1 + (-1)**scipy.asarray(ells))*pycute.counts
		testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)
	
def test_scos():
	sedges = [0,100]
	position1,weight1,position2,weight2 = load_catalogues()
	#weight1[:] = 1.; weight2[:] = 1.
	pycute = PyCute()
	pycute.set_2pcf_scos(sedges,muedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
	countsref = reference_s(sedges,position1,weight1,position2=position2,weight2=weight2)
	testing.assert_allclose(countsref,scipy.sum(pycute.counts,axis=1),rtol=1e-7,atol=1e-7)
	
def save_reference_3pcf_multi():
	position1,weight1,position2,weight2 = load_catalogues()
	pycute = PyCute()
	pycute.set_3pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads)
	scipy.save('ref_3pcf_multi.npy',pycute.counts)
		
def reference_3pcf_multi():
	return scipy.load('ref_3pcf_multi.npy')

def test_3pcf_multi():
	position1,weight1,position2,weight2 = load_catalogues()
	pycute = PyCute()
	pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
	countsref = reference_3pcf_multi()
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
	pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,position3=position2,weight3=weight2,ells=[ells,ells[::2]],nthreads=nthreads)
	testing.assert_allclose(countsref[...,::2],pycute.counts,rtol=1e-7,atol=1e-7)
	
def test_3pcf_multi_double_los():
	position1,weight1,position2,weight2 = load_catalogues()
	pycute = PyCute()
	pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
	countsref = 2.*pycute.counts
	countsref[...,scipy.mod(ells,2)==1,:] = 0.
	pycute.set_3pcf_multi_double_los(sedges,position1,weight1,position2,weight2,ells=ells,los='midpoint',nthreads=nthreads)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
	pycute.set_3pcf_multi_double_los(sedges,position1,weight1,position2,weight2,position3=position2,weight3=weight2,ells=ells,los='midpoint',nthreads=nthreads)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
	
def save_reference_2pcf_multi_radial():
	position1,weight1,position2,weight2 = load_catalogues()
	pycute = PyCute()
	pycute.set_2pcf_multi_radial(sedges,radialedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
	scipy.save('ref_2pcf_multi_radial.npy',pycute.counts)
		
def reference_2pcf_multi_radial():
	return scipy.load('ref_2pcf_multi_radial.npy')

def test_2pcf_multi_radial():
	position1,weight1,position2,weight2 = load_catalogues()
	pycute = PyCute()
	pycute.set_2pcf_multi_radial(sedges,radialedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
	countsref = reference_2pcf_multi_radial()
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

def test_4pcf_multi_radial():
	position1,weight1,position2,weight2,position3,weight3,position4,weight4 = load_catalogues(4)
	pycute = PyCute()
	pycute.set_2pcf_multi_radial(sedges,radialedges,position1,weight1,position2,weight2,ells=ells,normalize=False,nthreads=nthreads)
	counts1 = pycute.counts
	pycute.set_2pcf_multi_radial(sedges,radialedges,position3,weight3,position4,weight4,ells=ells,normalize=True,nthreads=nthreads)
	counts2 = pycute.counts
	countsref = scipy.einsum('irj,krl->ikjl',counts1,counts2)
	pycute.set_4pcf_multi_radial(sedges,radialedges,position1,weight1,position2,weight2,position3,weight3,position4,weight4,ells=ells,normalize=True,nthreads=nthreads)
	testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

def comp_3pcf_multi_radial():
	position1,weight1,position2,weight2 = load_catalogues()
	position3,weight3 = position2.copy(),weight2.copy()
	weight3 /= scipy.sum(weight3)
	pycute = PyCute()
	size = survey_size(position2)*1.2
	sedges = scipy.linspace(0.,size,10)
	radialedges = scipy.linspace(-size,size,10+1)
	for ells in [[0,1,2]]:
		pycute.set_2pcf_multi(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
		countsref = pycute.counts
		pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,position3=position3,weight3=weight3,ells=ells,nthreads=nthreads)
		countstest = scipy.sum(pycute.counts[...,0,:],axis=0)
		testing.assert_allclose(countsref,countstest,rtol=1e-7,atol=1e-7)


"""
save_catalogues()
save_reference_3pcf_multi()
save_reference_2pcf_multi_radial()
save_reference_3pcf_multi_radial()
"""

test_s()
test_angular()
test_smu()
test_multi()
test_scos()
test_3pcf_multi()
test_3pcf_multi_double_los()
test_2pcf_multi_radial()
test_4pcf_multi_radial()
comp_3pcf_multi_radial()


