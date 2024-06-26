import os

import numpy as np
from numpy import testing
from kdcount import sphere, correlate

from pycute import *
from pycute.pycute import wrap_phi

nthreads = 8
ns = 10
nspar = 6
nmu = 7
ells = [0,1,2,3,4]
maxsize = 20.
sedges = np.linspace(0.001,10.,ns+1)
radialedges = np.linspace(-4.,4.,nspar+1)
muedges = np.linspace(-1.,1.,nmu+1)
thetaedges = np.linspace(0.01,10.,ns+1)
binsize = 10
verbose = False

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def distance(pos):
    return np.sqrt((pos**2).sum(axis=-1))


def legendre(ell):
    coeffs = [0.]*(ell+1)
    coeffs[ell] = 1.
    return np.polynomial.legendre.Legendre(coeffs)

def save_catalogues(nrand=1000,ncat=4):
    rmax = maxsize
    rmin = -rmax
    for icat in range(1,ncat+1):
        cat = []
        for j in range(3): cat.append(np.random.uniform(rmin,rmax,nrand))
        cat.append(np.random.uniform(0.5,1.,nrand))
        cat.append(np.random.randint(binsize,size=nrand).astype(float))
        cat = np.asarray(cat).T
        np.savetxt('cat_{:d}.txt'.format(icat),cat)


def load_catalogues(ncat=2,return_bin=False,mask=None):
    res = []
    for icat in range(1,ncat+1):
        cat = np.loadtxt('cat_{:d}.txt'.format(icat),unpack=True)
        if return_bin: cols = [cat[:3].T,cat[-2],cat[-1].astype(int)]
        else: cols = [cat[:3].T,cat[-2]]
        if mask is not None: res += [c[mask] for c in cols]
        else: res += cols
    return res


def survey_size(pos):
    diff = pos.max(axis=0)-pos.min(axis=0)
    return np.sqrt((diff**2).sum())


def cartesian_to_angular(pos):
    pos = pos.T
    dist = np.sqrt((pos**2).sum(axis=0))
    costheta = pos[2]/dist
    theta = 90.-np.arccos(costheta)*180./np.pi
    phi = wrap_phi(np.arctan2(pos[1],pos[0]))*180./np.pi
    return np.asarray([phi,theta]).T


def angular_to_cartesian(pos):
    x = np.sin(pos[:,0]*np.pi/180.)
    cos_theta = np.cos(pos[:,0]*np.pi/180.)
    y = cos_theta*np.cos(pos[:,1]*np.pi/180.)
    z = cos_theta*np.cos(pos[:,2]*np.pi/180.)
    return np.asarray([x,y,z]).T


'''
def reference_2pcf_s(sedges,position1,weight1,position2=None,weight2=None):
    from Corrfunc.theory.DD import DD
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


def reference_2pcf_angular(thetaedges,position1,weight1,position2=None,weight2=None):
    from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
    """Reference pair counting via corrfunc"""
    RA2=None if position2 is None else position2[:,0]
    DEC2=None if position2 is None else position2[:,1]
    ref = DDtheta_mocks(position2 is None,nthreads,np.asarray(thetaedges),position1[:,0],position1[:,1],weights1=weight1,RA2=RA2,DEC2=DEC2,weights2=weight2,output_thetaavg=False,weight_type='pair_product',verbose=verbose)
    return ref['npairs']*ref['weightavg']
'''

def reference_2pcf_s(sedges,position1,weight1,position2=None,weight2=None):
    """Reference pair counting via kdcount"""
    tree1 = correlate.points(position1,boxsize=None,weights=weight1)
    factor = 1.
    if position2 is None:
        tree2 = tree1
        factor = 1./2.
    else: tree2 = correlate.points(position2,boxsize=None,weights=weight2)
    bins = correlate.RBinning(np.asarray(sedges))
    pc = correlate.paircount(tree1,tree2,bins,np=0,usefast=False,compute_mean_coords=True)
    return factor*pc.sum1


def reference_2pcf_angular(thetaedges,position1,weight1,position2=None,weight2=None):
    """Reference pair counting via kdcount"""
    tree1 = sphere.points(position1[:,0],position1[:,1],weights=weight1)
    factor = 1.
    if position2 is None:
        tree2 = tree1
        factor = 1./2.
    else: tree2 = sphere.points(position2[:,0],position2[:,1],weights=weight2)
    bins = sphere.AngularBinning(np.asarray(thetaedges))
    pc = correlate.paircount(tree1,tree2,bins)
    return factor*pc.sum1


def reference_2pcf_smu(sedges,muedges,position1,weight1,position2=None,weight2=None,los='midpoint'):
    """Reference pair counting via kdcount"""
    tree1 = correlate.points(position1,boxsize=None,weights=weight1)
    if position2 is None: tree2 = tree1
    else: tree2 = correlate.points(position2,boxsize=None,weights=weight2)
    if los=='midpoint':
        bins = correlate.RmuBinning(np.asarray(sedges),(len(muedges)-1),observer=(0,0,0),mu_min=muedges[0],mu_max=muedges[-1],absmu=False)
    else:
        bins = correlate.FlatSkyBinning(np.asarray(sedges),(len(muedges)-1),los='xyz'.index(los),mu_min=muedges[0],mu_max=muedges[-1],absmu=False)
    pc = correlate.paircount(tree2,tree1,bins,np=0,usefast=False,compute_mean_coords=True)
    return pc.sum1


def reference_2pcf_multi(sedges,position1,weight1,position2=None,weight2=None,ells=[0,1,2,3,4],los='midpoint'):
    """Reference pair counting via kdcount"""
    tree1 = correlate.points(position1,boxsize=None,weights=weight1)
    if position2 is None: tree2 = tree1
    else: tree2 = correlate.points(position2,boxsize=None,weights=weight2)
    if los=='midpoint':
        bins = correlate.MultipoleBinning(np.asarray(sedges),ells)
    else:
        bins = correlate.FlatSkyMultipoleBinning(np.asarray(sedges),ells,los='xyz'.index(los))
    pc = correlate.paircount(tree2,tree1,bins,np=0,usefast=False,compute_mean_coords=True)
    norm = (-1)**np.asarray(ells)*1./(2*np.asarray(ells)+1)
    return pc.sum1.T*norm


def test_2pcf_s():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()

    pycute.set_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    countsref = reference_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

    pycute.set_2pcf_s(sedges,position1,weight1,nthreads=nthreads)
    countsref = reference_2pcf_s(sedges,position1,weight1)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

    pycute.set_2pcf_s(sedges,position1,weight1,position2=position2,nthreads=nthreads)
    countsref = reference_2pcf_s(sedges,position1,weight1,position2=position2)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def test_2pcf_angular():
    position1,weight1,position2,weight2 = load_catalogues()
    position1 = cartesian_to_angular(position1)
    position2 = cartesian_to_angular(position2)

    pycute = PyCute()
    pycute.set_2pcf_angular(thetaedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads,celestial=True)
    countsref = reference_2pcf_angular(thetaedges,position1,weight1,position2=position2,weight2=weight2)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

    pycute.set_2pcf_angular(thetaedges,position1,weight1,nthreads=nthreads,celestial=True)
    countsref = reference_2pcf_angular(thetaedges,position1,weight1)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

    thetaedges_ = [1.,3.,5.,6.]
    pycute.set_2pcf_angular(thetaedges_,position1,weight1,nthreads=nthreads,celestial=True,thetabinning='custom')
    countsref = reference_2pcf_angular(thetaedges_,position1,weight1)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

    thetaedges_ = np.logspace(-2,1,10,base=10)
    pycute.set_2pcf_angular(thetaedges_,position1,weight1,nthreads=nthreads,celestial=True,thetabinning='log')
    countsref = reference_2pcf_angular(thetaedges_,position1,weight1)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def test_2pcf_smu():
    position1,weight1,position2,weight2 = load_catalogues()
    #weight1[:] = 1.; weight2[:] = 1.
    pycute = PyCute()
    """
    pycute.set_2pcf_smu(sedges,muedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    countsref = reference_2pcf_smu(sedges,muedges,position1,weight1,position2=position2,weight2=weight2)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

    pycute.set_2pcf_smu(sedges,muedges,position1,weight1,nthreads=nthreads)
    countsref = reference_2pcf_smu(sedges,muedges,position1,weight1,position2=position1,weight2=weight1)
    counts = pycute.counts + pycute.counts[:,::-1]
    testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)
    """
    sedges_ = [0.,1.,3.,5.,6.]
    pycute.set_2pcf_smu(sedges_,muedges,position1,weight1,nthreads=nthreads,sbinning='custom')
    countsref = reference_2pcf_smu(sedges_,muedges,position1,weight1)
    counts = pycute.counts + pycute.counts[:,::-1]
    testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)

    sedges_ = np.logspace(-2,1,10,base=10)
    pycute.set_2pcf_smu(sedges_,muedges,position1,weight1,nthreads=nthreads,sbinning='log')
    countsref = reference_2pcf_smu(sedges_,muedges,position1,weight1)
    counts = pycute.counts + pycute.counts[:,::-1]
    testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)
    """
    pycute.set_2pcf_smu(sedges,muedges,position1,weight1,nthreads=nthreads,los='x')
    countsref = reference_2pcf_smu(sedges,muedges,position1,weight1,position2=position1,weight2=weight1,los='x')
    counts = pycute.counts + pycute.counts[:,::-1]
    testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)
    """
    #from Corrfunc.theory.DDsmu import DDsmu
    #ref = DDsmu(0,nthreads,pycute.sedges,muedges[-1],nmu,X1=position1[0],Y1=position1[1],Z1=position1[2],weights1=weight1,X2=position2[0],Y2=position2[1],Z2=position2[2],weights2=weight2,weight_type='pair_product',periodic=False,verbose=True,output_savg=True)
    #countsref = (ref['npairs']*ref['weightavg']).reshape((len(sedges)-1,nmu//2))


def test_2pcf_multi():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    for ells in [[0,1,2,3,4,5,6],[0,2,4,6,8],[1,3,5,7,9]]:

        pycute.set_2pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,muedges=[-1.,1.])
        countsref = reference_2pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells)
        testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

        pycute.set_2pcf_multi(sedges,position1,weight1,ells=ells,nthreads=nthreads,muedges=[-1.,1.])
        countsref = reference_2pcf_multi(sedges,position1,weight1,ells=ells)
        counts = (1 + (-1)**np.asarray(pycute.ells))*pycute.counts
        testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)

        pycute.set_2pcf_multi(sedges,position1,weight1,ells=ells,nthreads=nthreads,muedges=[-1.,1.],los=[1.,0.,0.])
        countsref = reference_2pcf_multi(sedges,position1,weight1,ells=ells,los='x')
        counts = (1 + (-1)**np.asarray(pycute.ells))*pycute.counts
        testing.assert_allclose(countsref,counts,rtol=1e-7,atol=1e-7)


def test_2pcf_scos():
    sedges = [0,100]
    position1,weight1,position2,weight2 = load_catalogues()
    #weight1[:] = 1.; weight2[:] = 1.
    pycute = PyCute()
    pycute.set_2pcf_scos(sedges,muedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    countsref = reference_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2)
    testing.assert_allclose(countsref,np.sum(pycute.counts,axis=1),rtol=1e-7,atol=1e-7)


def save_reference_3pcf_multi():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_3pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads)
    np.save('ref_3pcf_multi.npy',pycute.counts)


def reference_3pcf_multi():
    return np.load('ref_3pcf_multi.npy')


def test_3pcf_multi():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
    countsref = reference_3pcf_multi()
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
    pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,position3=position2,weight3=weight2,ells=[ells,ells[::2]],nthreads=nthreads)
    testing.assert_allclose(countsref[...,::2],pycute.counts,rtol=1e-7,atol=1e-7)
    pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,position3=position2,weight3=weight2,ells=[ells,ells[::2]],los=['endpoint','firstpoint'],losn=[1,2],nthreads=nthreads)
    assert np.isnan(pycute.counts).sum() == 0
    assert np.any(pycute.counts)


def test_3pcf_multi_double_los():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
    countsref = 2.*pycute.counts
    countsref[...,np.mod(ells,2)==1,:] = 0.
    pycute.set_3pcf_multi_double_los(sedges,position1,weight1,position2,weight2,ells=ells,los='midpoint',nthreads=nthreads)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
    pycute.set_3pcf_multi_double_los(sedges,position1,weight1,position2,weight2,position3=position2,weight3=weight2,ells=ells,los=['midpoint','midpoint'],nthreads=nthreads)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def save_reference_2pcf_multi_binned():
    position1,weight1,bin1,position2,weight2,bin2 = load_catalogues(return_bin=True)
    pycute = PyCute()
    pycute.set_2pcf_multi_binned(sedges,binsize,position1,weight1,position2=position2,weight2=weight2,bin2=bin2,ells=ells,nthreads=nthreads)
    np.save('ref_2pcf_multi_binned.npy',pycute.counts)


def reference_2pcf_multi_binned():
    return np.load('ref_2pcf_multi_binned.npy')


def test_2pcf_multi_binned():
    position1,weight1,bin1,position2,weight2,bin2 = load_catalogues(return_bin=True)
    pycute = PyCute()
    pycute.set_2pcf_multi_binned(sedges,binsize,position1,weight1,position2=position2,weight2=weight2,bin2=bin2,ells=ells,nthreads=nthreads)
    countsref = reference_2pcf_multi_binned()
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def test_4pcf_multi_binned():
    position1,weight1,bin1,position2,weight2,bin2,position3,weight3,bin3,position4,weight4,bin4 = load_catalogues(4,return_bin=True)
    pycute = PyCute()
    pycute.set_2pcf_multi_binned(sedges,binsize,position1,weight1,position2=position2,weight2=weight2,bin2=bin2,ells=ells,nthreads=nthreads)
    counts1 = pycute.counts
    pycute.set_2pcf_multi_binned(sedges,binsize,position3,weight3,position2=position4,weight2=weight4,bin2=bin4,ells=ells,nthreads=nthreads)
    counts2 = pycute.counts
    countsref = np.einsum('irj,krl->ikjl',counts1,counts2)
    pycute.set_4pcf_multi_binned(sedges,binsize,position1,weight1,position2=position2,weight2=weight2,bin2=bin2,position3=position3,weight3=weight3,position4=position4,weight4=weight4,bin4=bin4,ells=ells,nthreads=nthreads)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def test_comp_3pcf_multi():
    position1,weight1,position2,weight2 = load_catalogues()
    position3,weight3 = position2.copy(),weight2.copy()
    weight3 /= np.sum(weight3)
    pycute = PyCute()
    size = survey_size(position2)*1.2
    sedges = np.linspace(0.,size,10)
    for ells in [[0,1,2]]:
        pycute.set_2pcf_multi(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
        countsref = pycute.counts
        pycute.set_3pcf_multi(sedges,position1,weight1,position2,weight2,position3=position3,weight3=weight3,ells=ells,nthreads=nthreads)
        countstest = np.sum(pycute.counts[...,0,:],axis=0)
        testing.assert_allclose(countsref,countstest,rtol=1e-7,atol=1e-7)


def test_los():

    position1,weight1,position2,weight2 = load_catalogues()
    offset = 100.
    position1 = position1 + offset
    position2 = position2 + offset
    pycute = PyCute()
    sedges = np.linspace(0.001,100.,ns+1)
    for ells in [[0,1,2,3,4,5,6],[0,2,4,6,8],[1,3,5,7,9]]:

        pycute.set_2pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,los='midpoint')
        middle = np.array(pycute.counts)
        pycute.set_2pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,los='endpoint')
        #print np.absolute(pycute.counts/middle-1.)[middle>0.].max()

    ells = [1,3,5,7,9]
    pycute.set_2pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,los='endpoint')
    countsref = pycute.counts
    pycute.set_2pcf_multi(sedges,position1,weight1*distance(position1),position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,los='endpoint',losn=1)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
    pycute.set_2pcf_multi(sedges,position1,weight1*distance(position1)**2,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,los='endpoint',losn=2)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def test_weight():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    countsref = pycute.counts
    weight1_ = np.array([weight1,3.*np.ones_like(weight1)]).T
    weight2_ = np.array([weight2,np.ones_like(weight2)]).T
    pycute.set_2pcf_s(sedges,position1,weight1_,position2=position2,weight2=weight2_,nthreads=nthreads,weighttype='prodsum')
    testing.assert_allclose(countsref,pycute.counts/4.,rtol=1e-7,atol=1e-7)
    weight1_ = np.array([weight1,np.zeros_like(weight1)]).T
    weight2_ = np.array([np.ones_like(weight2),weight2]).T
    pycute.set_2pcf_s(sedges,position1,weight1_,position2=position2,weight2=weight2_,nthreads=nthreads,weighttype='prodsum')
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)

    pycute.set_2pcf_s(sedges,position1,weight1,position2=position1,weight2=weight1,nthreads=nthreads)
    countsref = pycute.counts
    pycute.set_2pcf_s(sedges,position1,weight1,nthreads=nthreads,weighttype='prod')
    testing.assert_allclose(countsref,pycute.counts*2.,rtol=1e-7,atol=1e-7)
    countsref = pycute.counts
    weight1_ = np.array([weight1,np.ones_like(weight1)]).T
    pycute.set_2pcf_s(sedges,position1,weight1_,nthreads=nthreads,weighttype='prodsum')
    testing.assert_allclose(countsref,pycute.counts/2.,rtol=1e-7,atol=1e-7)


def save_reference_2pcf_multi_radial_legendre():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_2pcf_multi_radial_legendre(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
    np.save('ref_2pcf_multi_radial_legendre.npy',pycute.counts)


def reference_2pcf_multi_radial_legendre():
    return np.load('ref_2pcf_multi_radial_legendre.npy')


def test_2pcf_multi_radial_legendre():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_2pcf_multi_radial_legendre(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
    countsref = reference_2pcf_multi_radial_legendre()
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def save_reference_2pcf_multi_angular_legendre():
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_2pcf_multi_angular_legendre(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=nthreads)
    np.save('ref_2pcf_multi_angular_legendre.npy',pycute.counts)


def reference_2pcf_multi_angular_legendre():
    return np.load('ref_2pcf_multi_angular_legendre.npy')


def test_2pcf_multi_angular_legendre():
    #ells = [0,2,4,6,8,10]
    position1,weight1,position2,weight2 = load_catalogues()
    #position1 += 1000.
    #position2 += 1000.
    pycute = PyCute()
    pycute.set_2pcf_multi_angular_legendre(sedges,position1,weight1,position2,weight2,ells=ells,nthreads=1)
    countsref = reference_2pcf_multi_angular_legendre()
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)


def ref_integrate_legendre(self,ells=[0]):

    midmu = (self.muedges[:-1]+self.muedges[1:])/2.
    res = [(self.counts*legendre(ell)(midmu)).sum(axis=-1) for ell in ells]
    return np.asarray(res).T


def ref_integrate_legendre_bis(self,ells=[0]):

    res = []
    for ell in ells:
        leg = legendre(ell)(self.muedges)
        res.append((self.counts*(leg[:-1]+leg[1:])/2.).sum(axis=-1))
    return np.asarray(res).T


def test_integrate_legendre():

    muedges_ = np.linspace(-1.,1.,1000)
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_2pcf_smu(sedges,muedges_,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    pycute.integrate_legendre(ells,nthreads=nthreads)
    integralref = ref_integrate_legendre(pycute,ells=ells)
    integralbis = ref_integrate_legendre_bis(pycute,ells=ells)
    pycute.set_2pcf_multi(sedges,position1,weight1,position2=position2,weight2=weight2,ells=ells,nthreads=nthreads,muedges=[-1.,1.])
    testing.assert_allclose(integralref,pycute.integral,rtol=1e-7,atol=1e-7)
    print('Double v.s. simple (abs,rel) = {:.4g}, {:.4g}'.format(np.absolute(integralbis-integralref).max(),np.absolute((integralbis-integralref)/integralref).max()))
    print('Double v.s. truth (abs,rel) = {:.4g}, {:.4g}'.format(np.absolute(integralbis-pycute.counts).max(),np.absolute((integralbis-pycute.counts)/pycute.counts).max()))
    print('Simple v.s. truth (abs,rel) = {:.4g}, {:.4g}'.format(np.absolute(integralref-pycute.counts).max(),np.absolute((integralref-pycute.counts)/pycute.counts).max()))


def ref_integrate_radial_legendre(self,ells=[[0],[0]]):

    midmu = (self.muedges[:-1]+self.muedges[1:])/2.
    mids = (self.sedges[:-1]+self.sedges[1:])/2.
    meshd,meshs = np.meshgrid(mids,mids,sparse=False,indexing='ij')

    res = []
    for ell1 in ells[0]:
        res.append([])
        weight = (self.counts*legendre(ell1)(midmu)).T
        for ell2 in ells[1]:
            midmu2 = midmu[:,None,None]*meshd/meshs #mu,d,s
            tmp = weight[:,:,None]*legendre(ell2)(midmu2)
            tmp /= mids[None,None,:]
            tmp[np.absolute(midmu2)>1.] = 0.
            res[-1].append(tmp.sum(axis=0))

    return np.transpose(np.asarray(res),axes=(2,3,0,1))


def test_integrate_radial_legendre():

    muedges_ = np.linspace(-1.,1.,1000)
    sedges = np.linspace(0.001,10.,100)
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_2pcf_smu(sedges,muedges_,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    pycute.integrate_radial_legendre([ells,ells],nthreads=nthreads)
    integralref = ref_integrate_radial_legendre(pycute,ells=[ells,ells])
    testing.assert_allclose(integralref,pycute.integral,rtol=1e-7,atol=1e-7)


def ref_integrate_angular_legendre(self,ells=[[0],[0]]):

    midmu = (self.muedges[:-1]+self.muedges[1:])/2.
    mids = (self.sedges[:-1]+self.sedges[1:])/2.
    meshd,meshs = np.meshgrid(mids,mids,sparse=False,indexing='ij')

    res = []
    for ell1 in ells[0]:
        res.append([])
        weight = (self.counts*legendre(ell1)(midmu)).T
        for ell2 in ells[1]:
            midmu2 = 1.-(meshd/meshs)**2*(1.-midmu**2)[:,None,None]
            midmu2[midmu2<0.] = np.nan
            midmu2 = np.sqrt(midmu2)
            tmp = weight[:,:,None]*(legendre(ell2)(midmu2)+legendre(ell2)(-midmu2))
            tmp /= (mids[None,None,:])**2*midmu2
            tmp[np.isnan(midmu2)] = 0.
            res[-1].append(tmp.sum(axis=0))

    return np.transpose(np.asarray(res),axes=(2,3,0,1))


def test_integrate_angular_legendre():

    muedges_ = np.linspace(-1.,1.,201)
    sedges = np.linspace(0.001,10.,101)
    position1,weight1,position2,weight2 = load_catalogues()
    pycute = PyCute()
    pycute.set_2pcf_smu(sedges,muedges_,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    pycute.integrate_angular_legendre([ells,ells],nthreads=nthreads)
    integralref = ref_integrate_angular_legendre(pycute,ells=[ells,ells])
    testing.assert_allclose(integralref,pycute.integral,rtol=1e-7,atol=1e-7)


def test_verbosity():
    position1,weight1,position2,weight2 = load_catalogues()
    print('No output in between <<')
    pycute = PyCute()
    pycute.set_verbosity('quiet')
    pycute.set_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    print('>>')
    countsref = reference_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2)
    testing.assert_allclose(countsref,pycute.counts,rtol=1e-7,atol=1e-7)
    """
    pycute = PyCute()
    pycute.set_verbosity('debug')
    pycute.set_2pcf_s(sedges,position1,weight1,position2=position2,weight2=weight2,nthreads=nthreads)
    """

if __name__ == '__main__':

    """
    save_catalogues()
    save_reference_3pcf_multi()
    save_reference_2pcf_multi_binned()
    save_reference_2pcf_multi_radial_legendre()
    save_reference_2pcf_multi_angular_legendre()
    """
    test_2pcf_s()
    test_2pcf_angular()
    test_2pcf_smu()
    test_2pcf_scos()
    test_los()
    test_weight()
    test_2pcf_multi()
    test_3pcf_multi()
    test_3pcf_multi_double_los()
    test_2pcf_multi_binned()
    test_4pcf_multi_binned()
    test_comp_3pcf_multi()
    test_2pcf_multi_radial_legendre()
    test_2pcf_multi_angular_legendre()
    test_integrate_legendre()
    test_integrate_radial_legendre()
    test_integrate_angular_legendre()

    test_verbosity()
