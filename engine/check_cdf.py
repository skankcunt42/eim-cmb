import numpy as np, healpy as hp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rng=np.random.default_rng(99)
NSIDE=64; LMAX=32; npix=hp.nside2npix(NSIDE)
VEC=np.array(hp.pix2vec(NSIDE,np.arange(npix)))
try:
    import camb
    pars=camb.set_params(H0=67.36,ombh2=0.02237,omch2=0.1200,ns=0.9649,As=2.1e-9,tau=0.0544,mnu=0.06)
    pars.set_for_lmax(LMAX+40); res=camb.get_results(pars)
    cl=np.asarray(res.get_cmb_power_spectra(pars,CMB_unit='muK',raw_cl=True)['total'][:LMAX+1,0])
except Exception:
    l=np.arange(LMAX+1); cl=np.zeros(LMAX+1); cl[2:]=2*np.pi*1000/(l[2:]*(l[2:]+1))
cl[:2]=0
def lset_map(alm,lset):
    a=np.zeros_like(alm)
    for L in lset:
        for m in range(L+1): a[hp.Alm.getidx(LMAX,L,m)]=alm[hp.Alm.getidx(LMAX,L,m)]
    return hp.alm2map(a,NSIDE,verbose=False)
def axis_of(alm,lset):
    m=lset_map(alm,lset); w=m*m
    M=np.einsum('p,ip,jp->ij',w,VEC,VEC); ev,evec=np.linalg.eigh(M)
    n=evec[:,0]; return n/np.linalg.norm(n)
def ang(a,b): return float(np.degrees(np.arccos(np.clip(abs(np.dot(a,b)),0,1))))

for N,tag in [(250,'low'),(8000,'high')]:
    a23=np.zeros(N)
    for i in range(N):
        alm=hp.synalm(cl,lmax=LMAX); a23[i]=ang(axis_of(alm,[2]),axis_of(alm,[3]))
    if tag=='low': a_low=a23
    else: a_high=a23
# exact percentile of observed ~10 deg
pct10=(a_high<=10).mean()
print("N=%d null:  median=%.2f deg (iso 60)"%(len(a_high),np.median(a_high)))
print("MC percentile of alpha<=10deg = %.4f   analytic 1-cos10 = %.4f"%(pct10,1-np.cos(np.radians(10))))
print("counts below thresholds (analytic expectation in []):")
for th in [2,5,10,20]:
    exp=(1-np.cos(np.radians(th)))*len(a_high)
    print("   alpha<%2ddeg: %4d  [analytic %.0f]"%(th,(a_high<th).sum(),exp))
# max deviation empirical CDF vs analytic (KS-style) => tests for any floor/selection
xs=np.sort(a_high); F=np.arange(1,len(xs)+1)/len(xs); Fan=1-np.cos(np.radians(xs))
ksd=np.max(np.abs(F-Fan)); print("max|CDF_emp - (1-cos a)| = %.4f (small => no floor/selection; matches isotropic)"%ksd)

fig,ax=plt.subplots(1,2,figsize=(12,4.4))
# left: empirical CDF vs analytic
ax[0].plot(xs,F,color='#4C78A8',lw=2,label='empirical CDF (N=8000)')
aa=np.linspace(0,90,300); ax[0].plot(aa,1-np.cos(np.radians(aa)),'k--',lw=1.5,label='analytic 1−cos α')
ax[0].axvline(10,color='#E45756',lw=2,label='observed ~10°  (pct=%.3f)'%pct10)
ax[0].set_xlabel('α between even & odd axis (deg)'); ax[0].set_ylabel('P(≤ α)'); ax[0].set_title('Empirical CDF vs isotropic null'); ax[0].legend(fontsize=8)
# right: histogram gaps low-N vs high-N
bins=np.linspace(0,90,31)
ax[1].hist(a_low,bins=bins,density=True,color='#B0C4DE',alpha=.9,label='histogram N=250 (gappy)')
ax[1].hist(a_high,bins=bins,density=True,histtype='step',color='#1f3a5f',lw=2,label='histogram N=8000 (smooth)')
ax[1].plot(aa,np.sin(np.radians(aa))*(np.pi/180),'k--',lw=1.5,label='analytic sin α')
ax[1].set_xlabel('α (deg)'); ax[1].set_ylabel('density'); ax[1].set_title('Low-angle “gaps” are binning, not structure'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('coaxiality_cdf.png',dpi=95); print("saved coaxiality_cdf.png")
