"""Pairs-bootstrap reference calculations for Umediation simulation studies.

Within-group sufficient statistics reproduce whole-cluster resampling.
Common random numbers permit exact reuse across structural path parameters.
"""
import numpy as np

STATNAMES=['a1','a2','bpath','c1','c2','eta','kappa_m','kappa_d','kappa_t','ie_l','ie_u','delta_ie','ie1','ie2','ie3']
def rng_for(st,rep,stream):
    return np.random.default_rng(np.random.SeedSequence([20261021,st['id'],rep,stream]))

def generate(st,rep):
    rng=rng_for(st,rep,0); g=st['g']
    if g:
        sizes=rng.integers(4,13,size=g) if st['t']==0 else np.full(g,st['t'])
        firm=np.repeat(np.arange(g),sizes)
        year=np.concatenate([np.arange(k) for k in sizes])
    else:
        firm=np.arange(st['n']); year=np.zeros(st['n'],int); sizes=np.ones(st['n'],int)
    n=len(firm); x=rng.uniform(-2,2,n); x[0]=-1.9; x[-1]=1.9
    z=rng.normal(size=n)
    if not g:
        u=rng.normal(size=n); v=rng.normal(size=n)
        if st['error']=='hetero':
            u*=.45+.55*np.abs(x); v*=.45+.55*np.abs(x)
        fe_m=np.zeros(n); fe_y=np.zeros(n)
    else:
        t3=st['error']=='t3_ar1'
        innovations=(rng.standard_t(3,(2,n))/np.sqrt(3)) if t3 else rng.normal(size=(2,n))
        scale=np.exp(.35*rng.normal(size=g)-.35**2/2)
        uv=np.empty((2,n)); offset=0; rho=.65
        for j,k in enumerate(sizes):
            uv[:,offset]=innovations[:,offset]*scale[j]
            for t in range(1,k):
                uv[:,offset+t]=rho*uv[:,offset+t-1]+np.sqrt(1-rho*rho)*innovations[:,offset+t]*scale[j]
            offset+=k
        u,v=uv
        means=np.bincount(firm,weights=x)/sizes
        fe_m=(.4*rng.normal(size=g)+.3*means)[firm]
        fe_y=(.4*rng.normal(size=g)-.2*means)[firm]
        if st['year_fe']:
            fe_m=fe_m+.12*np.sin(year); fe_y=fe_y+.12*np.cos(year)
    return dict(x=x,z=z,u=u,v=v,firm=firm,year=year,sizes=sizes,fe_m=fe_m,fe_y=fe_y,n=n)

def responses(data,case):
    x,z=data['x'],data['z']
    m=case['a1']*x+case['a2']*x*x+.3*z+data['fe_m']+data['u']
    y=case['c1']*x+case['c2']*x*x+case['b']*m+.2*z+data['fe_y']+data['v']
    return m,y

def sufficient(data,st):
    x=data['x']; n=len(x)
    columns=[x,x*x,data['z']]
    if st['g'] and st['year_fe']:
        columns.extend([(data['year']==j).astype(float) for j in range(1,st['t'])])
    if not st['g']: columns.append(np.ones(n))
    q=np.column_stack(columns+[data['u'],data['v']]); p=q.shape[1]-2
    if st['g']:
        starts=np.r_[0,np.cumsum(data['sizes'])[:-1]]
        means=np.add.reduceat(q,starts,axis=0)/data['sizes'][:,None]
        q=q-means[data['firm']]
    tri=np.triu_indices(q.shape[1]); products=q[:,tri[0]]*q[:,tri[1]]
    if st['g']: products=np.add.reduceat(products,starts,axis=0)
    return products,p,tri

def projections(products,p,tri,weights):
    flat=weights@products
    if flat.ndim==1: flat=flat[None,:]
    d=p+2; s=np.zeros((len(flat),d,d))
    s[:,tri[0],tri[1]]=flat; s[:,tri[1],tri[0]]=flat
    xx=s[:,:p,:p]; xy=s[:,:p,p:]
    # Balanced designs can become singular when only one firm is selected.
    eig=np.linalg.eigvalsh(xx)
    identified=(eig[:,0]>np.maximum(1.,eig[:,-1])*1e-11)
    co=np.full((len(s),p,2),np.nan)
    co[identified]=np.linalg.solve(xx[identified],xy[identified])
    ru=s[:,p,p]-np.sum(xy[:,:,0]*co[:,:,0],axis=1)
    ruv=s[:,p,p+1]-np.sum(xy[:,:,0]*co[:,:,1],axis=1)
    valid=identified & np.isfinite(ru) & (ru>1e-11*np.maximum(1.,s[:,p,p]))
    db=np.full(len(s),np.nan); db[valid]=ruv[valid]/ru[valid]
    out=np.column_stack([co[:,0,0],co[:,1,0],db,co[:,0,1],co[:,1,1]])
    out[~valid]=np.nan
    return out,valid

def evaluate(proj,case):
    a1=case['a1']+proj[:,0]; a2=case['a2']+proj[:,1]; b=case['b']+proj[:,2]
    c1=case['c1']+proj[:,3]-proj[:,2]*a1
    c2=case['c2']+proj[:,4]-proj[:,2]*a2
    eta=a1*b; km=2*a2*b; kd=2*c2; kt=km+kd
    left=eta-1.5*km; right=eta+1.5*km
    return np.column_stack([a1,a2,b,c1,c2,eta,km,kd,kt,left,right,3*km,left,eta,right])

def summarize(point,boot,case,valid_count,requested,n):
    q=np.quantile(boot,[.025,.05,.95,.975],axis=0,method='linear')
    se=boot.std(axis=0,ddof=1)
    row={'n':n,'requested':requested,'valid':valid_count,'failed':requested-valid_count,'outer_ok':int(valid_count>=max(50,np.ceil(.8*requested)))}
    for j,name in enumerate(STATNAMES):
        row[name+'_hat']=point[j]; row[name+'_lo']=q[0,j]; row[name+'_hi']=q[3,j]; row[name+'_se']=se[j]
    pos_either=(q[0,9]>0 and q[3,10]<0); neg_either=(q[3,9]<0 and q[0,10]>0)
    row.update(either_posneg=int(pos_either),either_negpos=int(neg_either),
        prespec_posneg=int(q[1,9]>0 and q[2,10]<0),
        prespec_negpos=int(q[2,9]<0 and q[1,10]>0))
    # Since IE is linear, use its 2x2 covariance to obtain exactly the same grid SDs.
    grid=np.linspace(-1.5,1.5,101); est=point[5]+point[6]*grid
    cov=np.cov(boot[:,[5,6]],rowvar=False,ddof=1)
    sd=np.sqrt(np.maximum(0.,cov[0,0]+2*grid*cov[0,1]+grid*grid*cov[1,1]))
    maxdev=np.zeros(len(boot))
    for x,s in zip(grid,sd):
        if s>0: maxdev=np.maximum(maxdev,np.abs((boot[:,5]+boot[:,6]*x-point[5]-point[6]*x)/s))
    critical=np.quantile(maxdev,.95,method='linear')
    truth=case['b']*(case['a1']+2*case['a2']*grid)
    row.update(sup_critical=critical,
        band_cover=int(np.all((truth>=est-critical*sd)&(truth<=est+critical*sd))),
        band_same_positive=int(np.min(est-critical*sd)>0),
        band_same_negative=int(np.max(est+critical*sd)<0),
        band_mean_width=float(np.mean(2*critical*sd)))
    if pos_either: classification='support_reversal_posneg'
    elif neg_either: classification='support_reversal_negpos'
    elif row['band_same_positive']: classification='support_same_sign_positive'
    elif row['band_same_negative']: classification='support_same_sign_negative'
    else: classification='uncertain'
    roots=-boot[:,0]/(2*boot[:,1]); finite=np.isfinite(roots)&(boot[:,1]!=0)&(boot[:,2]!=0)
    row['root_hat']=-point[0]/(2*point[1]) if point[1]!=0 and point[2]!=0 else np.nan
    rci=np.quantile(roots[finite],[.025,.975],method='linear') if finite.any() else [np.nan,np.nan]
    row.update(root_lo=rci[0],root_hi=rci[1],root_finite_share=float(finite.mean()),
               root_denom_clear=int(q[0,1]>0 or q[3,1]<0),root_path_clear=int(q[0,2]>0 or q[3,2]<0),classification=classification)
    return row
