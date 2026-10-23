"""Reference implementation of Umediation inference for Monte Carlo replication.

The component rectangle is projected by interval arithmetic. WCR-C and WCR-S
test linear contrasts; reversal decisions combine their one-sided p-values.
"""
import numpy as np
from scipy.stats import t as student_t
import ols
L,U=-1.5,1.5
NAMES=ols.STATNAMES

def interval_product(a,b):
    corners=np.array([a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]])
    return np.array([corners.min(),corners.max()])

def interval_linear(a,b,c):
    d=np.asarray(b)*c
    if c<0: d=d[::-1]
    return np.asarray(a)+d

def root_set(s,k,mid=0.):
    """Exact image of midpoint-s/k, excluding k=0; extended-real intervals."""
    sl,su=s; kl,ku=k
    if kl==ku==0: return [],'undefined'
    if kl>0 or ku<0:
        values=[mid-a/b for a in s for b in k]
        return [[min(values),max(values)]],'bounded'
    if sl<=0<=su:
        if sl==su==0: return [[mid,mid]],'bounded'
        if kl==0 and sl==0: return [[-np.inf,mid]],'unbounded'
        if kl==0 and su==0: return [[mid,np.inf]],'unbounded'
        if ku==0 and sl==0: return [[mid,np.inf]],'unbounded'
        if ku==0 and su==0: return [[-np.inf,mid]],'unbounded'
        return [[-np.inf,np.inf]],'all_real'
    pieces=[]
    if kl<0:
        if sl>0: pieces.append([mid-sl/kl,np.inf])
        else: pieces.append([-np.inf,mid-su/kl])
    if ku>0:
        if sl>0: pieces.append([-np.inf,mid-sl/ku])
        else: pieces.append([mid-su/ku,np.inf])
    pieces.sort(key=lambda v:v[0])
    return pieces,('disjoint' if len(pieces)>1 else 'unbounded')

def region_from_deletions(point,deletions,level=95,L=L,U=U,at=(L,0.,U)):
    g=len(deletions); mid=(L+U)/2; alpha=1-level/100
    centered=deletions-deletions.mean(axis=0)
    V=(g-1)/g*(centered.T@centered)
    def components(z):
        return np.column_stack([z[:,0]+2*mid*z[:,1],2*z[:,1],z[:,2],2*z[:,4]])
    theta=components(point[None,:])[0]; tj=components(deletions)
    c=tj-tj.mean(axis=0); vc=(g-1)/g*c.T@c
    crit=student_t.isf(alpha/8,g-1)
    ci=np.column_stack([theta-crit*np.sqrt(np.diag(vc)),theta+crit*np.sqrt(np.diag(vc))])
    s,k,b,d=ci
    a1=interval_linear(s,k,-mid); a2=k/2
    eta=interval_product(a1,b); km=interval_product(k,b)
    # c1 is the only ordinary marginal interval, outside the joint family.
    ccrit=student_t.isf(alpha/2,g-1)
    c1=point[3]+np.array([-1.,1.])*ccrit*np.sqrt(V[3,3])
    out=[a1,a2,b,c1,d/2,eta,km,d,km+d]
    def ie_at(x): return interval_product(interval_linear(s,k,x-mid),b)
    out += [ie_at(L),ie_at(U),km*(U-L)]
    out += [ie_at(x) for x in at]
    out=np.asarray(out)
    grid=np.linspace(L,U,101); band=np.asarray([ie_at(x) for x in grid])
    pos=out[9,0]>0 and out[10,1]<0
    neg=out[9,1]<0 and out[10,0]>0
    positive=np.min(band[:,0])>0; negative=np.max(band[:,1])<0
    classification=('support_reversal_posneg' if pos else 'support_reversal_negpos' if neg
       else 'support_same_sign_positive' if positive else 'support_same_sign_negative' if negative else 'uncertain')
    roots,root_type=root_set(s,k,mid)
    # Component deletion shifts, scaled by the same CV3J standard errors.
    denom=np.sqrt(np.diag(vc)); shifts=np.abs(tj-theta)
    shifts=np.divide(shifts,denom,out=np.zeros_like(shifts),where=denom>0)
    return dict(ci=out,V=V,component_ci=ci,component_v=vc,critical=crit,band=band,
                grid=grid,posneg=int(pos),negpos=int(neg),classification=classification,
                root_set=roots,root_type=root_type,influence=shifts.max(axis=1))

def wild_weights(g,B,rng,kind='webb'):
    values=np.array([-np.sqrt(1.5),-1.,-np.sqrt(.5),np.sqrt(.5),1.,np.sqrt(1.5)]) if kind=='webb' else np.array([-1.,1.])
    m=len(values); count=m**g if g<20 else B+1
    if count<=B and count<=100000:
        code=np.arange(count); W=np.empty((count,g))
        for j in range(g): W[:,j]=values[(code//m**j)%m]
        return W,True
    return values[rng.integers(m,size=(B,g))],False

class Equation:
    def __init__(self,X,y,groups):
        self.X=X; self.y=y; self.groups=groups; self.g=int(groups.max()+1)
        self.A=np.einsum('ni,nj->nij',X,X)
        self.A=np.add.reduceat(self.A,np.r_[0,np.flatnonzero(np.diff(groups))+1],axis=0)
        self.xy=np.add.reduceat(X*y[:,None],np.r_[0,np.flatnonzero(np.diff(groups))+1],axis=0)
        self.Q=self.A.sum(axis=0); self.inv=np.linalg.inv(self.Q)
        self.beta=self.inv@self.xy.sum(axis=0)
        self.scores=self.xy-np.einsum('gij,j->gi',self.A,self.beta)
        self.leave_inv=np.linalg.inv(self.Q[None,:,:]-self.A)
        self.leave_beta=np.einsum('gij,gj->gi',self.leave_inv,self.xy.sum(axis=0)-self.xy)

    def wild(self,r,W,enumerated=False,variant='classic'):
        g=self.g; v=self.inv@r; delta=float(r@self.beta)
        actual_se=np.sqrt(g/(g-1)*np.sum((self.scores@v)**2))
        if not np.isfinite(actual_se) or actual_se<=1e-14:
            raise ArithmeticError('zero or undefined original cluster standard error')
        beta_r=self.beta-v*(delta/(r@v))
        if variant=='classic': F=self.xy-np.einsum('gij,j->gi',self.A,beta_r)
        elif variant=='score':
            vr=np.einsum('gij,j->gi',self.leave_inv,r)
            br=self.leave_beta-vr*((self.leave_beta@r)/np.einsum('gi,i->g',vr,r))[:,None]
            F=self.xy-np.einsum('gij,gj->gi',self.A,br)
        else: raise ValueError(variant)
        fv=F@v; delta_star=W@F@self.inv
        Av=np.einsum('gij,j->gi',self.A,v)
        scores=W*fv[None,:]-delta_star@Av.T
        se=np.sqrt(g/(g-1)*np.sum(scores*scores,axis=1))
        ts=np.divide(W@fv,se,out=np.full(len(W),np.nan),where=se>1e-14)
        obs=delta/actual_se; invalid=~np.isfinite(ts)
        tolerance=1e-10*(1+abs(obs))
        # Degenerate draws count against rejection in both tails; never dropped.
        above=np.count_nonzero(invalid | (ts>=obs-tolerance))
        below=np.count_nonzero(invalid | (ts<=obs+tolerance))
        plus=0 if enumerated else 1
        return dict(estimate=delta,se=actual_se,t=obs,p_pos=(above+plus)/(len(W)+plus),
                    p_neg=(below+plus)/(len(W)+plus),invalid=int(invalid.sum()),tstar=ts)

def combine_reversal(results,alpha=.05):
    left,right,b=results
    branches=np.array([max(b['p_pos'],left['p_pos'],right['p_neg']),
                       max(b['p_neg'],left['p_neg'],right['p_pos']),
                       max(b['p_pos'],left['p_neg'],right['p_pos']),
                       max(b['p_neg'],left['p_pos'],right['p_neg'])])
    pp=min(1.,2*branches[:2].min()); pn=min(1.,2*branches[2:].min())
    pe=min(1.,4*branches.min())
    return dict(p_posneg=pp,p_negpos=pn,p_either=pe,prespec_posneg=int(pp<=alpha),
                prespec_negpos=int(pn<=alpha),either_posneg=int(pe<=alpha and pp<=pn),
                either_negpos=int(pe<=alpha and pn<pp),branches=branches)

def make_equations(data,st):
    x=data['x']; X=np.column_stack([x,x*x,data['z']])
    if st['g'] and st['year_fe']:
        X=np.column_stack([X]+[(data['year']==j).astype(float) for j in range(1,st['t'])])
    if not st['g']: X=np.column_stack([X,np.ones(len(x))])
    q=np.column_stack([X,data['u'],data['v']]); groups=data['firm']
    if st['g']:
        starts=np.r_[0,np.cumsum(data['sizes'])[:-1]]
        means=np.add.reduceat(q,starts,axis=0)/data['sizes'][:,None]
        q=q-means[groups]
    return q[:,:-2],q[:,-2],q[:,-1],groups

def generate(st,rep):
    if not st.get('extension'): return ols.generate(st,rep)
    rng=ols.rng_for(st,rep,0); g=st['g']
    sizes=np.full(g,st['t'])
    if st.get('imbalance'): sizes=np.array([220]+[20]*9)
    firm=np.repeat(np.arange(g),sizes); n=len(firm)
    year=np.concatenate([np.arange(k) for k in sizes])
    x=rng.uniform(-2,2,n); x[0]=-1.9; x[-1]=1.9
    if st.get('leverage'):
        # All X values remain in [-2,2]; most X^2 information lies in cluster 0.
        x[firm!=0]*=.15
        x[0]=-1.9; x[sizes[0]-1]=1.9
    z=rng.normal(size=n); innovations=rng.normal(size=(2,n))
    scale=np.exp(.35*rng.normal(size=g)-.35**2/2); uv=np.empty((2,n)); offset=0; rho=.65
    for j,k in enumerate(sizes):
        uv[:,offset]=innovations[:,offset]*scale[j]
        for t in range(1,k): uv[:,offset+t]=rho*uv[:,offset+t-1]+np.sqrt(1-rho*rho)*innovations[:,offset+t]*scale[j]
        offset+=k
    means=np.bincount(firm,weights=x)/sizes
    return dict(x=x,z=z,u=uv[0],v=uv[1],firm=firm,year=year,sizes=sizes,n=n,
                fe_m=(.4*rng.normal(size=g)+.3*means)[firm],
                fe_y=(.4*rng.normal(size=g)-.2*means)[firm])

def analyze_one(data,st,cases,B=1999,do_wild=True,do_pairs=False):
    products,p,tri=ols.sufficient(data,st); g=len(products)
    point,ok=ols.projections(products,p,tri,np.ones(g))
    leave,valid=ols.projections(products,p,tri,np.ones((g,g))-np.eye(g))
    if not ok[0] or not valid.all():
        return [dict(scenario=c['id'],method='jackknife',outer_ok=0,failed_deletions=int((~valid).sum())) for c in cases]
    X,u,v,groups=make_equations(data,st)
    wild_cache={}; pairs=None
    if do_pairs:
        weights=ols.rng_for(st,data['rep'],1).multinomial(g,np.full(g,1/g),size=B)
        bp,bv=ols.projections(products,p,tri,weights); pairs=(bp[bv],int(bv.sum()))
    if do_wild:
        W,enumerated=wild_weights(g,B,ols.rng_for(st,data['rep'],30))
        # X and u reparameterize observed M exactly, leaving the b contrast invariant.
        Z=np.column_stack([X,u])
    rows=[]
    for case in cases:
        truth=ols.evaluate(np.zeros((1,5)),case)[0]
        ph=ols.evaluate(point,case)[0]; jk=ols.evaluate(leave,case)
        region=region_from_deletions(ph,jk)
        ci=region['ci']; grid=region['grid']; band=region['band']
        truecurve=case['b']*(case['a1']+2*case['a2']*grid)
        truecomp=np.array([case['a1'],2*case['a2'],case['b'],2*case['c2']])
        component_cover=(truecomp>=region['component_ci'][:,0])&(truecomp<=region['component_ci'][:,1])
        trueroot=-case['a1']/(2*case['a2']) if case['a2']!=0 and case['b']!=0 else np.nan
        rcover=any(a<=trueroot<=b for a,b in region['root_set']) if np.isfinite(trueroot) else np.nan
        row=dict(scenario=case['id'],case=case['case'],method='jackknife',outer_ok=1,n=data['n'],G=g,
            either_posneg=region['posneg'],either_negpos=region['negpos'],
            prespec_posneg=region['posneg'],prespec_negpos=region['negpos'],
            band_cover=int(np.all((truecurve>=band[:,0])&(truecurve<=band[:,1]))),
            band_mean_width=float(np.mean(band[:,1]-band[:,0])),classification=region['classification'],
            component_joint_cover=int(component_cover.all()),component_min_cover=int(component_cover.min()),
            root_cover=rcover,root_type=region['root_type'],max_influence=float(region['influence'].max()),
            size_effective_G=float(data['n']**2/np.sum(data['sizes']**2)),failed_deletions=0)
        for j,name in enumerate(NAMES):
            row[name+'_hat']=ph[j]; row[name+'_lo']=ci[j,0]; row[name+'_hi']=ci[j,1]
            row[name+'_cover']=int(ci[j,0]<=truth[j]<=ci[j,1]); row[name+'_width']=ci[j,1]-ci[j,0]
        for j,name in enumerate(['slope_mid','curvature_xm','bpath','kappa_d']):
            row['component_'+name+'_cover']=int(component_cover[j])
        rows.append(row)
        if do_pairs:
            pb,pvalid=pairs
            s=ols.summarize(ph,ols.evaluate(pb,case),case,pvalid,B,data['n'])
            s.update(scenario=case['id'],case=case['case'],method='pairs',G=g)
            for j,name in enumerate(NAMES):
                s[name+'_cover']=int(s[name+'_lo']<=truth[j]<=s[name+'_hi'])
                s[name+'_width']=s[name+'_hi']-s[name+'_lo']
            rows.append(s)
        if do_wild:
            for variant in ['classic','score']:
                keys=[('left',case['a1'],case['a2'],variant),('right',case['a1'],case['a2'],variant),('b',case['b'],variant)]
                results=[]
                for target,key in zip(['left','right','b'],keys):
                    if key not in wild_cache:
                        if target=='b':
                            eq=Equation(Z,case['b']*u+v,groups)
                            r=np.zeros(Z.shape[1]); r[-1]=1
                        else:
                            eq=Equation(X,case['a1']*X[:,0]+case['a2']*X[:,1]+u,groups)
                            r=np.zeros(X.shape[1]); r[0]=1; r[1]=2*(L if target=='left' else U)
                        wild_cache[key]=eq.wild(r,W,enumerated,variant)
                    results.append(wild_cache[key])
                comb=combine_reversal(results); branches=comb.pop('branches')
                rows.append(dict(scenario=case['id'],case=case['case'],method='wild_'+variant,outer_ok=1,n=data['n'],G=g,
                    requested=B,wild_valid=len(W),wild_enumerated=int(enumerated),
                    wild_invalid=max(z['invalid'] for z in results),**comb,
                    p_left_positive=results[0]['p_pos'],p_left_negative=results[0]['p_neg'],
                    p_right_positive=results[1]['p_pos'],p_right_negative=results[1]['p_neg'],
                    p_b_positive=results[2]['p_pos'],p_b_negative=results[2]['p_neg']))
    return rows
