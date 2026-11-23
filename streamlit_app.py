#!/usr/bin/env python3
"""FrameLab 2D OpenSees Live. Self-contained source, Python 3.12.
Install: python -m pip install openseespy streamlit plotly numpy scipy pandas
Run GUI: python -m streamlit run FrameLab2D_OpenSees_Live.py
Run solver only: python FrameLab2D_OpenSees_Live.py --verify
OpenSees is the structural solver; Streamlit/Plotly provide the interface.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MPIR_CVAR_CH3_INTERFACE_HOSTNAME','127.0.0.1')
os.environ.setdefault('MPICH_INTERFACE_HOSTNAME','127.0.0.1')
os.environ.setdefault('HWLOC_COMPONENTS','-linux')

# ---------- Independent matrix solver ----------
"""Independent direct-stiffness and corotational calculations; no OpenSees calls."""
import numpy as np
from scipy.linalg import eigh, lu_factor, lu_solve

def ref_elt(m, e):
    i, j, E, A, I = e
    i = int(i)
    j = int(j)
    v = m['xy'][j] - m['xy'][i]
    L = np.linalg.norm(v)
    c, s = v / L
    a = E * A / L
    b = 12 * E * I / L ** 3
    cc = 6 * E * I / L ** 2
    d = 4 * E * I / L
    ee = d / 2
    k = np.array([[a, 0, 0, -a, 0, 0], [0, b, cc, 0, -b, cc], [0, cc, d, 0, -cc, ee], [-a, 0, 0, a, 0, 0], [0, -b, -cc, 0, b, -cc], [0, cc, ee, 0, -cc, d]])
    T = np.array([[c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1]])
    ids = np.r_[np.arange(3 * i, 3 * i + 3), np.arange(3 * j, 3 * j + 3)]
    return (k, T, ids, L)

def ref_linear(m, alpha=None, F=None):
    if F is None:
        F = m['F']
    if alpha is None:
        alpha = np.zeros(len(m['el']))
    K = np.zeros((len(F), len(F)))
    data = []
    for e, row in enumerate(m['el']):
        k, T, ids, L = ref_elt(m, row)
        k = k * (1 - alpha[e])
        K[np.ix_(ids, ids)] += T.T @ k @ T
        data.append((k, T, ids, L))
    free = m['free']
    u = np.zeros(len(F))
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    forces = np.array([k @ T @ u[ids] for k, T, ids, L in data])
    return dict(u=u, K=K, forces=forces, R=K @ u - F, data=data)
D = np.array([[-1.0, 0, 0, 1, 0, 0], [0, -1.0, 0, 0, 1, 0]])

def ref_corot(m, U):
    ff = np.zeros(len(U))
    K = np.zeros((len(U), len(U)))
    actions = []
    energy = 0.0
    for row in m['el']:
        i, j, E, A, I = row
        i = int(i)
        j = int(j)
        ids = np.r_[np.arange(3 * i, 3 * i + 3), np.arange(3 * j, 3 * j + 3)]
        ue = U[ids]
        v0 = m['xy'][j] - m['xy'][i]
        L = np.linalg.norm(v0)
        dv = D @ ue
        v = v0 + dv
        ell = np.linalg.norm(v)
        n = v / ell
        c, s = n
        angle = np.arctan2(v0[0] * v[1] - v0[1] * v[0], v0 @ v)
        ext = (2 * v0 @ dv + dv @ dv) / (ell + L)
        qb = np.array([ext, ue[2] - angle, ue[5] - angle])
        kb = np.array([[E * A / L, 0, 0], [0, 4 * E * I / L, 2 * E * I / L], [0, 2 * E * I / L, 4 * E * I / L]])
        B = np.zeros((3, 6))
        B[0] = n @ D
        grad = np.array([-s, c]) @ D / ell
        B[1:] = -grad
        B[1, 2] += 1
        B[2, 5] += 1
        pb = kb @ qb
        fe = B.T @ pb
        Hl = D.T @ (np.eye(2) - np.outer(n, n)) @ D / ell
        Hb = D.T @ np.array([[2 * c * s, s * s - c * c], [s * s - c * c, -2 * c * s]]) @ D / ell ** 2
        ke = B.T @ kb @ B + pb[0] * Hl - (pb[1] + pb[2]) * Hb
        ff[ids] += fe
        K[np.ix_(ids, ids)] += ke
        actions.append([-pb[0], (pb[1] + pb[2]) / ell, pb[1], pb[0], -(pb[1] + pb[2]) / ell, pb[2]])
        energy += 0.5 * qb @ pb
    return (ff, K, np.array(actions), energy)

def ref_newton(m, F):
    free = m['free']
    span = max(1, np.ptp(m['xy'][:, 0]), np.ptp(m['xy'][:, 1]))
    w = np.ones(len(free))
    w[free % 3 == 2] = 1 / span
    u = np.zeros(len(F))
    fraction = 0.0
    step = 0.125
    history = []
    while fraction < 1 - 1e-12:
        target = min(1, fraction + step)
        Ft = target * F
        den = max(1, np.linalg.norm(w * Ft[free]))
        trial = u.copy()
        ok = False
        for it in range(36):
            fi, K, _, _ = ref_corot(m, trial)
            res = Ft - fi
            rel = np.linalg.norm(w * res[free]) / den
            if rel <= 1e-09:
                ok = True
                break
            if it == 35:
                break
            try:
                np.linalg.cholesky(K[np.ix_(free, free)])
            except np.linalg.LinAlgError:
                break
            du = np.zeros(len(F))
            du[free] = np.linalg.solve(K[np.ix_(free, free)], res[free])
            a = 1.0
            accepted = False
            for _ in range(14):
                cand = trial + a * du
                fc = ref_corot(m, cand)[0]
                rc = np.linalg.norm(w * (Ft - fc)[free]) / den
                if rc < (1 - 0.0001 * a) * rel or rc <= 1e-09:
                    accepted = True
                    break
                a /= 2
            if not accepted:
                break
            trial = cand
        history.append([target, it, rel])
        if ok:
            u = trial
            fraction = target
            if it <= 6:
                step = min(0.5, 1.5 * step)
        else:
            step /= 2
            if step < 1 / 2048:
                raise ValueError('Independent Newton solve did not converge on a stable branch.')
    return (u, history)

def ref_dynamics(m, mass, zeta, frequency, dt, duration, height, ncols):
    fs = m['free']
    K = ref_linear(m)['K'][np.ix_(fs, fs)]
    md = np.array([mass / ncols * (height ** 2 / 12 if i % 3 == 2 else 1) for i in fs])
    M = np.diag(md)
    w1 = np.sqrt(eigh(K, M, eigvals_only=True)[0])
    C = 2 * zeta * w1 * M
    lu = lu_factor(K + 4 / dt ** 2 * M + 2 / dt * C)
    u = np.zeros(len(fs))
    v = u.copy()
    a = u.copy()
    states = [np.zeros(len(m['F']))]
    for k in range(1, int(np.floor(duration / dt + 1e-09)) + 1):
        pe = m['F'][fs] * np.sin(2 * np.pi * frequency * k * dt) + M @ (4 / dt ** 2 * u + 4 / dt * v + a) + C @ (2 / dt * u + v)
        un = lu_solve(lu, pe)
        an = 4 / dt ** 2 * (un - u) - 4 / dt * v - a
        vn = v + dt / 2 * (a + an)
        u, v, a = (un, vn, an)
        full = np.zeros(len(m['F']))
        full[fs] = u
        states.append(full)
    return (np.array(states), w1 / (2 * np.pi))
# ---------- OpenSees structural solver ----------
"""Live, parameter-driven OpenSees RC frame analysis. SI units internally."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MPIR_CVAR_CH3_INTERFACE_HOSTNAME','127.0.0.1')
os.environ.setdefault('MPICH_INTERFACE_HOSTNAME','127.0.0.1')
os.environ.setdefault('HWLOC_COMPONENTS','-linux')
import numpy as np,json,hashlib,uuid,time,threading
from datetime import datetime,timezone
import openseespy.opensees as ops

import builtins
LOCK=vars(builtins).setdefault('_framelab_opensees_lock',threading.RLock()) # OpenSees owns one global domain per Python process.
DEFAULT=dict(stories=2,bays=2,width=6.,height=3.5,Ec_GPa=30.,Es_GPa=200.,b_mm=400.,h_mm=400.,bars_per_face=2,diameter_mm=16.,cover_mm=50.,roof_kN=100.,mass_tonne=80.,damping_pct=5.,frequency_Hz=1.2,dt=.02,duration=6.,gravity_kN=1000.,gravity_factor=5.,segments=4,scope='All members',custom_members='1,6',Ec_loss_pct=10.,As_loss_pct=10.,run_dynamic=True,run_nonlinear=True,compare_healthy=True)
def fingerprint(c):return hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()
def section(c,ec_loss=0.,as_loss=0.):
    b=c['b_mm']/1000;h=c['h_mm']/1000;d=c['diameter_mm']/1000;a=c['cover_mm']/1000;n=int(c['bars_per_face']);Ec=c['Ec_GPa']*1e9;Es=c['Es_GPa']*1e9
    if not(Ec>0 and Es>0 and b>0 and h>0 and n>=1 and d>0 and a>=d/2 and h-2*a>=d and b>=2*a):raise ValueError('Invalid section or bar-centre cover: reinforcement must fit inside the rectangle.')
    if n>1 and (b-2*a)/(n-1)<d:raise ValueError('Bars overlap within a face.')
    Ab=np.pi*d*d/4;Ib=np.pi*d**4/64;z=h/2-a;As0=2*n*Ab;Is0=2*n*(Ib+Ab*z*z);Ac=b*h-As0;Ic=b*h**3/12-Is0
    if Ac<=0 or Ic<=0:raise ValueError('Reinforcement leaves nonpositive concrete area or inertia.')
    ret=1-as_loss;As=ret*As0;Is=2*n*(ret**2*Ib+ret*Ab*z*z);EA=Ec*(1-ec_loss)*Ac+Es*As;EI=Ec*(1-ec_loss)*Ic+Es*Is
    return dict(Ag=b*h,Ig=b*h**3/12,As0=As0,As=As,Ac=Ac,Ic=Ic,EA=EA,EI=EI,Aeq=EA/Ec,Ieq=EI/Ec)
def make_model(c,healthy=False):
    ns=int(c['stories']);nb=int(c['bays']);nc=nb+1
    if not(1<=ns<=4 and 1<=nb<=3 and c['width']>0 and c['height']>0):raise ValueError('Use 1-4 stories, 1-3 bays and positive dimensions.')
    xy=np.array([(col*c['width'],f*c['height'])for f in range(ns+1)for col in range(nc)])
    edges=[];meta=[]
    for f in range(1,ns+1):
        for col in range(nc):edges.append(((f-1)*nc+col,f*nc+col));meta.append(('Column',f,col))
        for col in range(nb):edges.append((f*nc+col,f*nc+col+1));meta.append(('Beam',f,col))
    scope=c['scope'];selected=[]
    if scope=='All members':selected=list(range(len(edges)))
    elif scope=='Bottom columns':selected=[i for i,v in enumerate(meta)if v[0]=='Column'and v[1]==1]
    elif scope=='Left columns':selected=[i for i,v in enumerate(meta)if v[0]=='Column'and v[2]==0]
    elif scope=='Custom member IDs':
        try:selected=[int(x.strip())-1 for x in c['custom_members'].split(',')if x.strip()]
        except ValueError:raise ValueError('Custom member IDs must be comma-separated integers.')
        if any(x<0 or x>=len(edges)for x in selected):raise ValueError('A custom member ID is outside this frame. Inspect the member map.')
    if healthy:selected=[]
    if not(0<=c['Ec_loss_pct']<=80 and 0<=c['As_loss_pct']<=80):raise ValueError('Loss percentages must be between 0 and 80.')
    el=[];audit=[]
    for e,(i,j)in enumerate(edges):
        ce=c['Ec_loss_pct']/100 if e in selected else 0;ar=c['As_loss_pct']/100 if e in selected else 0;p=section(c,ce,ar)
        el.append((i,j,c['Ec_GPa']*1e9,p['Aeq'],p['Ieq']))
        audit.append(dict(Element=e+1,NodeI=i+1,NodeJ=j+1,Type=meta[e][0],Story=meta[e][1],Ec_loss_pct=ce*100,As_loss_pct=ar*100,EA_N=p['EA'],EI_Nm2=p['EI'],Aeq_m2=p['Aeq'],Ieq_m4=p['Ieq']))
    F=np.zeros(3*len(xy));roof=np.arange(ns*nc,(ns+1)*nc);F[3*roof]=c['roof_kN']*1000/nc
    return dict(xy=xy,el=el,F=F,free=np.arange(3*nc,len(F)),roof=roof,ncols=nc,stories=ns,audit=audit)
def refine(m,n):
    xy=m['xy'].tolist();el=[]
    for i,j,E,A,I in m['el']:
        ids=[i]
        for k in range(1,n):ids.append(len(xy));xy.append((m['xy'][i]+k/n*(m['xy'][j]-m['xy'][i])).tolist())
        ids.append(j);el.extend((ids[k],ids[k+1],E,A,I)for k in range(n))
    F=np.zeros(3*len(xy));F[:len(m['F'])]=m['F']
    return {**m,'xy':np.array(xy),'el':el,'F':F,'free':np.arange(3*m['ncols'],len(F))}
def gravity(m,c):
    F=np.zeros(len(m['F']));nc=m['ncols']
    for f in range(1,m['stories']+1):F[3*np.arange(f*nc,(f+1)*nc)+1]=-c['gravity_kN']*1000/nc*c['gravity_factor']
    return F
def build(m,nonlinear=False):
    ops.wipe();ops.model('basic','-ndm',2,'-ndf',3)
    for i,(x,y)in enumerate(m['xy'],1):ops.node(i,float(x),float(y))
    for i in range(m['ncols']):ops.fix(i+1,1,1,1)
    ops.geomTransf('Corotational'if nonlinear else'Linear',1)
    for e,(i,j,E,A,I)in enumerate(m['el'],1):ops.element('elasticBeamColumn',e,int(i)+1,int(j)+1,float(A),float(E),float(I),1)
def setup():ops.constraints('Plain');ops.numberer('RCM');ops.system('BandGeneral')
def op_static(m,F=None,nonlinear=False):
    build(m,nonlinear);F=m['F']if F is None else F
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1)
    for i in range(len(m['xy'])):
        f=F[3*i:3*i+3]
        if np.any(f):ops.load(i+1,*map(float,f))
    setup();ops.test('NormUnbalance',1e-5,80);ops.algorithm('Newton'if nonlinear else'Linear');ops.integrator('LoadControl',1/32 if nonlinear else 1);ops.analysis('Static')
    if ops.analyze(32 if nonlinear else 1)!=0:raise ValueError('OpenSees could not converge. Reduce load/gravity or refine the stability mesh. No unconverged result is accepted.')
    u=np.array([ops.nodeDisp(i+1,d)for i in range(len(m['xy']))for d in [1,2,3]])
    ops.reactions();R=np.array([ops.nodeReaction(i+1,d)for i in range(len(m['xy']))for d in [1,2,3]])
    forces=np.array([ops.eleResponse(e+1,'localForce')for e in range(len(m['el']))])
    return u,R,forces
def eigen(m,c):
    build(m);mass=c['mass_tonne']*1000/m['ncols']
    for i in range(m['ncols'],len(m['xy'])):ops.mass(i+1,mass,mass,mass*c['height']**2/12)
    vals=np.array(ops.eigen('-fullGenLapack',min(6,len(m['free']))));freq=np.sqrt(vals)/2/np.pi
    phi=np.array([ops.nodeEigenvector(i+1,1,d)for i in range(len(m['xy']))for d in [1,2,3]])
    return freq,phi

def op_dynamic(m,c,w1):
    build(m);mass=c['mass_tonne']*1000/m['ncols'];dt=c['dt'];nt=int(np.floor(c['duration']/dt+1e-9))
    for i in range(m['ncols'],len(m['xy'])):ops.mass(i+1,mass,mass,mass*c['height']**2/12)
    ops.rayleigh(2*c['damping_pct']/100*w1,0,0,0)
    values=np.sin(2*np.pi*c['frequency_Hz']*np.arange(nt+2)*dt)
    ops.timeSeries('Path',1,'-dt',dt,'-values',*map(float,values));ops.pattern('Plain',1,1)
    for i in range(len(m['xy'])):
        f=m['F'][3*i:3*i+3]
        if np.any(f):ops.load(i+1,*map(float,f))
    setup();ops.algorithm('Linear');ops.integrator('Newmark',.5,.25);ops.analysis('Transient');states=[np.zeros(len(m['F']))]
    for k in range(nt):
        if ops.analyze(1,dt)!=0:raise ValueError('OpenSees dynamic analysis failed.')
        states.append(np.array([ops.nodeDisp(i+1,d)for i in range(len(m['xy']))for d in [1,2,3]]))
    return np.arange(nt+1)*dt,np.array(states)
def rel(a,b):return float(np.linalg.norm(np.array(a)-b)/max(np.linalg.norm(a),1e-12))
def mesh_dict(m):return dict(nodes=m['xy'].tolist(),elements=[[int(e[0]),int(e[1])]for e in m['el']],ncols=m['ncols'],roof=m['roof'].tolist())
def analyze_case(m,c):
    rows=[];static=ref_linear(m);u,R,actions=op_static(m);rho=rel(static['u'],u);rows.append(dict(Check='Static displacement vector',RelativeError=rho,Tolerance=1e-7,Pass=rho<1e-7))
    res=static['K']@u-m['F'];free_res=np.linalg.norm(res[m['free']])/max(1,np.linalg.norm(m['F']))
    balance=(R+m['F']).reshape(-1,3);fx=float(np.sum(balance[:,0]));fy=float(np.sum(balance[:,1]));moment=float(np.sum(balance[:,2]+m['xy'][:,0]*balance[:,1]-m['xy'][:,1]*balance[:,0]))
    force_scale=max(1,np.linalg.norm(m['F']));span=max(1,np.ptp(m['xy'][:,0]),np.ptp(m['xy'][:,1]));eq=max(abs(fx),abs(fy),abs(moment)/span)/force_scale
    rows.append(dict(Check='Global force/moment equilibrium',RelativeError=eq,Tolerance=1e-8,Pass=eq<1e-8))
    freq,phi=eigen(m,c);fs=m['free'];md=np.array([c['mass_tonne']*1000/m['ncols']*(c['height']**2/12 if i%3==2 else 1)for i in fs]);K=static['K'][np.ix_(fs,fs)]
    from scipy.linalg import eigh
    wref=np.sqrt(eigh(K,np.diag(md),eigvals_only=True)[0]);ef=abs(freq[0]-wref/2/np.pi)/(wref/2/np.pi)
    rows.append(dict(Check='First natural frequency',RelativeError=float(ef),Tolerance=1e-7,Pass=bool(ef<1e-7)))
    modal_res=K@phi[fs]-(2*np.pi*freq[0])**2*md*phi[fs];modal_err=np.linalg.norm(modal_res)/max(1e-12,(np.linalg.norm(K)+((2*np.pi*freq[0])**2)*np.linalg.norm(md))*np.linalg.norm(phi[fs]))
    rows.append(dict(Check='OpenSees modal equation residual',RelativeError=float(modal_err),Tolerance=1e-8,Pass=bool(modal_err<1e-8)))
    out=dict(mesh=mesh_dict(m),sections=m['audit'],u=u.tolist(),reactions=R.tolist(),local_end_forces=actions.tolist(),roof_static_mm=float(np.mean(u[3*m['roof']])*1000),freq_Hz=freq.tolist(),mode1=phi.tolist(),equilibrium=dict(sum_Fx_N=fx,sum_Fy_N=fy,sum_M_Nm=moment,free_residual=free_res),checks=rows)
    if c['run_dynamic']:
        t,ud=op_dynamic(m,c,2*np.pi*freq[0]);ur,fr=ref_dynamics(m,c['mass_tonne']*1000,c['damping_pct']/100,c['frequency_Hz'],c['dt'],c['duration'],c['height'],m['ncols'])
        err=rel(ur,ud);rows.append(dict(Check='Full Newmark displacement history',RelativeError=err,Tolerance=1e-7,Pass=err<1e-7))
        floors=[]
        for f in range(1,m['stories']+1):floors.append(np.mean(ud[:,3*np.arange(f*m['ncols'],(f+1)*m['ncols'])],axis=1)*1000)
        out['dynamic']=dict(time_s=t.tolist(),floor_mm=np.array(floors).tolist(),roof_peak_mm=float(np.max(abs(floors[-1]))),u_history=ud.tolist())
    if c['run_nonlinear']:
        try:
            mesh=refine(m,int(c['segments']));F=mesh['F']+gravity(mesh,c);un,rn,an=op_static(mesh,F,True);unref,hist=ref_newton(mesh,F);en=rel(unref,un)
            rows.append(dict(Check='Corotational Newton displacement',RelativeError=en,Tolerance=1e-7,Pass=en<1e-7))
            resn=ref_corot(mesh,un)[0]-F;rr=np.linalg.norm(resn[mesh['free']])/max(1,np.linalg.norm(F[mesh['free']]))
            rows.append(dict(Check='Nonlinear free-DOF equilibrium',RelativeError=float(rr),Tolerance=1e-8,Pass=bool(rr<1e-8)))
            first=ref_linear(mesh,F=F)['u'];nominal=max(abs(an[e,3])/mesh['el'][e][3]for e in range(len(an)))/1e6
            out['nonlinear']=dict(mesh=mesh_dict(mesh),u=un.tolist(),roof_mm=float(np.mean(un[3*m['roof']])*1000),first_order_combined_roof_mm=float(np.mean(first[3*m['roof']])*1000),nominal_N_over_Aeq_MPa=float(nominal))
        except Exception as exc:
            out['nonlinear']=dict(error=str(exc));rows.append(dict(Check='Nonlinear solution completed',RelativeError=None,Tolerance=1e-7,Pass=False))
    out['all_checks_pass']=all(x['Pass']for x in rows);return out

def run(c):
    c={**DEFAULT,**c}
    if c['mass_tonne']<=0 or not(0<=c['damping_pct']<=30):raise ValueError('Mass must be positive; damping 0-30%.')
    if c['dt']<=0 or c['duration']<c['dt'] or round(c['duration']/c['dt'])>1500:raise ValueError('Use positive duration/dt and at most 1500 time steps.')
    if c['gravity_kN']<0 or c['gravity_factor']<0:raise ValueError('Gravity inputs must be nonnegative.')
    model=make_model(c);n=int(c['segments'])
    if n<1 or len(model['el'])*n>160:raise ValueError('Use 1 or more segments and at most 160 total stability elements.')
    start=time.perf_counter()
    with LOCK:
        try:
            active=analyze_case(model,c);healthy=analyze_case(make_model(c,True),c)if c['compare_healthy']else None;version=ops.version()
        finally:ops.wipe()
    result=dict(run_id=uuid.uuid4().hex[:12],timestamp_utc=datetime.now(timezone.utc).isoformat(),elapsed_s=time.perf_counter()-start,input_hash=fingerprint(c),inputs=c,opensees_version=version,active=active,healthy=healthy,matlab_runtime_compared=False)

    return json.loads(json.dumps(result,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))

from types import SimpleNamespace
engine=SimpleNamespace(DEFAULT=DEFAULT,run=run,fingerprint=fingerprint,section=section)

# ---------- Comparison and visualization helpers ----------
def is_screenshot_preset(c):
    keys=['stories','bays','width','height','Ec_GPa','Es_GPa','b_mm','h_mm','bars_per_face','diameter_mm','cover_mm','roof_kN','mass_tonne','damping_pct','frequency_Hz','dt','duration','gravity_kN','gravity_factor','segments','scope','Ec_loss_pct','As_loss_pct']
    return all(c[k]==DEFAULT[k]for k in keys)

def comparison_rows(r):
    c=r['inputs'];preset=is_screenshot_preset(c);records=[]
    for key,label in [('healthy','Healthy RC'),('active','Deteriorated RC' if c['scope']=='All members'and c['Ec_loss_pct']==10 and c['As_loss_pct']==10 else 'Active RC')]:
        a=r.get(key)
        if a is None:continue
        m=make_model(c,healthy=key=='healthy');ind=ref_linear(m);u=np.array(a['u']);K=ind['K'];fs=m['free'];md=np.array([c['mass_tonne']*1000/m['ncols']*(c['height']**2/12 if i%3==2 else 1) for i in fs]);w=np.sqrt(eigh(K[np.ix_(fs,fs)],np.diag(md),eigvals_only=True)[0]);cases=[]
        cases.append(('Static roof X','mm',a['roof_static_mm'],float(np.mean(ind['u'][3*m['roof']])*1000),8.6465 if key=='healthy'else 9.6072))
        nc=m['ncols'];dd=[]
        for level in range(1,m['stories']+1):
            up=3*np.arange(level*nc,(level+1)*nc);lo=3*np.arange((level-1)*nc,level*nc);dd.extend(abs(u[up]-u[lo])/c['height']*100)
        dri=max(dd);ddi=[]
        for level in range(1,m['stories']+1):
            up=3*np.arange(level*nc,(level+1)*nc);lo=3*np.arange((level-1)*nc,level*nc);ddi.extend(abs(ind['u'][up]-ind['u'][lo])/c['height']*100)
        cases.append(('Maximum column drift','%',float(dri),float(max(ddi)),.1431 if key=='healthy'else .1590))
        cases.append(('First natural frequency','Hz',a['freq_Hz'][0],float(w/(2*np.pi)),None if key=='healthy'else 1.64))
        if 'dynamic'in a:
            ud,_=ref_dynamics(m,c['mass_tonne']*1000,c['damping_pct']/100,c['frequency_Hz'],c['dt'],c['duration'],c['height'],nc)
            cases.append(('Peak roof response','mm',a['dynamic']['roof_peak_mm'],float(np.max(abs(np.mean(ud[:,3*m['roof']],axis=1)))*1000),None if key=='healthy'else 28.498))
        if 'nonlinear'in a and 'error'not in a['nonlinear']:
            mesh=refine(m,int(c['segments']));un,_=ref_newton(mesh,mesh['F']+gravity(mesh,c));cases.append(('P-Delta roof X','mm',a['nonlinear']['roof_mm'],float(np.mean(un[3*m['roof']])*1000),9.5084 if key=='healthy'else 10.6837))
        for quantity,unit,val,calc,gui in cases:
            records.append({'Model':label,'Result':quantity,'Unit':unit,'OpenSees':val,'Independent':calc,'Solver difference %':100*abs(val-calc)/max(abs(calc),1e-14),'MATLAB screenshot':gui if preset else None,'Screenshot difference %':100*abs(val-gui)/max(abs(val),1e-14)if preset and gui is not None else None})
    return records

def dynamic_animation(a,c):
    mesh=a['mesh'];xy=np.array(mesh['nodes']);history=np.array(a['dynamic']['u_history']);times=np.array(a['dynamic']['time_s']);span=max(c['bays']*c['width'],c['stories']*c['height']);mag=np.max(np.hypot(history[:,::3],history[:,1::3]));sc=.12*span/max(mag,1e-14)
    def segments(points):
        xx=[];yy=[]
        for i,j in mesh['elements']:xx.extend([points[i,0],points[j,0],None]);yy.extend([points[i,1],points[j,1],None])
        return xx,yy
    bx,by=segments(xy);fig=go.Figure();fig.add_trace(go.Scatter(x=bx,y=by,mode='lines',line=dict(color='#8b98a5',dash='dash'),name='Undeformed'));fig.add_trace(go.Scatter(x=bx,y=by,mode='lines',line=dict(color='#e76d32',width=4),name='Dynamic (scaled)'))
    fig.add_trace(go.Scatter(x=xy[:mesh['ncols'],0],y=xy[:mesh['ncols'],1],mode='markers+text',text=['FIXED']*mesh['ncols'],textposition='bottom center',marker=dict(symbol='square',size=12,color='#18354a'),name='Fixed bases'))
    picks=np.unique(np.linspace(0,len(times)-1,min(101,len(times))).astype(int));frames=[]
    for k in picks:
        pt=xy+sc*history[k].reshape(-1,3)[:,:2];xx,yy=segments(pt);frames.append(go.Frame(name=str(k),data=[go.Scatter(x=xx,y=yy)],traces=[1],layout=go.Layout(title=f'OpenSees response | t = {times[k]:.2f} s | display ×{sc:.1f}')))
    fig.frames=frames
    fig.update_layout(height=480,margin=dict(l=20,r=20,t=50,b=90),title='OpenSees response | use Play',xaxis=dict(title='X (m)',range=[-.18*span,c['bays']*c['width']+.18*span]),yaxis=dict(title='Y (m)',range=[-.18*span,c['stories']*c['height']+.18*span],scaleanchor='x',scaleratio=1),updatemenus=[dict(type='buttons',buttons=[dict(label='Play',method='animate',args=[None,dict(frame=dict(duration=60,redraw=True),transition=dict(duration=0),fromcurrent=True)]),dict(label='Pause',method='animate',args=[[None],dict(frame=dict(duration=0,redraw=False),mode='immediate')])],direction='left',x=0,y=-.10)],sliders=[dict(y=-.03,currentvalue=dict(prefix='Frame '),steps=[dict(label=f'{times[k]:.1f}s',method='animate',args=[[str(k)],dict(mode='immediate',frame=dict(duration=0,redraw=True),transition=dict(duration=0))])for k in picks])])
    return fig

def export_html(result,path):
    import html
    c=result['inputs'];parts=[]
    for key,label in [('healthy','Healthy RC'),('active','Deteriorated RC')]:
        a=result.get(key)
        if a is None:continue
        parts.append('<h2>'+label+'</h2>')
        parts.append(draw(a,'Static deformation',100,True,c).to_html(full_html=False,include_plotlyjs=True if not parts[:-1]else False))
        parts.append(dynamic_animation(a,c).to_html(full_html=False,include_plotlyjs=False))
        f=go.Figure()
        for j,y in enumerate(a['dynamic']['floor_mm']):f.add_trace(go.Scatter(x=a['dynamic']['time_s'],y=y,name=f'Floor {j+1}'))
        f.update_layout(title=label+' | dynamic floor histories',xaxis_title='Time (s)',yaxis_title='Mean floor X (mm)');parts.append(f.to_html(full_html=False,include_plotlyjs=False))
    from pathlib import Path
    Path(path).write_text('<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>FrameLab OpenSees actual results</title><style>body{font-family:Georgia,serif;margin:24px;color:#143047}table{border-collapse:collapse;font-size:13px}th,td{padding:8px;border:1px solid #ccd5dd}pre{background:#f4f7fa;padding:16px;white-space:pre-wrap}</style></head><body><h1>FrameLab 2D OpenSees results</h1><p>Actual saved OpenSees run '+html.escape(result['run_id'])+' | version '+html.escape(result['opensees_version'])+'</p><p>This file replays saved results and animations. For fresh analysis with editable inputs, run FrameLab2D_OpenSees_Live.py with Streamlit.</p>'+pd.DataFrame(comparison_rows(result)).to_html(index=False,float_format=lambda x:f'{x:.7g}')+'<h2>Inputs actually used</h2><pre>'+html.escape(json.dumps(c,indent=2))+'</pre>'+''.join(parts)+'</body></html>',encoding='utf-8')

# ---------- Streamlit live interface ----------
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,io,csv,hashlib
import numpy as np,pandas as pd
import streamlit as st
import plotly.graph_objects as go

if '--verify' not in __import__('sys').argv:
    st.set_page_config(page_title='OpenSees Live | FrameLab 2D',page_icon='🏗️',layout='wide')

def draw(result,view,scale,labels,c):
    data=result;mesh=data['mesh'];u=np.array(data['u']);isforce=view in ['Bending moment','Axial force']
    if view=='P-Delta deformation':
        if 'nonlinear'not in data or 'error'in data['nonlinear']:return None
        mesh=data['nonlinear']['mesh'];u=np.array(data['nonlinear']['u'])
    if view=='First vibration mode':
        u=np.array(data['mode1']);amp=max(np.max(abs(u[::3])),np.max(abs(u[1::3])),1e-15);scale=.12*max(c['bays']*c['width'],c['stories']*c['height'])/amp
    if view=='Model and supports':scale=0
    xy=np.array(mesh['nodes']);fig=go.Figure();els=mesh['elements'];nc=mesh['ncols'];norig=(c['stories']+1)*nc
    for e,(i,j)in enumerate(els):
        v=xy[j]-xy[i];L=np.linalg.norm(v);cc,ss=v/L;ids=[3*i,3*i+1,3*i+2,3*j,3*j+1,3*j+2];ue=u[ids];a=cc*ue[0]+ss*ue[1];b=cc*ue[3]+ss*ue[4];vi=-ss*ue[0]+cc*ue[1];vj=-ss*ue[3]+cc*ue[4]
        xi=np.linspace(0,1,21);axial=a*(1-xi)+b*xi;vv=(1-3*xi**2+2*xi**3)*vi+L*(xi-2*xi**2+xi**3)*ue[2]+(3*xi**2-2*xi**3)*vj+L*(-xi**2+xi**3)*ue[5]
        pg=xy[i]+np.c_[cc*(xi*L+scale*axial)-ss*scale*vv,ss*(xi*L+scale*axial)+cc*scale*vv]
        fig.add_trace(go.Scatter(x=xy[[i,j],0],y=xy[[i,j],1],mode='lines',line=dict(color='#a6b4bf',dash='dash',width=1.5),name='Original',showlegend=e==0,hoverinfo='skip'))
        if not isforce:fig.add_trace(go.Scatter(x=pg[:,0],y=pg[:,1],mode='lines',line=dict(color='#c25132'if scale else'#173d59',width=3),name='Deformed'if scale else'Model',showlegend=e==0,hoverinfo='skip'))
    if isforce:
        actions=np.array(data['local_end_forces']);mx=max(abs(actions[:,[2,5]]).ravel())if view=='Bending moment'else max(abs(actions[:,[0,3]]).ravel());sf=.10*max(c['bays']*c['width'],c['stories']*c['height'])/max(mx,1e-15)
        for e,(i,j)in enumerate(els):
            v=xy[j]-xy[i];normal=np.array([-v[1],v[0]])/np.linalg.norm(v);xi=np.linspace(0,1,21);val=(-actions[e,2])*(1-xi)+actions[e,5]*xi if view=='Bending moment'else np.repeat(actions[e,3],21)
            pg=xy[i]+xi[:,None]*v+sf*val[:,None]*normal
            fig.add_trace(go.Scatter(x=pg[:,0],y=pg[:,1],mode='lines',line=dict(color='#c25132',width=2),customdata=val/1000,hovertemplate=('M: %{customdata:.4f} kN m'if view=='Bending moment'else'N: %{customdata:.4f} kN')+'<extra>E'+str(e+1)+'</extra>',name=view,showlegend=e==0))
    for n in range(nc):
        x,y=xy[n];fig.add_shape(type='line',x0=x-.22,x1=x+.22,y0=y,y1=y,line=dict(width=5,color='#173d59'))
        for h in np.linspace(x-.2,x+.2,5):fig.add_shape(type='line',x0=h-.1,x1=h,y0=y-.16,y1=y,line=dict(width=1,color='#173d59'))
    if view=='Model and supports' or isforce:
        showxy=xy[:norig];custom=np.array(data['u']).reshape(-1,3)[:norig]
    else:showxy=xy[:norig]+scale*u.reshape(-1,3)[:norig,:2];custom=u.reshape(-1,3)[:norig]
    fig.add_trace(go.Scatter(x=showxy[:,0],y=showxy[:,1],mode='markers+text'if labels else'markers',text=[f'N{i+1}'for i in range(norig)],textposition='top right',marker=dict(size=6,color='#173d59'),customdata=custom,hovertemplate='Ux %{customdata[0]:.8f} m<br>Uy %{customdata[1]:.8f} m<br>Rz %{customdata[2]:.8f} rad<extra></extra>',name='Nodes',showlegend=False))
    if view=='Model and supports' and c['roof_kN']!=0:
        for node in mesh['roof']:
            x,y=xy[node];fig.add_annotation(x=x+.7*np.sign(c['roof_kN']),y=y+.5,ax=x,ay=y+.5,xref='x',yref='y',axref='x',ayref='y',text='',showarrow=True,arrowhead=2,arrowcolor='#b64f2b')
    suffix='arbitrary mode amplitude'if view=='First vibration mode'else'diagram offsets scaled, not displacement'if isforce else f'displacement display scale x{scale:g}'
    fig.update_layout(height=510,margin=dict(l=20,r=20,t=45,b=20),title=view+' | '+suffix,paper_bgcolor='white',plot_bgcolor='#f8fbfd',legend=dict(orientation='h',y=1.05),xaxis_title='X (m)',yaxis_title='Y (m)')
    fig.update_yaxes(scaleanchor='x',scaleratio=1);return fig

def main():
    st.title('OpenSees Live Analysis')
    st.caption('FrameLab 2D independent RC model | actual OpenSees solver on every Run | inputs, results and verification together')
    st.info('Screenshot preset loaded: 2 stories × 2 bays; RC 400 mm; losses 10% / 10% on all members. Run solves both healthy and deteriorated cases. Gravity factor = 5.')
    c=dict(engine.DEFAULT)
    with st.sidebar:
        st.header('Analysis inputs')
        with st.expander('Geometry and supports',expanded=True):
            c['stories']=st.number_input('Stories',1,4,2,key='stories');c['bays']=st.number_input('Bays',1,3,2,key='bays')
            c['width']=st.number_input('Bay width (m)',1.,15.,6.,key='width');c['height']=st.number_input('Story height (m)',1.,8.,3.5,key='height')
            st.caption('Fixed bases: Ux = Uy = Rz = 0. Three DOFs per joint; no rigid-floor constraint.')
        with st.expander('RC section and reinforcement'):
            c['Ec_GPa']=st.number_input('Concrete Ec (GPa)',1.,100.,30.,key='Ec_GPa');c['Es_GPa']=st.number_input('Steel Es (GPa)',1.,300.,200.,key='Es_GPa')
            c['b_mm']=st.number_input('Gross section width b (mm)',100.,1500.,400.,key='b_mm');c['h_mm']=st.number_input('Gross section depth h (mm)',100.,1500.,400.,key='h_mm')
            c['bars_per_face']=st.number_input('Bars per opposite face',1,12,2,key='bars_per_face');c['diameter_mm']=st.number_input('Bar diameter (mm)',4.,50.,16.,key='diameter_mm');c['cover_mm']=st.number_input('Bar-centre cover (mm)',5.,200.,50.,key='cover_mm')
        with st.expander('Deterioration assumptions',expanded=True):
            scopes=['No loss','Left columns','Bottom columns','All members','Custom member IDs'];c['scope']=st.selectbox('Members receiving losses',scopes,index=3,key='scope')
            c['custom_members']=st.text_input('Custom IDs, comma-separated','1,6',key='custom_members',disabled=c['scope']!='Custom member IDs')
            c['Ec_loss_pct']=st.number_input('Concrete modulus loss (%)',0.,80.,10.,key='Ec_loss_pct');c['As_loss_pct']=st.number_input('Steel area loss (%)',0.,80.,10.,key='As_loss_pct')
            st.caption('Losses are assumptions. A chloride/sulfate label does not calculate deterioration.')
        with st.expander('Lateral load and dynamics',expanded=True):
            c['roof_kN']=st.number_input('Total lateral roof force (kN)',-500.,500.,100.,key='roof_kN')
            c['run_dynamic']=st.checkbox('Run dynamic response',True,key='run_dynamic');c['mass_tonne']=st.number_input('Mass per floor (tonnes)',1.,1000.,80.,key='mass_tonne');c['damping_pct']=st.number_input('First-mode damping (%)',0.,30.,5.,key='damping_pct');c['frequency_Hz']=st.number_input('Applied harmonic force frequency (Hz)',.05,10.,1.2,key='frequency_Hz')
            c['dt']=st.number_input('Time step (s)',.005,.1,.02,step=.005,format='%.3f',key='dt');c['duration']=st.number_input('Duration (s)',.1,20.,6.,key='duration')
        with st.expander('P-Delta inputs'):
            c['run_nonlinear']=st.checkbox('Run corotational Newton analysis',True,key='run_nonlinear');c['gravity_kN']=st.number_input('Gravity per floor (kN)',0.,10000.,1000.,key='gravity_kN');c['gravity_factor']=st.number_input('Gravity multiplier',0.,50.,5.,key='gravity_factor');c['segments']=st.number_input('Segments per member',1,8,4,key='segments')
        c['compare_healthy']=st.checkbox('Also solve matching healthy RC baseline',True,key='compare_healthy')
        run=st.button('RUN LIVE OPENSEES ANALYSIS',type='primary',use_container_width=True,key='run_analysis')
    if run:
        try:
            with st.spinner('Building models and running fresh OpenSees + independent matrix calculations...'):
                result=engine.run(c)
            st.session_state['result']=result
        except Exception as exc:
            st.session_state.pop('result',None);st.error(str(exc));return
    r=st.session_state.get('result')
    if r is None:
        st.info('Enter or review the inputs, then click RUN LIVE OPENSEES ANALYSIS. No result is shown before a real solve.')
        st.subheader('Current input draft');st.json(c);return
    if r['input_hash']!=engine.fingerprint(c):
        st.warning('Inputs changed. The previous results are hidden. Click RUN LIVE OPENSEES ANALYSIS to solve these new inputs.')
        st.subheader('New input draft');st.json(c);return
    st.caption(f"OpenSees {r['opensees_version']} | Run ID {r['run_id']} | UTC {r['timestamp_utc']} | Solve time {r['elapsed_s']:.2f} s")
    a=r['active'];labels={'Active RC':a}
    if r['healthy']is not None:labels['Healthy RC baseline']=r['healthy']
    which=st.radio('Display results for',list(labels),horizontal=True,key='case_display');a=labels[which]
    if a['all_checks_pass']:st.success('PASS: this run agrees with the independent matrix calculation and the listed equilibrium/residual checks.')
    else:st.error('Some requested checks did not pass or a solve was incomplete. Inspect verification details; do not accept this run as verified.')
    col=st.columns(4);col[0].metric('Static mean roof X',f"{a['roof_static_mm']:.6f} mm");col[1].metric('First natural frequency',f"{a['freq_Hz'][0]:.6f} Hz")
    col[2].metric('Dynamic roof peak',f"{a['dynamic']['roof_peak_mm']:.6f} mm"if'dynamic'in a else'Not run');col[3].metric('P-Delta mean roof',f"{a['nonlinear']['roof_mm']:.6f} mm"if'nonlinear'in a and'error'not in a['nonlinear']else'Not completed')
    st.subheader('Live OpenSees versus independent calculation')
    st.dataframe(pd.DataFrame(comparison_rows(r)),hide_index=True,use_container_width=True)
    st.caption('The stored MATLAB screenshot column is a reference, not a new MATLAB solve. It is shown only for the exact matching preset. Sensor diagnosis and placement are not OpenSees built-in verifications in this companion.')
    tabs=st.tabs(['Model and diagrams','Inputs actually used','Verification','Tables and export'])
    with tabs[0]:
        view=st.selectbox('Diagram',['Model and supports','Static deformation','First vibration mode','P-Delta deformation','Bending moment','Axial force'],index=1,key='diagram')
        scale=st.slider('Deformation magnification (display only)',1,300,100,key='magnification');labels_on=st.checkbox('Show original node IDs',True,key='labels')
        fig=draw(a,view,scale,labels_on,c)
        if fig is not None:st.plotly_chart(fig,use_container_width=True)
        elif'nonlinear'in a and'error'in a['nonlinear']:st.error(a['nonlinear']['error'])
        else:st.info('Enable and run P-Delta to display this diagram.')
        if view=='First vibration mode':st.caption('This is a vibration eigenmode with arbitrary amplitude, not a buckling mode or a physical displacement prediction.')
        if'dynamic'in a:
            if st.checkbox('Show animated displacement playback',True,key='show_animation'):
                st.plotly_chart(dynamic_animation(a,c),use_container_width=True)
                st.caption('Play replays the freshly computed time history. Use Run Analysis to calculate again after changing inputs.')
            d=a['dynamic'];fig=go.Figure()
            for j,y in enumerate(d['floor_mm']):fig.add_trace(go.Scatter(x=d['time_s'],y=y,name=f'Floor {j+1}'))
            fig.update_layout(title='Actual OpenSees harmonic-force response',xaxis_title='Time (s)',yaxis_title='Mean floor X (mm)',height=350);st.plotly_chart(fig,use_container_width=True)
    with tabs[1]:
        st.subheader('Frozen input record for this exact run');st.json(r['inputs']);p=engine.section(c)
        st.write(f"Gross A = {p['Ag']:.9g} m²; gross I = {p['Ig']:.9g} m⁴. Each roof node receives {c['roof_kN']/(c['bays']+1):.6f} kN.")
        st.write('Time histories use t=0, dt, 2dt, ... up to the last grid point not exceeding the requested duration.');st.write('Element IDs, selected losses and actual solver stiffnesses:');st.dataframe(pd.DataFrame(a['sections']),hide_index=True,use_container_width=True)
        st.info('Uncracked rectangular RC, symmetric bars, perfect bond, removed steel remains void. No cracking, yield, bond-slip or exposure-time model.')
    with tabs[2]:
        st.subheader('Independent checks computed during this run');st.dataframe(pd.DataFrame(a['checks']),hide_index=True,use_container_width=True)
        st.json(a['equilibrium'])
        st.caption('Two implementations share the stated assumptions but use separate structural solvers. Agreement checks equation implementation; it does not validate real material capacity or omitted physics.')
        st.latex(r'\mathbf F=\mathbf K\mathbf u,\quad\mathbf k_e^g=\mathbf T^T\mathbf k_e^\ell\mathbf T,\quad\mathbf K_T\Delta\mathbf u=\mathbf F_{ext}-\mathbf f_{int}')
        st.write('Dynamics: Newmark β=1/4, γ=1/2. Damping C=2ζω₁M. Translational nodal mass equals floor mass/(bays+1); rotary inertia is an assumed nodal mass × story height²/12.')
        if'nonlinear'in a and'error'not in a['nonlinear']:st.write(f"Nominal max |N|/Aeq = {a['nonlinear']['nominal_N_over_Aeq_MPa']:.3f} MPa. This is not a constituent stress or strength check.")
        st.warning('Reference buckling multipliers and noisy inverse diagnosis are outside this live app. A single mesh or time step does not establish convergence.')
    with tabs[3]:
        u=np.array(a['u']).reshape(-1,3);df=pd.DataFrame(dict(Node=np.arange(1,len(u)+1),Ux_mm=u[:,0]*1000,Uy_mm=u[:,1]*1000,Rz_rad=u[:,2]));st.dataframe(df,hide_index=True,use_container_width=True)
        force=pd.DataFrame(a['local_end_forces'],columns=['Ni_N','Vi_N','Mi_Nm','Nj_N','Vj_N','Mj_Nm']);force.insert(0,'Element',np.arange(1,len(force)+1));st.dataframe(force,hide_index=True,use_container_width=True)
        st.download_button('Download this run: inputs + actual results + checks (JSON)',json.dumps(r,indent=2),file_name=f"OpenSees_run_{r['run_id']}.json",mime='application/json')
        st.download_button('Download nodal results CSV',df.to_csv(index=False),file_name=f"OpenSees_nodes_{r['run_id']}.csv",mime='text/csv')
    st.caption('Local live engineering app: fresh OpenSees solves on Run. Plot zoom, selected view and magnification do not change the analysis inputs. No MATLAB GUI run is implied.')


if __name__=='__main__':
    import sys
    if '--verify' in sys.argv:
        from pathlib import Path
        result=run(DEFAULT)
        Path('OpenSees_actual_results.json').write_text(json.dumps(result,indent=2))
        pd.DataFrame(comparison_rows(result)).to_csv('OpenSees_comparison.csv',index=False)
        export_html(result,'OpenSees_Actual_Results.html')
        print(json.dumps({'OpenSees_version':result['opensees_version'],'run_id':result['run_id'],'all_checks_pass':result['active']['all_checks_pass'] and result['healthy']['all_checks_pass'],'comparison':comparison_rows(result)},indent=2))
    else:
        main()
