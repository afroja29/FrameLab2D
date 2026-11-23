#!/usr/bin/env python3
"""FrameLab 2D OpenSees Live. Self-contained source, Python 3.12.
Install: python -m pip install openseespy streamlit plotly numpy scipy pandas
Run GUI: python -m streamlit run streamlit_app.py
Run solver only: python FrameLab2D_OpenSees_Live.py --verify
OpenSees is the structural solver; Streamlit/Plotly provide the interface.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MPIR_CVAR_CH3_INTERFACE_HOSTNAME','127.0.0.1')
os.environ.setdefault('MPICH_INTERFACE_HOSTNAME','127.0.0.1')
os.environ.setdefault('HWLOC_COMPONENTS','-linux')

# ---------- Independent matrix solver ----------

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
    if healthy or c.get('active_state')=='Healthy RC':selected=[]
    if not(0<=c['Ec_loss_pct']<=80 and 0<=c['As_loss_pct']<=80):raise ValueError('Loss percentages must be between 0 and 80.')
    el=[];audit=[]
    for e,(i,j)in enumerate(edges):
        ce=c['Ec_loss_pct']/100 if e in selected else 0;ar=c['As_loss_pct']/100 if e in selected else 0
        if 'member_losses' in c and not healthy and c.get('active_state')!='Healthy RC':
            if len(c['member_losses'])!=len(edges):raise ValueError('Member loss table does not match the current frame.')
            entry=c['member_losses'][e];ce=float(entry['Ec_loss_pct'])/100;ar=float(entry['As_loss_pct'])/100
            if not(np.isfinite(ce)and np.isfinite(ar)and 0<=ce<=.8 and 0<=ar<=.8):raise ValueError('Member loss must be finite and between 0 and 80%.')
        p=section(c,ce,ar);p['EI']*=c.get('ei_factor',1.);p['Ieq']*=c.get('ei_factor',1.)
        el.append((i,j,c['Ec_GPa']*1e9,p['Aeq'],p['Ieq']))
        audit.append(dict(Element=e+1,NodeI=i+1,NodeJ=j+1,Type=meta[e][0],Story=meta[e][1],Ec_loss_pct=ce*100,As_loss_pct=ar*100,EA_N=p['EA'],EI_Nm2=p['EI'],Aeq_m2=p['Aeq'],Ieq_m4=p['Ieq'],Steel_area_m2=p['As'],Steel_area_original_m2=p['As0']))
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
    if c['run_nonlinear']:
        try:
            out['buckling']=buckling_reference(m,c)
            if out['buckling'].get('residual') is not None:
                br=out['buckling']['residual'];rows.append(dict(Check='Reference buckling eigen-equation',RelativeError=br,Tolerance=1e-8,Pass=br<1e-8))
        except Exception as exc:
            out['buckling']=dict(error=str(exc));rows.append(dict(Check='Reference buckling completed',RelativeError=None,Tolerance=1e-8,Pass=False))
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
    if c.get('active_state')=='Healthy RC' or c.get('ei_factor',1.)!=1.:return False
    if 'member_losses' in c and any(x['Ec_loss_pct']!=10 or x['As_loss_pct']!=10 for x in c['member_losses']):return False
    return all(c[k]==DEFAULT[k]for k in keys)

def comparison_rows(r):
    c=r['inputs'];preset=is_screenshot_preset(c);records=[]
    for key,label in [('healthy','Healthy RC'),('active',c.get('active_state','Deteriorated RC' if c['scope']=='All members'and c['Ec_loss_pct']==10 and c['As_loss_pct']==10 else 'Active RC'))]:
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
    st.set_page_config(page_title='FrameLab 2D | Live structural analysis',page_icon='🏗️',layout='wide')

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

# ---------- Connected seven-module workflow ----------
import copy
from itertools import combinations
from scipy.linalg import eigh as symmetric_eigh

DIRECTIONS = {'X displacement': 0, 'Y displacement': 1, 'Z rotation': 2}
TAB_NAMES = ['1. Frame model', '2. Sensor network', '3. Damage diagnosis', '4. Dynamics', '5. Sensor placement', '6. P-Delta & buckling', '7. Exposure scenarios']

def serial(v):
    return json.loads(json.dumps(v,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))

def scaled_model(m, element, loss):
    q={**m,'el':list(m['el'])}
    i,j,E,A,I=q['el'][element]
    q['el'][element]=(i,j,E,A*(1-loss),I*(1-loss))
    return q

def static_trials(m, scenarios):
    """One actual OpenSees static solve per member/loss trial, under the shared lock."""
    out=[]
    with LOCK:
        try:
            for element,loss in scenarios:
                q=m if element is None else scaled_model(m,int(element),float(loss))
                out.append(op_static(q)[0])
        finally:ops.wipe()
    return np.asarray(out)

def sensor_dofs(sensors):
    return np.array([3*(int(s['Node'])-1)+DIRECTIONS[s['Direction']] for s in sensors],dtype=int)

def model_id(c):
    return fingerprint({k:c[k] for k in ['stories','bays','width','height','Ec_GPa','Es_GPa','b_mm','h_mm','bars_per_face','diameter_mm','cover_mm','roof_kN','active_state','member_losses','ei_factor']})

def sensitivity_report(m,sensors):
    ds=sensor_dofs(sensors)
    if not len(ds):raise ValueError('Add at least one sensor.')
    uu=static_trials(m,[(None,0)]+[(e,.001)for e in range(len(m['el']))])
    H=((uu[1:,ds]-uu[0,ds])/.001).T
    scales=np.maximum(np.abs(uu[0,ds]),1e-6)
    noise=np.array([s['Noise_percent']for s in sensors])/100
    sigma=np.maximum(noise*scales,1e-9)
    W=H/sigma[:,None];norms=np.linalg.norm(W,axis=0)
    N=W/np.maximum(norms,1e-30)
    sv=np.linalg.svd(N,compute_uv=False);rank=int(np.sum(sv>max(sv[0]if len(sv)else 0,1e-30)*1e-7))
    return serial(dict(H=H,singular_values=sv,rank=rank,unresolved=len(m['el'])-rank,sigma=sigma))

def simulate_measurements(m,sensors,hidden,loss,noise,seed):
    ds=sensor_dofs(sensors)
    if not len(ds):raise ValueError('Add at least one sensor before simulation.')
    uu=static_trials(m,[(None,0),(hidden-1,loss)])
    base,truth=uu[:,ds];percent=np.hypot([s['Noise_percent']for s in sensors],noise)
    sigma=np.maximum(percent/100*np.maximum(np.maximum(abs(truth),abs(base)),1e-6),1e-9)
    draws=np.random.default_rng(int(seed)).standard_normal(len(ds));draws[percent==0]=0
    measured=truth+sigma*draws
    rows=[]
    for j,s in enumerate(sensors):
        rows.append({**s,'DOF':int(ds[j]+1),'Unit':'rad'if ds[j]%3==2 else'm','Baseline':base[j],'Truth':truth[j],'Measured':measured[j],'Sigma':sigma[j],'Effective_noise_percent':percent[j],'Noise_Z':draws[j]})
    return serial(dict(y=measured,sigma=sigma,truth_u=uu[1],rows=rows,hidden=hidden,loss=loss,seed=seed))

def diagnose_model(m,sensors,data):
    ds=sensor_dofs(sensors);y=np.array(data['y']);sigma=np.array(data['sigma']);nm=len(m['el'])
    scenarios=[(e,k/100)for e in range(nm)for k in range(81)]
    predictions=static_trials(m,scenarios)
    residual=predictions[:,ds]-y
    objective=np.sum((residual/sigma)**2,axis=1);rmse=np.sqrt(np.mean(residual**2,axis=1))
    scan=[dict(Element=e+1,Alpha=loss,Objective=float(objective[i]),RMSE_raw=float(rmse[i]))for i,(e,loss)in enumerate(scenarios)]
    best_indices=[e*81+int(np.argmin(objective[e*81:(e+1)*81]))for e in range(nm)]
    best_indices.sort(key=lambda i:(objective[i],scenarios[i][0]))
    best=best_indices[0];ind=ref_linear(scaled_model(m,*scenarios[best]))['u']
    return serial(dict(candidates=[scan[i]for i in best_indices],scan=scan,best_u=predictions[best],predictions=predictions[best_indices],independent_relative_error=rel(ind,predictions[best]),fit_rows=[dict(Node=s['Node'],Direction=s['Direction'],Measured=float(y[j]),Predicted=float(predictions[best,ds[j]]),Sigma=float(sigma[j]),Weighted_residual=float(residual[best,j]/sigma[j]))for j,s in enumerate(sensors)]))

def layout_score(pred,ds,sigma):
    if not len(ds):return 0.
    p,q=np.triu_indices(len(pred),1)
    if not len(p):raise ValueError('At least two scenarios are required.')
    delta=(pred[p][:,ds]-pred[q][:,ds])/sigma
    return float(np.min(np.linalg.norm(delta,axis=1)))

def greedy_layout(pred,existing,pool,count,sigma):
    selected=[];scores=[]
    for _ in range(min(int(count),len(pool))):
        best=-np.inf;chosen=None
        for d in pool:
            if d in selected:continue
            v=layout_score(pred,list(existing)+selected+[d],sigma)
            if v>best+1e-12:best=v;chosen=int(d)
        selected.append(chosen);scores.append(float(best))
    return dict(dofs=selected,nodes=[d//3+1 for d in selected],scores=scores,final_score=layout_score(pred,list(existing)+selected,sigma))

def optimize_sensors(c,sensors,count,sigma,loss,mode,diagnosis=None):
    m=make_model(c);ds=sensor_dofs(sensors);allx=list(range(3*m['ncols'],len(m['F']),3));existing=[int(d)for d in ds if d%3==0]
    pool=[d for d in allx if d not in existing]
    if not pool:raise ValueError('All available horizontal channels already have sensors.')
    if mode=='Diagnosis scenarios':
        if diagnosis is None:raise ValueError('Run diagnosis on the current active model first, or select the controlled RC comparison.')
        rows=diagnosis['candidates'][:5];pred=static_trials(m,[(x['Element']-1,x['Alpha'])for x in rows])
        # This reproduces the MATLAB diagnosis-layout criterion: NEW channels only.
        layout=greedy_layout(pred,[],pool,count,sigma)
        return serial(dict(mode=mode,active=layout,scenarios=rows,predictions=pred,sigma=sigma,existing=existing,scenario_count=len(pred),criterion='New channels only; top five fitted single-member scenarios.'))
    hm=make_model(c,True);scenarios=[(None,0)]+[(e,loss)for e in range(len(m['el']))]
    ph=static_trials(hm,scenarios);pa=static_trials(m,scenarios)
    h=greedy_layout(ph,existing,pool,count,sigma);a=greedy_layout(pa,existing,pool,count,sigma)
    old_on_active=layout_score(pa,existing+h['dofs'],sigma)
    return serial(dict(mode=mode,healthy=h,active=a,existing=existing,scenario_count=len(pa),sigma=sigma,additional_loss=loss,old_on_active=old_on_active,set_changed=set(h['dofs'])!=set(a['dofs']),order_changed=h['dofs']!=a['dofs'],predictions=pa,healthy_predictions=ph,criterion='Existing X channels plus selected channels; baseline and one additional-loss trial per member.'))

def buckling_reference(m,c):
    mesh=refine(m,int(c['segments']));unitc={**c,'gravity_factor':1.};Fg=gravity(mesh,unitc)
    K=ref_linear(mesh)['K'];Kg=np.zeros_like(K)
    if not np.any(Fg):return dict(factors=[],mode=[],residual=None,description='No reference gravity pattern.')
    _,_,actions=op_static(mesh,Fg)
    for e,row in enumerate(mesh['el']):
        k,T,ids,L=ref_elt(mesh,row);N=actions[e,3];g=np.zeros((6,6));js=[1,2,4,5]
        g[np.ix_(js,js)]=N/(30*L)*np.array([[36,3*L,-36,3*L],[3*L,4*L*L,-3*L,-L*L],[-36,-3*L,36,-3*L],[3*L,-L*L,-3*L,4*L*L]])
        Kg[np.ix_(ids,ids)]+=T.T@g@T
    fs=mesh['free'];kf=K[np.ix_(fs,fs)];gf=Kg[np.ix_(fs,fs)]
    vals,V=symmetric_eigh(-gf,kf);keep=np.flatnonzero(vals>max(1e-14,np.max(abs(vals))*1e-10))[::-1][:5]
    factors=1/vals[keep];modes=[];res=[]
    for j,fac in zip(keep,factors):
        u=np.zeros(len(K));v=V[:,j];v=v/max(abs(v));u[fs]=v;modes.append(u)
        res.append(float(np.linalg.norm((kf+fac*gf)@v)/max(1e-30,(np.linalg.norm(kf)+fac*np.linalg.norm(gf))*np.linalg.norm(v))))
    return serial(dict(factors=factors,mode=modes[0]if modes else[],residual=res[0]if res else None,mesh=mesh_dict(mesh),gravity_kN=c['gravity_kN']))

def story_table(m,u1,u2,height):
    result=[];nc=m['ncols']
    for f in range(1,m['stories']+1):
        hi=3*np.arange(f*nc,(f+1)*nc);lo=hi-3*nc;d1=(np.array(u1)[hi]-np.array(u1)[lo])*1000;d2=(np.array(u2)[hi]-np.array(u2)[lo])*1000
        result.append(dict(Story=f,MeanDrift1_mm=float(np.mean(d1)),MeanDrift2_mm=float(np.mean(d2)),MeanDrift1_percent=float(np.mean(d1)/height/10),MeanDrift2_percent=float(np.mean(d2)/height/10),MaxColumnDrift2_percent=float(np.max(abs(d2))/height/10)))
    return result

def summary_metrics(m,a,c):
    u=np.array(a['u']);nc=m['ncols'];dr=[abs(u[3*(f*nc+j)]-u[3*((f-1)*nc+j)])/c['height']*100 for f in range(1,m['stories']+1)for j in range(nc)]
    return {'Roof mean X (mm)':a['roof_static_mm'],'Max column drift (%)':max(dr),'Compliance Fᵀu (N m)':float(m['F']@u)}

def plan_figure(m,sensors,selected=None,healthy=None,affected=None):
    # Fixed front-plane positions keep sensors independent of a scaled deformation plot.
    dummy=dict(mesh=mesh_dict(m),u=np.zeros(len(m['F'])).tolist())
    conf={'stories':m['stories'],'bays':m['ncols']-1,'width':np.ptp(m['xy'][:,0])/(m['ncols']-1),'height':np.ptp(m['xy'][:,1])/m['stories'],'roof_kN':0}
    fig=draw(dummy,'Model and supports',0,True,conf)
    for e,(i,j,*_) in enumerate(m['el']):
        mid=(m['xy'][i]+m['xy'][j])/2;fig.add_annotation(x=mid[0],y=mid[1],text=f'E{e+1}',showarrow=False,xshift=10,yshift=8,font=dict(size=10,color='#65798a'))
    if affected is not None:
        i,j=map(int,m['el'][affected][:2]);fig.add_trace(go.Scatter(x=m['xy'][[i,j],0],y=m['xy'][[i,j],1],mode='lines',line=dict(color='#d8342f',width=6),name='Estimated member'))
    for rows,symbol,color,name in [(sensors,'circle-open','#0075db','Existing sensors'),([dict(Node=d//3+1,Direction='X displacement')for d in (healthy or[])],'diamond-open','#008e65','Healthy recommendation'),([dict(Node=d//3+1,Direction='X displacement')for d in (selected or[])],'star','#ef8b00','Active recommendation')]:
        if rows:
            pts=np.array([m['xy'][int(s['Node'])-1]for s in rows]);text=[f"N{s['Node']}: {s['Direction'][0]}"for s in rows]
            fig.add_trace(go.Scatter(x=pts[:,0],y=pts[:,1],mode='markers+text',marker=dict(symbol=symbol,size=18,line=dict(width=2,color=color),color=color),text=text,textposition='top left'if name=='Existing sensors'else'bottom right',name=name))
    fig.update_layout(title='2D front view | fixed bases | sensor locations',height=540,legend=dict(orientation='h',y=1.18));return fig

def clear_inference():
    for k in ['measurements','diagnosis','placement','adequacy','load_path','mesh_check']:st.session_state.pop(k,None)

def apply_exposure():
    draft=st.session_state.exposure_draft
    try:
        for row in draft:
            for k in ['Ec_loss_pct','As_loss_pct']:
                v=float(row[k])
                if not np.isfinite(v) or not 0<=v<=80:raise ValueError('Enter finite member losses from 0 to 80%.')
    except (TypeError,ValueError) as exc:
        st.session_state.exposure_error='Enter valid numeric Ec and steel-area losses (0–80%) in every row.'
        return
    st.session_state.pop('exposure_error',None)
    st.session_state.member_losses=copy.deepcopy(draft)
    st.session_state.active_state='Deteriorated RC';st.session_state.pop('result',None);clear_inference();st.session_state.pending_run=True

def restore_healthy():
    st.session_state.active_state='Healthy RC';st.session_state.pop('result',None);clear_inference();st.session_state.pending_run=True

def load_demo():
    for k,v in DEFAULT.items():
        if k not in ['scope','custom_members','compare_healthy','Ec_loss_pct','As_loss_pct']:st.session_state[k]=v
    for k in ['geometry_id','member_losses','exposure_draft','sensors','result','editor_revision']:st.session_state.pop(k,None)
    st.session_state.active_state='Healthy RC';st.session_state.ei_factor=1.;clear_inference()

def main():
    st.title('FrameLab 2D')
    st.caption('Structural Analysis, Deterioration Assessment and Sensor Placement · By Afroja Sultana')
    st.caption('Live OpenSees structural responses · FrameLab diagnosis and greedy sensor-placement methods')
    c=dict(DEFAULT)
    with st.sidebar:
        st.header('Shared model inputs')
        st.button('LOAD RC DEMO',on_click=load_demo,key='load_demo',width='stretch')
        with st.expander('Geometry and supports',expanded=True):
            c['stories']=st.number_input('Stories',1,4,2,key='stories');c['bays']=st.number_input('Bays',1,3,2,key='bays')
            c['width']=st.number_input('Bay width (m)',1.,15.,6.,key='width');c['height']=st.number_input('Story height (m)',1.,8.,3.5,key='height')
            st.caption('Fixed bases: Ux = Uy = Rz = 0. Online limit: 4 stories × 3 bays.')
        with st.expander('RC section and reinforcement'):
            c['Ec_GPa']=st.number_input('Concrete Ec (GPa)',1.,100.,30.,key='Ec_GPa');c['Es_GPa']=st.number_input('Steel Es (GPa)',1.,300.,200.,key='Es_GPa')
            c['b_mm']=st.number_input('Gross section width b (mm)',100.,1500.,400.,key='b_mm');c['h_mm']=st.number_input('Gross section depth h (mm)',100.,1500.,400.,key='h_mm')
            c['bars_per_face']=st.number_input('Bars per opposite face',1,12,2,key='bars_per_face');c['diameter_mm']=st.number_input('Bar diameter (mm)',4.,50.,16.,key='diameter_mm');c['cover_mm']=st.number_input('Bar-centre cover (mm)',5.,200.,50.,key='cover_mm')
        with st.expander('Load and dynamics'):
            c['roof_kN']=st.number_input('Total lateral roof force (kN)',-500.,500.,100.,key='roof_kN')
            c['run_dynamic']=st.checkbox('Run dynamic response',True,key='run_dynamic');c['mass_tonne']=st.number_input('Mass per floor (tonnes)',1.,1000.,80.,key='mass_tonne');c['damping_pct']=st.number_input('First-mode damping (%)',0.,30.,5.,key='damping_pct');c['frequency_Hz']=st.number_input('Applied harmonic force frequency (Hz)',.05,10.,1.2,key='frequency_Hz')
            c['dt']=st.number_input('Time step (s)',.005,.1,.02,step=.005,format='%.3f',key='dt');c['duration']=st.number_input('Duration (s)',.1,20.,6.,key='duration')
        with st.expander('P-Delta settings'):
            c['run_nonlinear']=st.checkbox('Run corotational Newton analysis',True,key='run_nonlinear');c['gravity_kN']=st.number_input('Gravity per floor (kN)',0.,10000.,1000.,key='gravity_kN');c['gravity_factor']=st.number_input('Gravity multiplier',0.,50.,5.,key='gravity_factor');c['segments']=st.number_input('Segments per member',1,8,4,key='segments');c['ei_factor']=st.number_input('Effective EI factor (all modules)',.1,1.,1.,key='ei_factor')
        st.radio('Active model',['Healthy RC','Deteriorated RC'],key='active_state')
        st.caption('Exposure losses affect every response module only when Deteriorated RC is active. Diagnosis tests additional loss.')
    geom=(c['stories'],c['bays'])
    if st.session_state.get('geometry_id')!=geom:
        st.session_state.geometry_id=geom
        preview=make_model({**c,'active_state':'Healthy RC'})
        rows=[dict(Element=r['Element'],Type=r['Type'],Story=r['Story'],Exposure='Combined',Ec_loss_pct=10.,As_loss_pct=10.)for r in preview['audit']]
        st.session_state.member_losses=copy.deepcopy(rows);st.session_state.exposure_draft=copy.deepcopy(rows)
        nc=c['bays']+1;nodes=sorted(set([2*nc,(c['stories']+1)*nc]))
        st.session_state.sensors=[dict(Node=n,Direction='X displacement',Noise_percent=1.)for n in nodes]
        st.session_state.editor_revision=st.session_state.get('editor_revision',0)+1
        clear_inference()
    c['active_state']=st.session_state.active_state;c['member_losses']=st.session_state.member_losses;c['compare_healthy']=True
    try:m=make_model(c);healthy_m=make_model(c,True)
    except Exception as exc:st.error(str(exc));return
    key=model_id(c);ds=sensor_dofs(st.session_state.sensors);sk=fingerprint([key,st.session_state.sensors])
    current_hash=fingerprint(c);r=st.session_state.get('result')
    if r is not None and r['input_hash']!=current_hash:r=None
    st.info(f"Active: {c['active_state']} | {c['stories']} stories × {c['bays']} bays | {len(m['el'])} members | roof load {c['roof_kN']:g} kN | {len(ds)} sensors")
    run_pressed=st.button('RUN LIVE OPENSEES ANALYSIS',type='primary',key='run_analysis',width='stretch')
    if run_pressed or st.session_state.pop('pending_run',False):
        try:
            with st.spinner('Solving active and healthy RC in OpenSees; checking independent equations...'):r=run(c)
            st.session_state.result=r
        except Exception as exc:st.error(str(exc));r=None
    if r is None:st.caption('Inputs are ready. Click Run for static, modal, dynamic and stability results. Sensor studies can also be run from their tabs.')
    else:st.caption(f"OpenSees {r['opensees_version']} | Run ID {r['run_id']} | solved in {r['elapsed_s']:.2f} s | all tabs use this shared active RC model")
    tabs=st.tabs(TAB_NAMES)
    with tabs[0]:
        if r:
            a=r['active'];metrics=summary_metrics(m,a,c);cols=st.columns(3)
            for col,(name,value)in zip(cols,metrics.items()):col.metric(name,f'{value:.6f}')
            if a['all_checks_pass']:st.success('OpenSees results pass the independent matrix and equilibrium checks for this run.')
            else:st.error('One or more requested checks failed. Inspect verification below.')
            view=st.selectbox('Diagram',['Model and supports','Static deformation','First vibration mode','Bending moment','Axial force'],index=1,key='diagram')
            scale=st.slider('Deformation magnification (display only)',1,300,100,key='magnification');st.plotly_chart(draw(a,view,scale,True,c),width='stretch',key='static_plot')
            with st.expander('Verification, actual inputs and downloads'):
                st.dataframe(pd.DataFrame(comparison_rows(r)),hide_index=True,width='stretch')
                st.dataframe(pd.DataFrame(a['checks']),hide_index=True,width='stretch');st.json(c)
                st.download_button('Download full run JSON',json.dumps(r,indent=2),file_name=f"FrameLab_run_{r['run_id']}.json",mime='application/json')
                nodal=pd.DataFrame(np.array(a['u']).reshape(-1,3),columns=['Ux_m','Uy_m','Rz_rad']);nodal.insert(0,'Node',np.arange(1,len(nodal)+1))
                st.download_button('Download node CSV',nodal.to_csv(index=False),file_name='FrameLab_nodes.csv')
                st.latex(r'\mathbf F=\mathbf K\mathbf u,\qquad \mathbf K=\sum_e\mathbf A_e^T\mathbf T_e^T\mathbf k_e\mathbf T_e\mathbf A_e')
        else:st.plotly_chart(plan_figure(m,[]),width='stretch',key='preview')
        st.caption('Uncracked elastic RC, perfect bond, fixed bases, no diaphragm. Online limits keep live calculations manageable; no material capacity check is implied.')
    with tabs[1]:
        st.subheader('Existing sensor network')
        ca,cb,cc=st.columns(3)
        node=ca.number_input('Sensor node',m['ncols']+1,len(m['xy']),len(m['xy']),key='sensor_node')
        direction=cb.selectbox('Measurement direction',list(DIRECTIONS),key='sensor_direction');snoise=cc.number_input('Channel noise (%)',0.,20.,1.,key='sensor_noise')
        if st.button('ADD SENSOR',key='add_sensor'):
            s=dict(Node=int(node),Direction=direction,Noise_percent=float(snoise))
            if any(x['Node']==s['Node']and x['Direction']==s['Direction']for x in st.session_state.sensors):st.warning('This channel already exists. Remove it before changing its noise.')
            else:st.session_state.sensors.append(s);clear_inference();st.rerun()
        sensor_rows=[{**s,'Channel':i+1,'DOF':int(sensor_dofs([s])[0]+1)}for i,s in enumerate(st.session_state.sensors)]
        st.dataframe(pd.DataFrame(sensor_rows),hide_index=True,width='stretch')
        if sensor_rows:
            remove=st.selectbox('Channel to remove',list(range(1,len(sensor_rows)+1)),key='remove_channel')
            if st.button('REMOVE SELECTED',key='remove_sensor'):st.session_state.sensors.pop(remove-1);clear_inference();st.rerun()
        st.plotly_chart(plan_figure(m,st.session_state.sensors),width='stretch',key='network_plot')
        if st.button('CHECK SENSOR LIMITATIONS',key='sensor_adequacy'):
            try:
                with st.spinner('Computing OpenSees finite-difference sensitivity...'):rep=sensitivity_report(m,st.session_state.sensors)
                st.session_state.adequacy=dict(key=sk,data=rep)
            except Exception as exc:st.error(str(exc))
        adeq=st.session_state.get('adequacy')
        if adeq and adeq['key']==sk:
            rep=adeq['data'];st.write(f"{len(ds)} channels; sensitivity rank {rep['rank']}/{len(m['el'])}; {rep['unresolved']} unresolved parameter combinations.")
            if rep['unresolved']:st.warning('Simultaneous loss in all members cannot be independently identified. A single-member hypothesis can still be tested; ranking is not confidence.')
            fig=go.Figure(go.Bar(x=list(range(1,len(rep['singular_values'])+1)),y=rep['singular_values']));fig.update_layout(xaxis_title='Singular-value index',yaxis_title='Column-normalized weighted sensitivity',height=300);st.plotly_chart(fig,key='sensitivity_plot',width='stretch')
            st.download_button('Download sensitivity CSV',pd.DataFrame(rep['H']).to_csv(index=False),file_name='FrameLab_sensitivity.csv')
    with tabs[2]:
        st.subheader('Additional single-member stiffness loss')
        ca,cb=st.columns(2);hidden=ca.number_input('Hidden member',1,len(m['el']),1,key='hidden_member');loss=cb.slider('Additional stiffness loss (%)',0,80,30,key='damage_loss')/100
        ca,cb=st.columns(2);noise=ca.number_input('Extra simulation noise (%)',0.,20.,0.,key='simulation_noise');seed=cb.number_input('Random seed',0,999999,7405,key='noise_seed')
        st.caption('Additional loss multiplies both EA and EI of one member. Channel and extra noise combine in quadrature. Set BOTH noise inputs to zero for an exact recovery test. Python and MATLAB random draws need not match.')
        measurement_key=fingerprint([sk,int(hidden),loss,noise,int(seed)])
        if st.button('SIMULATE DATA',key='simulate_data'):
            try:
                data=simulate_measurements(m,st.session_state.sensors,int(hidden),loss,noise,int(seed));st.session_state.measurements=dict(key=measurement_key,data=data);st.session_state.pop('diagnosis',None);st.session_state.pop('placement',None)
            except Exception as exc:st.error(str(exc))
        data_state=st.session_state.get('measurements');data=data_state['data']if data_state and data_state['key']==measurement_key else None
        if data:
            st.dataframe(pd.DataFrame(data['rows']),hide_index=True,width='stretch')
            st.download_button('Download actual measurements CSV',pd.DataFrame(data['rows']).to_csv(index=False),file_name='FrameLab_measurements.csv')
        else:st.info('Simulate measurements after selecting the model, channels and noise settings.')
        if st.button('RUN DIAGNOSIS',key='run_diagnosis',disabled=data is None):
            try:
                with st.spinner('OpenSees is solving every member at 0–80% additional loss, in 1% steps...'):diag=diagnose_model(m,st.session_state.sensors,data)
                st.session_state.diagnosis=dict(key=measurement_key,model_sensor_key=sk,data=diag);st.session_state.pop('placement',None)
            except Exception as exc:st.error(str(exc))
        diag_state=st.session_state.get('diagnosis');diag=diag_state['data']if diag_state and diag_state['key']==measurement_key else None
        if diag:
            best=diag['candidates'][0];st.success(f"Best fit: E{best['Element']}, {best['Alpha']*100:.0f}% additional loss | J = {best['Objective']:.6g} | raw RMSE = {best['RMSE_raw']:.6g}")
            st.plotly_chart(plan_figure(m,st.session_state.sensors,affected=best['Element']-1),width='stretch',key='diagnosis_plan')
            st.dataframe(pd.DataFrame(diag['candidates']),hide_index=True,width='stretch')
            chart=go.Figure(go.Bar(x=[x['Element']for x in diag['candidates']],y=[x['Alpha']*100 for x in diag['candidates']]));chart.update_layout(xaxis_title='Candidate element',yaxis_title='Best-fit additional loss (%)',height=300);st.plotly_chart(chart,width='stretch',key='diagnosis_chart')
            st.dataframe(pd.DataFrame(diag['fit_rows']),hide_index=True,width='stretch')
            gap=diag['candidates'][1]['Objective']-best['Objective']if len(diag['candidates'])>1 else 0
            st.caption('Raw RMSE uses the channel units; with mixed translations and rotations use the dimensionless weighted objective for ranking.')
            st.caption(f"Next-candidate objective gap: {gap:.6g}. Ranking is not probability or confidence. Independent best-fit displacement difference: {diag['independent_relative_error']:.3g}.")
            st.latex(r'J(e,\alpha)=\sum_i[(y_i-u_i(e,\alpha))/\sigma_i]^2,\qquad (\hat e,\hat\alpha)=\arg\min J')
            st.download_button('Download complete trial scan CSV',pd.DataFrame(diag['scan']).to_csv(index=False),file_name='FrameLab_diagnosis_scan.csv')
    with tabs[3]:
        st.subheader('Linear harmonic dynamic response')
        if r and 'dynamic'in r['active']:
            a=r['active'];cols=st.columns(2);cols[0].metric('First natural frequency',f"{a['freq_Hz'][0]:.6f} Hz");cols[1].metric('Peak mean roof response',f"{a['dynamic']['roof_peak_mm']:.6f} mm")
            st.plotly_chart(dynamic_animation(a,c),width='stretch',key='dynamic_animation')
            fig=go.Figure()
            for j,y in enumerate(a['dynamic']['floor_mm']):fig.add_trace(go.Scatter(x=a['dynamic']['time_s'],y=y,name=f'Floor {j+1}'))
            fig.update_layout(xaxis_title='Time (s)',yaxis_title='Mean floor X (mm)',height=400);st.plotly_chart(fig,width='stretch',key='dynamic_history')
            dd=pd.DataFrame({'Time_s':a['dynamic']['time_s'],**{f'Floor_{j+1}_mm':v for j,v in enumerate(a['dynamic']['floor_mm'])}});st.download_button('Download dynamic histories CSV',dd.to_csv(index=False),file_name='FrameLab_dynamics.csv')
            st.latex(r'\mathbf M\ddot{\mathbf u}+\mathbf C\dot{\mathbf u}+\mathbf K\mathbf u=\mathbf F\sin(2\pi f t),\quad \mathbf C=2\zeta\omega_1\mathbf M')
            st.caption('Play replays this freshly solved response. Applied harmonic roof force, not base acceleration. Newmark β=¼, γ=½; nodal rotary inertia = nodal mass × story height²/12.')
        else:st.info('Enable dynamic response in the sidebar, then run the shared analysis.')
    with tabs[4]:
        st.subheader('Greedy sensor placement')
        mode=st.selectbox('Scenario mode',['Diagnosis scenarios','Healthy vs deteriorated RC'],key='optimization_mode')
        if mode=='Healthy vs deteriorated RC' and c['active_state']=='Healthy RC':st.info('The active model is healthy. These layouts will be identical until an exposure scenario is applied in tab 7.')
        ca,cb,cc=st.columns(3);count=ca.number_input('New X sensors',1,12,2,key='new_sensors');sigma=cb.number_input('Comparison sigma X (mm)',.001,10.,.1,format='%.3f',key='comparison_sigma')/1000;design_loss=cc.number_input('Design additional loss (%)',1.,80.,30.,key='design_loss')/100
        opt_key=fingerprint([sk,mode,count,sigma,design_loss,measurement_key if mode=='Diagnosis scenarios'else None])
        ds_now=st.session_state.get('diagnosis');opt_diag=ds_now['data']if ds_now and ds_now['key']==measurement_key else None
        opt_clicked=st.button('OPTIMIZE SENSOR LOCATIONS',key='optimize_sensors')
        reopt=st.session_state.pop('request_reoptimization',False)
        if reopt:
            mode='Healthy vs deteriorated RC';opt_key=fingerprint([sk,mode,count,sigma,design_loss,None])
        if opt_clicked or reopt:
            try:
                with st.spinner('OpenSees trial responses and greedy maximin scenario separation...'):placement=optimize_sensors(c,st.session_state.sensors,count,sigma,design_loss,mode,opt_diag)
                st.session_state.placement=dict(key=opt_key,model_sensor_key=sk,data=placement)
            except Exception as exc:st.error(str(exc))
        p_state=st.session_state.get('placement');p=p_state['data']if p_state and p_state['key']==opt_key else None
        if p:
            a=p['active'];h=p.get('healthy');st.plotly_chart(plan_figure(m,st.session_state.sensors,a['dofs'],h['dofs']if h else[]),width='stretch',key='placement_plot')
            rows=[]
            for i,(d,score)in enumerate(zip(a['dofs'],a['scores'])):
                row=dict(Pick=i+1,Active_node=d//3+1,DOF=d+1,Direction='X',Separation_score=score)
                if h:row.update(Healthy_node=h['nodes'][i],Pick_changed=h['dofs'][i]!=d,Healthy_score=h['scores'][i])
                rows.append(row)
            st.dataframe(pd.DataFrame(rows),hide_index=True,width='stretch')
            st.write(f"Recommended nodes: {a['nodes']} | minimum cumulative separation {a['final_score']:.6g} sigma | {p['scenario_count']} scenarios")
            if h:
                metrics=[{'Metric':'Existing X channels only','Healthy':layout_score(np.array(p['healthy_predictions']),p['existing'],sigma),'Active':layout_score(np.array(p['predictions']),p['existing'],sigma)},{'Metric':'Re-optimized layout score','Healthy':h['final_score'],'Active':a['final_score']},{'Metric':'Healthy layout evaluated on each model','Healthy':h['final_score'],'Active':p['old_on_active']}]
                st.dataframe(pd.DataFrame(metrics),hide_index=True,width='stretch')
                st.info(f"Layout set changed: {p['set_changed']}. Pick order changed: {p['order_changed']}. Old layout on active RC: {p['old_on_active']:.6g}; new layout: {a['final_score']:.6g}. Uniform losses can change scores without changing locations.")
                st.caption('Green diamonds: healthy recommendation. Orange stars: active recommendation. Blue circles: existing sensors. Same scenarios, sensor budget and absolute noise in both models.')
            st.caption(p['criterion']+' Scores below 1 sigma mean weak minimum separation. Greedy selection is not a guaranteed global optimum.')
            st.latex(r'j^*=\arg\max_j\min_{p<q}\sqrt{\sum_{i\in S\cup\{j\}}[(u_i^{(p)}-u_i^{(q)})/\sigma]^2}')
            st.download_button('Download placement audit JSON',json.dumps(p,indent=2),file_name='FrameLab_sensor_placement.json')
            def accept_layout(dofs=tuple(a['dofs'])):
                for d in dofs:st.session_state.sensors.append(dict(Node=d//3+1,Direction='X displacement',Noise_percent=1.))
                clear_inference()
            st.button('ADD RECOMMENDED SENSORS TO NETWORK',key='accept_layout',on_click=accept_layout)
        else:st.info('Run diagnosis first for fitted scenarios, or choose Healthy vs deteriorated RC for a controlled design comparison.')
    with tabs[5]:
        st.subheader('Geometric nonlinear equilibrium and reference buckling')
        if r and 'nonlinear'in r['active']:
            a=r['active'];nl=a['nonlinear'];bk=a.get('buckling',{})
            if 'error'in nl:st.error(nl['error'])
            else:
                first=nl['first_order_combined_roof_mm'];second=nl['roof_mm'];cols=st.columns(3);cols[0].metric('First-order mean roof',f'{first:.6f} mm');cols[1].metric('Newton mean roof',f'{second:.6f} mm');cols[2].metric('Amplification',f'{second/first:.6f}'if abs(first)>1e-12 else'Undefined')
                st.plotly_chart(draw(a,'P-Delta deformation',100,True,c),width='stretch',key='pdelta_plot')
                mesh=refine(m,int(c['segments']));u1=ref_linear(mesh,F=mesh['F']+gravity(mesh,c))['u']
                st.dataframe(pd.DataFrame(story_table(m,u1,nl['u'],c['height'])),hide_index=True,width='stretch')
                st.download_button('Download stability results JSON',json.dumps({'inputs':c,'nonlinear':nl,'buckling':bk},indent=2),file_name='FrameLab_stability.json')
            if bk.get('factors'):
                critical=bk['factors'][0];st.write(f"Reference elastic critical gravity factor: {critical:.6f}; applied/reference = {c['gravity_factor']/critical:.6f}; eigen-equation residual = {bk['residual']:.3g}.")
                bview=st.checkbox('Show first reference buckling mode',key='show_buckling')
                if bview:
                    dummy={'mesh':bk['mesh'],'u':bk['mode'],'mode1':bk['mode']};st.plotly_chart(draw(dummy,'First vibration mode',1,True,c).update_layout(title='Reference elastic buckling mode | arbitrary amplitude'),key='buckling_mode',width='stretch')
                st.dataframe(pd.DataFrame({'Mode':range(1,len(bk['factors'])+1),'Reference_factor':bk['factors']}),hide_index=True,width='stretch')
                if c['gravity_factor']>=critical:st.warning('Applied gravity reaches/exceeds the elastic reference. Do not interpret a converged algebraic solution as a stable capacity result.')
            if st.button('CALCULATE GRAVITY–SWAY PATH',key='load_path_button'):
                path=[]
                with st.spinner('Solving independent OpenSees states along the gravity path...'):
                    mesh=refine(m,int(c['segments']));upper=min(max(1.,c['gravity_factor']*1.1),.85*bk['factors'][0])if bk.get('factors')else max(1.,c['gravity_factor']*1.1)
                    with LOCK:
                        try:
                            for factor in np.unique(np.r_[np.linspace(0,upper,13),c['gravity_factor']]):
                                try:
                                    uu,_,_=op_static(mesh,mesh['F']+gravity(mesh,{**c,'gravity_factor':float(factor)}),True);path.append(dict(Gravity_factor=float(factor),Roof_X_mm=float(np.mean(uu[3*m['roof']])*1000)))
                                except Exception:break
                        finally:ops.wipe()
                st.session_state.load_path=dict(key=current_hash,rows=path)
            path=st.session_state.get('load_path')
            if path and path['key']==current_hash:
                fig=go.Figure(go.Scatter(x=[v['Roof_X_mm']for v in path['rows']],y=[v['Gravity_factor']for v in path['rows']],mode='lines+markers'));fig.update_layout(xaxis_title='Mean roof X (mm)',yaxis_title='Gravity multiplier',height=350);st.plotly_chart(fig,key='gravity_path',width='stretch')
            with st.expander('Equations and mesh convergence'):
                st.latex(r'\mathbf K_T\Delta\mathbf u=\mathbf F_{ext}-\mathbf f_{int},\qquad(\mathbf K+\lambda_{cr}\mathbf K_g)\boldsymbol\phi=0')
                st.caption('Reference buckling uses elastic stiffness and geometric stiffness from OpenSees gravity-only axial forces at multiplier 1. It is a custom matrix eigenvalue calculation, not a native OpenSees buckling command or collapse prediction.')
                if st.button('CHECK CURRENT MESH AGAINST FINER MESH',key='mesh_check_button'):
                    try:
                        finer=min(8,int(c['segments'])*2)
                        if finer==c['segments']:raise ValueError('Already at the maximum online subdivision. Reduce the current mesh to compare.')
                        cm={**c,'segments':finer};fm=refine(m,finer)
                        if len(fm['el'])>160:raise ValueError('The finer mesh exceeds the 160-element online limit. Use a smaller frame.')
                        with LOCK:
                            try:uf,_,_=op_static(fm,fm['F']+gravity(fm,cm),True)
                            finally:ops.wipe()
                        val=float(np.mean(uf[3*m['roof']])*1000);st.session_state.mesh_check=dict(key=current_hash,finer=finer,roof=val,difference_pct=100*abs(val-nl['roof_mm'])/max(abs(val),1e-12))
                    except Exception as exc:st.error(str(exc))
                ch=st.session_state.get('mesh_check')
                if ch and ch['key']==current_hash:st.write(f"{c['segments']} → {ch['finer']} segments/member: finer roof {ch['roof']:.6f} mm, relative difference {ch['difference_pct']:.5g}%.")
        else:st.info('Enable corotational Newton analysis in the sidebar and run the shared analysis.')
        st.caption('Elastic materials only; geometry is nonlinear. A small residual proves numerical equilibrium, not real-world safety. No yield, cracking, crushing or bond failure is simulated.')
    with tabs[6]:
        st.subheader('Assumed exposure-related deterioration')
        if st.session_state.get('exposure_error'):st.error(st.session_state.exposure_error)
        st.caption('Exposure labels are descriptions. Enter assumed Ec loss and steel AREA loss; steel Es stays constant. The model does not predict losses from chloride/sulfate exposure or elapsed time.')
        ca,cb,cc=st.columns(3);group=ca.selectbox('Assign to group',['All members','Columns','Beams','Bottom columns','Left columns','Selected story'],key='exposure_group');story=cb.number_input('Exposure story',1,c['stories'],1,key='exposure_story');label=cc.selectbox('Exposure label',['Normal','Sulfate','Chloride / seawater','Combined'],index=3,key='exposure_label')
        ca,cb=st.columns(2);ce=ca.number_input('Assumed Ec loss (%)',0.,80.,10.,key='assign_ec');ar=cb.number_input('Assumed steel area loss (%)',0.,80.,10.,key='assign_as')
        if st.button('ASSIGN TO GROUP',key='assign_exposure'):
            draft=copy.deepcopy(st.session_state.exposure_draft)
            for row in draft:
                i=row['Element']-1;col=int(m['el'][i][0])%m['ncols']
                use=group=='All members'or(group=='Columns'and row['Type']=='Column')or(group=='Beams'and row['Type']=='Beam')or(group=='Bottom columns'and row['Type']=='Column'and row['Story']==1)or(group=='Left columns'and row['Type']=='Column'and col==0)or(group=='Selected story'and row['Story']==story)
                if use:row.update(Exposure=label,Ec_loss_pct=ce,As_loss_pct=ar)
            st.session_state.exposure_draft=draft;st.session_state.editor_revision+=1
        edited=st.data_editor(pd.DataFrame(st.session_state.exposure_draft),hide_index=True,width='stretch',disabled=['Element','Type','Story'],column_config={'Exposure':st.column_config.SelectboxColumn(options=['Normal','Sulfate','Chloride / seawater','Combined']),'Ec_loss_pct':st.column_config.NumberColumn(min_value=0,max_value=80),'As_loss_pct':st.column_config.NumberColumn(min_value=0,max_value=80)},key=f"exposure_editor_{st.session_state.editor_revision}")
        st.session_state.exposure_draft=edited.to_dict('records')
        ca,cb,cc=st.columns(3);ca.button('APPLY AND COMPARE',on_click=apply_exposure,key='apply_exposure',width='stretch');cb.button('RESTORE HEALTHY RC',on_click=restore_healthy,key='restore_healthy',width='stretch')
        def queue_optimization():
            st.session_state.optimization_mode='Healthy vs deteriorated RC';st.session_state.request_reoptimization=True
        cc.button('RE-OPTIMIZE SENSORS',on_click=queue_optimization,key='exposure_reopt',width='stretch')
        st.caption('Apply commits the table and activates deteriorated RC. The shared analysis updates all response tabs automatically. Restore returns every module to healthy composite stiffness. Re-optimize calculates a fresh controlled layout in tab 5 without changing the active state.')
        if fingerprint(st.session_state.exposure_draft)!=fingerprint(st.session_state.member_losses):st.warning('Exposure table has unapplied edits. Click APPLY AND COMPARE to activate them.')
        st.dataframe(pd.DataFrame(m['audit']),hide_index=True,width='stretch')
        if r:
            hh=summary_metrics(healthy_m,r['healthy'],c);aa=summary_metrics(m,r['active'],c)
            if 'nonlinear'in r['active']and'error'not in r['active']['nonlinear']and'nonlinear'in r['healthy']and'error'not in r['healthy']['nonlinear']:
                hh['P-Delta roof X (mm)']=r['healthy']['nonlinear']['roof_mm'];aa['P-Delta roof X (mm)']=r['active']['nonlinear']['roof_mm']
            table=[{'Metric':name,'HealthyRC':value,'ActiveRC':aa[name],'Change_percent':100*(aa[name]-value)/abs(value)if abs(value)>1e-12 else None}for name,value in hh.items()]
            st.dataframe(pd.DataFrame(table),hide_index=True,width='stretch')
            st.download_button('Download exposure comparison CSV',pd.DataFrame(table).to_csv(index=False),file_name='FrameLab_exposure_comparison.csv')
        with st.expander('Section calculation details'):
            e=st.number_input('Calculation member',1,len(m['el']),1,key='section_member')-1;row=m['audit'][e];p0=section(c);p=section(c,row['Ec_loss_pct']/100,row['As_loss_pct']/100)
            st.latex(r'EA=E_{c0}(1-\ell_c)A_c+E_s A_s,\quad EI=E_{c0}(1-\ell_c)I_c+E_s I_s')
            st.latex(r'A_s=(1-r)A_{s0},\quad I_s=2n[(1-r)^2 I_b+(1-r)A_bz^2],\quad z=h/2-a')
            st.write(f"EA healthy = {p0['EA']:.8g} N; active = {p['EA']:.8g} N; ratio = {p['EA']/p0['EA']:.8f}. EI healthy = {p0['EI']:.8g} N m²; active = {p['EI']:.8g} N m²; ratio = {p['EI']/p0['EI']:.8f}.")
            k,T,_,L=ref_elt(m,m['el'][e]);st.write('Local element stiffness k (mixed SI force/displacement units):');st.dataframe(pd.DataFrame(k),width='stretch')
            st.caption(f"Shared effective-EI factor = {c.get('ei_factor',1.):g}; it multiplies the composite EI shown here in every module.")
            st.caption('Gross concrete area excludes the original bar holes. Removed steel remains void. Bar-centre cover is surface-to-bar-centre distance, not clear cover.')
            st.download_button('Download section and model audit JSON',json.dumps({'inputs':c,'sections':m['audit'],'healthy_section':p0,'selected_section':p},indent=2),file_name='FrameLab_section_audit.json')
    with st.expander('Method, limits and official references'):
        st.write('Static, dynamics and every diagnosis/design scenario are fresh OpenSees solves. Diagnosis grid search and greedy maximin sensor placement are FrameLab algorithms. Buckling is a separate reference matrix eigenproblem. Independent matrix checks compare solver implementations under the same assumptions.')
        st.write('Online: up to 4 stories, 3 bays, 160 stability elements and 1500 time steps. Sensors are displacement/rotation channels; recommendations add horizontal-displacement sensors. No real sensor hardware connection, plasticity, cracking or time-dependent corrosion. Results may reset when the browser session closes or the hosting service sleeps.')
        st.markdown('[OpenSees elastic beam-column](https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html) · [Corotational transformation](https://openseespydoc.readthedocs.io/en/latest/src/corotational.html) · [Source and MATLAB download](https://github.com/afroja29/FrameLab2D)')


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
