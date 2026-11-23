#!/usr/bin/env python3
"""FrameLab 2D independent OpenSeesPy solver, extracted from the deployed source.
Python 3.12. Install numpy scipy openseespy==3.8.0.0.
Run: python FrameLab2D_OpenSees.py --output OpenSees_results.json
The default run compares healthy RC with all-member 10% Ec / 10% steel-area loss.
This file performs actual OpenSees static, modal, dynamic and corotational solves.
The custom matrix checks and reference buckling calculation are separate helpers.
It has no MATLAB dependency and no Streamlit dependency.
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


from scipy.linalg import eigh as symmetric_eigh
def serial(v):
    return json.loads(json.dumps(v,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))

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

def main():
    import argparse
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',help='Optional JSON object overriding DEFAULT inputs.')
    parser.add_argument('--output',default='OpenSees_results.json')
    args=parser.parse_args()
    config=dict(DEFAULT)
    if args.input:
        override=json.loads(Path(args.input).read_text())
        if not isinstance(override,dict):raise ValueError('Input JSON must be an object.')
        unknown=set(override)-set(DEFAULT)-{'active_state','member_losses','ei_factor'}
        if unknown:raise ValueError('Unknown input keys: '+', '.join(sorted(unknown)))
        config.update(override)
    result=run(config)
    Path(args.output).write_text(json.dumps(result,indent=2))
    print('Actual OpenSees version:',result['opensees_version'])
    for key in ['healthy','active']:
        a=result.get(key)
        if a is None:continue
        print(key, 'static roof mm=',a['roof_static_mm'],'f1 Hz=',a['freq_Hz'][0],
              'dynamic peak mm=',a.get('dynamic',{}).get('roof_peak_mm'),
              'P-Delta roof mm=',a.get('nonlinear',{}).get('roof_mm'),
              'checks pass=',a['all_checks_pass'])
    print('Inputs, result vectors and checks written to:',args.output)
    if any(result[k] is not None and not result[k]['all_checks_pass'] for k in ['active','healthy']):
        raise SystemExit('One or more requested verification checks failed. Inspect the JSON.')

if __name__=='__main__':main()
