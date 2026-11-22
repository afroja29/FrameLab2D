const C=require('./solver.js'),assert=require('node:assert/strict');
const near=(a,b,tol=1e-7)=>assert.ok(Math.abs(a-b)<=tol*Math.max(Math.abs(b),1e-20),`${a} != ${b}`);
// Closed-form cantilever tip response and equilibrium (inclined/global transform covered below).
const E=200e9,A=.02,I=8e-5,L=3,F=1000;
const beam={nodes:[[0,0],[L,0]],elements:[{i:0,j:1,E,A,I}],F:[0,0,0,500,F,0],fixed:[0,1,2]};
const r=C.solve(beam);near(r.u[3],500*L/(E*A));near(r.u[4],F*L**3/(3*E*I));near(r.u[5],F*L**2/(2*E*I));near(r.reaction[0],-500);near(r.reaction[1],-F);near(r.reaction[2],-F*L);assert.ok(r.residual<1e-10);
const ang=.71,ca=Math.cos(ang),sa=Math.sin(ang),incl={...beam,nodes:[[0,0],[L*ca,L*sa]],F:[0,0,0,500*ca-F*sa,500*sa+F*ca,0]},ri=C.solve(incl);near(ri.u[3],r.u[3]*ca-r.u[4]*sa);near(ri.u[4],r.u[3]*sa+r.u[4]*ca);
// Section hand reference from the MATLAB package.
const s=C.section({E:30e9,Es:200e9,A:.16,I:.4**4/12,nFace:2,d:.016,a:.05},.1,.2);near(s.EA_ratio,.896741774395939);near(s.EI_ratio,.894594677002191);near(s.As/s.As0,.8);assert.throws(()=>C.section({E:200e9,Es:200e9,A:.02,I:8e-5,nFace:2,d:.016,a:.05}));
const m=C.frame({stories:2,bays:2,width:5,height:3.5,E,A,I,load:100e3});const base=C.solve(m),alpha=C.zeros(m.elements.length);alpha[0]=.3;const damage=C.solve(m,alpha);assert.ok(C.roof(m,damage.u)>C.roof(m,base.u));assert.ok(base.residual<1e-10);
const sensors=m.nodes.map((_,i)=>3*i).filter(d=>!m.fixed.includes(d));const dat=C.simulate(m,sensors,0,.3,0),diag=C.diagnose(m,sensors,dat);assert.equal(diag[0].element,0);near(diag[0].loss,.3);const rec=C.recommend(m,sensors.slice(0,2),diag,2,1e-6);assert.equal(rec.length,2);assert.ok(rec.every(x=>!sensors.slice(0,2).includes(x.dof)));
const sm=C.sensitivity(m,sensors.slice(0,2));assert.ok(sm.rank<=2);assert.ok(sm.unresolved>=8);
// Analytic corotational tangent vs independent central finite differences.
const u=[0,0,.01,.004,.008,-.012],co=C.corot(beam,beam.elements[0],u);let err=0,den=0;for(let j=0;j<6;j++){const up=u.slice(),um=u.slice(),h=1e-6;up[j]+=h;um[j]-=h;const p=C.corot(beam,beam.elements[0],up).f,n=C.corot(beam,beam.elements[0],um).f;for(let i=0;i<6;i++){err+=((p[i]-n[i])/(2*h)-co.K[i][j])**2;den+=co.K[i][j]**2;}}assert.ok(Math.sqrt(err/den)<1e-7);
// Small-load nonlinear solution tends to exact linear solution.
const small={...beam,F:beam.F.map(x=>x*1e-3),bays:1,width:L,stories:1,height:L};const nr=C.newton(small,small.F),lin=C.solve(small);near(nr.u[4],lin.u[4],1e-5);assert.ok(nr.residual<1e-8);
// Buckling: pinned-pinned reference is not used; cantilever fixed-free with 8 segments.
const column={nodes:[[0,0],[0,L]],elements:[{i:0,j:1,E,A,I,kind:'Column',story:1}],F:[0,0,0,1,0,0],fixed:[0,1,2],stories:1,bays:0,width:1,height:L};
const st=C.stability(column,1000,1,1,8),criticalLoad=st.critical*1000;near(criticalLoad,Math.PI**2*E*I/(4*L*L),1e-4);assert.ok(st.second.residual<1e-8);
const dyn=C.dynamics(m,{mass:100000,damping:.02,dt:.02,duration:2,frequency:1});assert.ok(dyn.f1>0&&dyn.peak>0);assert.equal(dyn.time.length,101);
console.log(JSON.stringify({passed:true,cantilever_mm:r.u[4]*1000,EA_ratio:s.EA_ratio,EI_ratio:s.EI_ratio,diagnosed_member:diag[0].element+1,tangent_relative_error:Math.sqrt(err/den),NR_residual:nr.residual,critical_load_N:criticalLoad,critical_exact_N:Math.PI**2*E*I/(4*L*L),dynamic_f1_Hz:dyn.f1},null,2));
