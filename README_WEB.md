# FrameLab 2D — Browser edition

**Structural Analysis, Deterioration Assessment and Sensor Placement**

By Afroja Sultana

## Open and run

Open the public homepage in a modern browser. No MATLAB, Python, account or installation is needed. The app starts with a solved two-story, two-bay example. Choose **Load demo** to reset it. Calculations run locally on your device; inputs are not uploaded. Refreshing the page resets the session.

For an offline backup, download the web package, extract it and open `index.html`. Keep `solver.js`, `app.js` and `styles.css` beside it. An Internet connection is used only for optional fonts and external links.

## One-minute audience example

1. Open the homepage and choose **Load demo**. Read the mean roof displacement (about **33.2056 mm**) and inspect the fixed supports and orange deformed shape.
2. Open **Damage diagnosis**. Keep member **1**, stiffness loss **30%** and noise **1%**.
3. Choose **Simulate data**, then **Run diagnosis**. Compare the inferred member and loss with the synthetic truth.
4. Open **Sensor placement** and choose **Optimize placement** to see suggested new X-displacement channels.

A seeded random generator makes web demonstrations reproducible. The web random sequence is not the MATLAB random sequence. With limited sensors and noise, the estimated member/loss may differ from the truth. Candidate bars are separate single-member hypotheses, not simultaneous damage.

## Modules

| Module | Calculation and outputs |
| --- | --- |
| Frame model | 2D Euler–Bernoulli direct stiffness analysis; displacements, end forces, reactions, story drifts and CSV export. |
| Sensor network | Add or remove X, Y or rotation channels at non-base nodes. Column-normalized finite-difference sensitivity spectrum and local rank. |
| Damage diagnosis | Synthetic measurements with per-channel noise; weighted least-squares scan of one damaged member at a time, with loss from 0% to 80% in 1% increments. |
| Dynamic simulation | Linear elastic harmonic roof loading; lumped translational and rotational mass; mass-proportional damping; Newmark average acceleration; frequency and all-floor histories. |
| Sensor placement | Greedy max–min separation among the top five fitted damage hypotheses; candidate channels measure horizontal displacement. |
| P-Delta & buckling | Corotational elastic full Newton–Raphson, analytic current tangent, adaptive load stepping, backtracking, equilibrium log and reference elastic eigenbuckling. |
| Deterioration | Member/group Ec and reinforcement area losses; healthy-composite versus deteriorated RC; section ratios, deformed response, local stiffness matrix, JSON audit. |

## Inputs, units and supports

Frame geometry uses meters. Elastic modulus uses GPa in the UI and Pa internally. Area is in m², second moment is in m⁴, roof load and gravity use kN in the UI and N internally. Displacements are displayed in mm, rotations in radians and drift in percent. Three DOFs are retained per node: [ux, uy, θz]. All bases are fixed, restraining all three DOFs. Roof lateral load is shared equally across roof nodes.

The browser interface supports 1–4 stories and 1–3 bays. Initial E, A and I are uniform across members. Deterioration inputs vary by member. The initial steel-like demo is not an RC section. Choose **RC 400 mm preset** to set Ec = 30 GPa, A = 0.16 m², I = 0.4⁴/12 m⁴ throughout the model. This resets loss assumptions and downstream results.

## Direct stiffness equations

For local end DOFs [ui, vi, θi, uj, vj, θj]ᵀ, local force is **fe = ke ue**. The axial coefficient is EA/L; bending coefficients are 12EI/L³, 6EI/L², 4EI/L and 2EI/L. Each element is rotated using **kg = Tᵀ ke T**, assembled globally and solved on the free DOFs using **Kff uf = Ff**. Fixed reactions are **R = K u − F**. Damage α scales an entire member's elastic stiffness by 1 − α.

## RC hand check

Use RC 400 mm preset, two 16 mm bars per face, Es = 200 GPa and 50 mm bar-centre cover. Cover is measured from the concrete face to the bar centre, not to the bar surface. The section is rectangular, uncracked, perfectly bonded and reinforced symmetrically on two faces.

Let q be the concrete modulus loss fraction, r the steel area loss fraction, n the bars per face, Ab = πd²/4, Ib = πd⁴/64 and z = h/2 − cover.

- As₀ = 2n Ab; Is₀ = 2n (Ib + Ab z²).
- Ac = Ag − As₀; Ic = Ig − Is₀. Concrete area/inertia stay fixed as steel is removed.
- As = (1 − r) As₀.
- Is = 2n [(1 − r)² Ib + (1 − r) Ab z²].
- EA₀ = Ec Ac + Es As₀; EI₀ = Ec Ic + Es Is₀.
- EA = (1 − q) Ec Ac + Es As; EI = (1 − q) Ec Ic + Es Is.
- Solver equivalents: Aeq = EA/Ec and Ieq = EI/Ec.

At **10% Ec loss and 20% As loss**, As₀ = 804.2477 mm², As = 643.3982 mm², **EA/EA₀ = 0.8967417744** and **EI/EI₀ = 0.8945946770**. Each member can be assigned separately; comparison uses healthy composite RC under the same loads and supports. Zero losses produce identical healthy/deteriorated responses.

Exposure labels describe assumptions. They do not calculate corrosion rates, exposure duration, chloride ingress or loss percentages. Steel loss remains void; Es is constant. Cracking, spalling, bond loss, plasticity, crushing and ultimate capacity are not represented.

## Newton–Raphson and buckling

The geometric nonlinear solver rebuilds internal force and the consistent tangent every iteration: **Kt Δu = Fext − fint(u)**. Basic deformations are chord extension and end rotations relative to the current chord. Dead global nodal loads are used. The weighted residual tolerance is 1e−8; moments are normalized by frame length. Nonconvergence or a nonpositive tangent produces an error rather than an unconverged answer.

Gravity factor λ scales gravity only, while the roof lateral load stays constant. EI factor changes bending stiffness. Reference buckling solves Kφ = −λ Kgφ using gravity-prestress geometric stiffness. The critical factor is relative to the specified per-floor reference gravity. With zero reference gravity it is undefined. It is not a nonlinear collapse load. Use 1, 2 and 4 segments per member to check mesh convergence; do not infer postbuckling behavior.

## Dynamics

Newmark β = 1/4 and γ = 1/2 integrate the linear elastic equations. Floor mass is shared equally among nodes for each translational direction; nodal rotational inertia is nodeMass × storyHeight²/12. Damping C = 2ζω1 M is calibrated to the first natural frequency. The graph is a harmonic force response, not an earthquake record. Time-step convergence must be checked for the chosen frequency and structural modes. This web release plots histories; it does not reproduce the MATLAB dynamic animation.

## Sensor interpretation

Adequacy is a local sensitivity rank check, not a guarantee of unique global damage identification. The weighted diagnosis objective accounts for each channel's noise scale. If rotation and displacement channels are mixed, raw RMSE combines rad and m, so rank candidates using the dimensionless weighted objective. Recommended sensors are scenario-dependent and are not a global optimum.

The browser noise setting is a single effective percentage for all channels. Unlike the MATLAB app, it does not combine an additional per-channel noise percentage in quadrature. Synthetic noise has a small numerical floor even at 0%.

## Verification

Run the included automated numerical checks from a local source checkout with Node.js:

```sh
node solver.test.cjs
```

Checks cover closed-form axial/flexural cantilever response, inclined coordinate transformation, support reactions and free equilibrium, RC section hand values, exact synthetic damage recovery with full X channels, sensor rank and recommendations, analytic corotational tangent versus central differences, the nonlinear small-load limit, refined cantilever Euler buckling and dynamic output sanity. The cantilever's 8-segment elastic critical load differs from π²EI/(4L²) by about 0.00021%. These are educational verification checks, not certification of general structural design accuracy. The browser JavaScript is independently translated from the MATLAB formulation; complete MATLAB/web output parity is not claimed.

## Files and audience sharing

- `index.html`: GUI homepage.
- `solver.js`: structural calculations without runtime dependencies.
- `app.js`: controls, graphics, tables and exports.
- `styles.css`: responsive layout.
- `readme.html` and `README_WEB.md`: browser guide and editable instructions.
- `WEB_VERIFICATION.json`: recorded numerical check results.
- The original MATLAB app remains available in the GitHub repository.

Present the live homepage link and its web QR code at the start of the demo. Test on another computer and a phone at least 24 hours before presenting. Keep an offline web folder and a short screen recording as backups.
