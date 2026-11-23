# FrameLab 2D: connected seven-tab Streamlit app

**Live app: https://framelab2d-opensees.streamlit.app/**

Use the existing `FrameLab2D_OpenSees_QR.png` on the presentation slide. The URL and QR have not changed. The browser app needs no MATLAB, Python installation or user account. On a phone, expand the sidebar for shared inputs and scroll the tab strip horizontally for the seven modules. Calculations execute on the server.

## One-minute audience example

1. Click **LOAD RC DEMO**, then **RUN LIVE OPENSEES ANALYSIS**. Healthy RC is the starting state.
2. Open **3. Damage diagnosis**, click **SIMULATE DATA**, then **RUN DIAGNOSIS**. The simulation has 30% additional loss in member 1. Existing channels have 1% noise even when extra simulation noise is zero.
3. Open **5. Sensor placement**, choose **Diagnosis scenarios**, and click **OPTIMIZE SENSOR LOCATIONS**. Inspect stars, existing sensors and separation scores.
4. Open **7. Exposure scenarios**. The draft has assumed 10% Ec loss and 10% steel-area loss for each member. Click **APPLY AND COMPARE**. This activates deteriorated RC and solves the shared analysis automatically.
5. Click **RE-OPTIMIZE SENSORS**, then open tab 5. Compare the healthy and active layouts under the same additional-loss scenarios, sensor budget and fixed absolute noise. Uniform losses may change scores without changing node choices.
6. **RESTORE HEALTHY RC** returns every response and inference module to healthy composite stiffness and reruns the shared response analysis.

For an exact diagnosis check, remove/re-add the two default channels (nodes 6 and 9, X displacement) with **Channel noise = 0%**, and set **Extra simulation noise = 0%**. The 1% search grid contains the simulated 30% loss, so member 1 and 30% are recovered with zero residual. Python's seeded Gaussian draws are reproducible here but do not reproduce MATLAB's random sequence.

## Seven connected modules

| Tab | Inputs and calculation | Results |
|---|---|---|
| 1. Frame model | Shared geometry, RC section, roof force, active state | OpenSees static displacement, drift, compliance, diagrams, independent equation checks, input record and exports |
| 2. Sensor network | Add/remove a non-base node's X, Y or rotation channel; channel noise | Front-plane sensor map, finite-difference sensitivity spectrum, rank and unresolved parameter combinations |
| 3. Damage diagnosis | Hidden member and additional loss for simulation, noise, seed | Actual measurements and noise audit; OpenSees scan of all members at 0–80% loss in 1% increments; ranked fits, objective, residuals and full scan export |
| 4. Dynamics | Floor mass, first-mode damping, harmonic roof force frequency, time step and duration | OpenSees modal frequency, Newmark histories, peak response and response animation |
| 5. Sensor placement | Scenario mode, new-X-sensor budget, fixed sigma, trial loss | Greedy maximin scenario separation, recommended nodes, healthy/active layouts, scores and complete response audit |
| 6. P-Delta & buckling | Gravity pattern, multiplier, mesh segments | Corotational OpenSees Newton response checked against independent Newton implementation; story drift, amplification, gravity-sway path, custom elastic buckling eigenproblem, optional mesh comparison |
| 7. Exposure scenarios | Group assignment or individual member table; assumed Ec/steel-area losses | Active section EA/EI and remaining steel area, healthy/active response differences, section matrix and exports, restore and fresh sensor re-optimization |

Shared **Active model** explicitly identifies Healthy RC or Deteriorated RC. Pending exposure-table edits do not change the model until Apply. Model/sensor/noise changes hide stale inference results. Exposure application does not compound previous losses: every section is rebuilt from the original healthy RC section. Sensor recommendations can be added to the actual network using the button in tab 5.

## Numerical methods and independence

Every structural response used by diagnosis and sensor-placement scenarios is obtained by building and solving an OpenSees `elasticBeamColumn` frame. FrameLab performs the grid search and greedy selection around those responses. These inverse/design routines are custom algorithms, not native OpenSees features. This is an independent implementation of the MATLAB workflow, not a fresh MATLAB execution.

Diagnosis minimizes `J = sum(((measured - predicted)/sigma)^2)` over element and trial additional loss. Both EA and EI of that trial member are reduced by `1-alpha`. Channel noise and extra simulation noise combine in quadrature. With zero actual noise, sigma is floored at 1e-9 solely to make weighting well defined. Ranking and separation scores are not probabilities or confidence levels.

For **Diagnosis scenarios**, the placement criterion reproduces the MATLAB method: minimum pairwise separation over the **new channels only**, using the five best fitted candidate scenarios and user-entered fixed sigma. Existing X channels are excluded from the candidate pool. For the **controlled RC comparison**, the criterion includes existing X channels plus new channels; scenarios are the unchanged baseline and one fixed additional-loss trial per member. Identical scenarios, candidate pool, budget and sigma are used in both RC models. Greedy placement is not guaranteed to be globally optimal.

The buckling reference solves `(K + lambda Kg) phi = 0`; K is the elastic matrix and Kg is assembled using OpenSees gravity-only axial forces at multiplier 1. This is a **custom linearized matrix eigenproblem**, not a built-in OpenSees eigenbuckling command and not a nonlinear collapse capacity. Actual second-order displacement uses the OpenSees corotational transformation and Newton algorithm. The shared EI factor applies to every module, not just the stability tab.

## Executed preset results

2 stories × 2 bays, 6 m width, 3.5 m height; fixed bases; 400 × 400 mm uncracked RC; Ec = 30 GPa, Es = 200 GPa, two 16 mm bars per opposite face, 50 mm bar-centre cover; 100 kN total roof force split equally across roof nodes.

Dynamics: 80 tonnes/floor, 5% first-mode damping, 1.2 Hz harmonic force, dt = 0.02 s, duration = 6 s. Translational mass is equally divided among floor nodes in X and Y; assumed rotary inertia is nodal mass × story height²/12. No rigid-floor constraint. Damping is C = 2 zeta omega1 M. Newmark beta = 1/4, gamma = 1/2. P-Delta uses 1000 kN/floor × multiplier 5, four segments/member, shared EI factor 1.

| Quantity | Healthy RC | All-member 10% Ec + 10% steel-area loss |
|---|---:|---:|
| Static mean roof X, mm | 8.646462936 | 9.607217788 |
| First frequency, Hz | 1.731469718 | 1.642613260 |
| Dynamic roof peak, mm | 22.699651345 | 28.497863724 |
| P-Delta mean roof X, mm | 9.508357682 | 10.683670747 |

Healthy reference elastic critical factor: **50.831683745**. The cantilever Euler buckling reference check with four segments differs from the analytical value by about **0.00328%**. Increasing this preset's deteriorated stability mesh from 4 to 8 segments/member changes roof response by approximately **0.1141%**. These checks do not establish convergence for every user configuration.

`STREAMLIT_NUMERICAL_TESTS.json` and `STREAMLIT_GUI_TESTS.json` record actual automated checks. They cover exact noiseless recovery, rank 2 for two channels/ten members, diagnosis-layout nodes 7 and 4, zero-loss identity, unaffected members, noncompounding application, Euler reference, GUI callbacks, exposure invalidation/recalculation, re-optimization, healthy restore, load edits, mesh comparison and sensor add/remove. Automated tests ran on Linux/Python 3.12; physical phone/Mac testing is not claimed. Hosted behavior is checked separately after deployment.

## Run locally

Use compatible 64-bit **Python 3.12**. Open a terminal in the extracted folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py
```

On Windows use `.venv\Scripts\activate` instead of `source ...`. On compatible Macs, `START_MAC.command` launches the same app. First launch installs dependencies. If macOS blocks the downloaded launcher, use the terminal commands above. The hosted app avoids local installation entirely.

Streamlit Community Cloud: repository `afroja29/FrameLab2D`, branch `main`, main file `streamlit_app.py`, Python **3.12**. No API secrets are required. OpenSees global solver access is protected by a process-wide lock; session results and input fingerprints prevent other visitors or edited inputs from reusing stale current results.

## Limits and assumptions

Online guardrails: **1–4 stories, 1–3 bays, at most 160 stability elements and 1500 time steps**. These are interactive resource limits. This is a planar RC moment-frame tool, not a general 3D modeler or code-design application. No yielding, concrete cracking/crushing, bond slip, exposure-time progression, real sensor hardware connection, seismic base acceleration or collapse-capacity prediction. Exposure percentages are user assumptions; labels do not calculate corrosion.

The section model has symmetric opposite-face reinforcement, perfect bond and uniform bar-area reduction at fixed bar-centre positions. Removed steel remains void. All members share the same original rectangle/reinforcement, with member-specific loss entries. Recommendation candidates are non-base X-displacement channels; Y/rotation channels can be used for diagnosis but are not selected by this placement method. Mobile layout supports normal browser use, but heavy calculations and hosting sleep/restarts can delay/reset a session.

## Official references

- https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html
- https://openseespydoc.readthedocs.io/en/latest/src/corotational.html
- https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/algorithm/Newton.html
- https://docs.streamlit.io/develop/api-reference/caching-and-state/st.session_state
