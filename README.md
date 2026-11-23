# FrameLab 2D

**Structural Analysis, Deterioration Assessment and Sensor Placement**  
By Afroja Sultana

## Start here

**[Download hub and web access](https://afroja29.github.io/FrameLab2D/)**

| Version | Open / download | Requirements |
|---|---|---|
| MATLAB GUI | [Complete MATLAB ZIP](https://github.com/afroja29/FrameLab2D/raw/refs/heads/main/FrameLab2D_Full_Package.zip) · [Main .m file](FrameLab2D.m) · [Run guide](README_MATLAB.md) | Desktop MATLAB |
| Full live web app, seven tabs | **[Run FrameLab 2D online](https://framelab2d-opensees.streamlit.app/)** · [Source and launcher ZIP](FrameLab2D_Streamlit_Full.zip) · [Web guide](README_OPENSEES.md) | Browser, no installation |
| Independent OpenSees solver code | [FrameLab2D_OpenSees.py](FrameLab2D_OpenSees.py) · [Solver ZIP](FrameLab2D_OpenSees_Code.zip) · [Instructions](README_OPENSEES_CODE.md) | Python 3.12 + OpenSeesPy, NumPy, SciPy |
| Existing JavaScript browser edition | [Open browser edition](https://afroja29.github.io/FrameLab2D/web/) · [Guide](web/README_WEB.md) | Browser, calculations execute locally in JavaScript |

The MATLAB source and package were synchronized to the latest connected version on **November 23, 2026**. The MATLAB demo now uses gravity multiplier **5** to match the verification example. MATLAB and OpenSees use separate structural solvers. The JavaScript browser edition uses its own local solver.

## MATLAB

1. Download and extract **FrameLab2D_Full_Package.zip**.
2. In MATLAB, select the extracted **FrameLab2D** folder as Current Folder.
3. Enter:

```matlab
FrameLab2D
```

Or open and run `START_FRAMELAB.m`. The main class and filename must stay **FrameLab2D**. The main `.m` file contains the calculation helpers. No Python or OpenSees installation is required for this MATLAB GUI.

Click **LOAD RC DEMO**, then **RUN STATIC ANALYSIS**. The healthy RC mean roof X is approximately **8.6465 mm**. The demo has 2 stories, 2 bays, 6 m spans, 3.5 m story height and fixed bases. It uses a 400 × 400 mm uncracked RC rectangle with Ec 30 GPa, Es 200 GPa, two 16 mm bars per opposite face and 50 mm bar-centre cover. A 100 kN total lateral roof force is distributed equally among roof nodes. Gravity is 1000 kN/floor × multiplier 5, with four segments/member and EI factor 1.

Optional MATLAB checks:

```matlab
FrameLab2D.verifyConnectedWorkflow
FrameLab2D.verifyNewtonRaphson
FrameLab2D.verifyEnvironment
verifyFrameLab2D
```

`verifyFrameLab2D` opens a temporary GUI, checks current callbacks and closes it. The preparation environment did not have MATLAB, so native MATLAB execution is not claimed. Syntax checks and independent numerical tests are documented separately. See [README_MATLAB.md](README_MATLAB.md).

## Live web app

Open **https://framelab2d-opensees.streamlit.app/**. On a phone, expand the sidebar for shared inputs and scroll the tab strip to see every module.

1. Choose **LOAD RC DEMO**, then **RUN LIVE OPENSEES ANALYSIS**.
2. In tab 3, choose **SIMULATE DATA**, then **RUN DIAGNOSIS**.
3. In tab 5, choose **OPTIMIZE SENSOR LOCATIONS**.
4. In tab 7, assign assumed losses and choose **APPLY AND COMPARE**. The active state changes to deteriorated RC and the response analysis runs again.
5. Choose **RE-OPTIMIZE SENSORS**, then inspect tab 5. It compares healthy and active layouts under the same trial scenarios, absolute noise and sensor budget.
6. **RESTORE HEALTHY RC** restores the healthy composite model for all modules.

OpenSees builds and solves a fresh structural model for each Run. Every diagnosis/design trial uses OpenSees structural responses. FrameLab's custom code performs the diagnosis grid search, greedy sensor placement and reference buckling eigenproblem. Those are not native OpenSees features. Animation Play replays the newly calculated time history. The MATLAB screenshot column is a recorded reference and not a fresh MATLAB solve.

## Seven connected modules

1. **Frame model:** static response, drift, compliance, deformation and force diagrams.
2. **Sensor network:** add/remove measurement channels and inspect sensitivity rank/limitations.
3. **Damage diagnosis:** simulate readings, audit noise, scan additional single-member EA/EI loss and rank candidates.
4. **Dynamics:** frequencies, harmonic-force histories and animation.
5. **Sensor placement:** greedy scenario separation, healthy/active layout comparison and recommendation scores.
6. **P-Delta and buckling:** geometric nonlinear equilibrium, reference elastic eigenvalue and optional mesh comparison.
7. **Exposure scenarios:** member/group losses, section EA/EI calculations, healthy/deteriorated comparison, restore and sensor re-optimization.

All response modules use the active stiffness model. Exposure percentages are assumptions, not losses predicted from a label. Green diamonds mark healthy recommendations, orange stars mark active recommendations, and blue circles mark existing channels. Uniform deterioration may change scores without changing node choices. Diagnosis estimates an additional loss relative to the active baseline. Noise can affect the estimated loss; set both channel and extra simulation noise to zero for an exact 30% grid-search recovery test.

## Audience links and QR codes

| Use | Link | QR image |
|---|---|---|
| Download MATLAB code | [MATLAB download hub](https://afroja29.github.io/FrameLab2D/#matlab) | [FrameLab2D_MATLAB_Download_QR.png](FrameLab2D_MATLAB_Download_QR.png) |
| Run the full live app | [Live web app](https://framelab2d-opensees.streamlit.app/) | [FrameLab2D_Live_Web_QR.png](FrameLab2D_Live_Web_QR.png) |

The MATLAB QR opens the download page; MATLAB desktop is needed to execute `.m` files. The live-web QR works through a browser without MATLAB/Python installation on the viewer's device. The earlier `FrameLab2D_OpenSees_QR.png` continues to open the same live app. Existing JavaScript web QR images still point to that separate edition.

## Verification and limits

For the matching 2-story/2-bay example:

| Result | Healthy RC | All-member 10% Ec + 10% steel-area loss |
|---|---:|---:|
| Mean static roof X, mm | 8.646462936 | 9.607217788 |
| First frequency, Hz | 1.731469718 | 1.642613260 |
| Dynamic roof peak, mm | 22.699651345 | 28.497863724 |
| P-Delta mean roof X, mm | 9.508357682 | 10.683670747 |

OpenSees version 3.8.0 calculated these values independently. Full standalone inputs, vectors and checks are in the solver ZIP. The web app's independent matrix checks run again with user inputs. Sensor-design and rank values are custom post-processing of response trials. A numerical match checks implementations under shared assumptions and is not a structural capacity validation.

The live web app has interactive limits of 4 stories, 3 bays, 160 stability elements and 1500 time steps. The model is planar and elastic. There is no material plasticity, concrete cracking/crushing, bond-slip, time-dependent corrosion or real sensor-hardware connection. Reference elastic buckling is not nonlinear collapse capacity. See each version's guide for its exact scope.
