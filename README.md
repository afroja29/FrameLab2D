# FrameLab 2D

**Structural Analysis, Deterioration Assessment and Sensor Placement**

By Afroja Sultana

> **Browser version:** [Open FrameLab 2D](https://afroja29.github.io/FrameLab2D/web/) · [Read the web guide](https://afroja29.github.io/FrameLab2D/web/readme.html)
>
> No MATLAB, installation or account is needed for the browser edition. Choose **Load demo** and start exploring. All structural calculations run locally in your browser. The web release covers all seven analysis areas with the scope and limitations documented in [README_WEB.md](web/README_WEB.md).

## Connected seven-tab Streamlit app with OpenSees

[Open the live OpenSees app](https://framelab2d-opensees.streamlit.app/) · [Run guide](README_OPENSEES.md) · [Audience QR code](FrameLab2D_OpenSees_QR.png)

No installation is needed for this hosted app. On a phone, open the sidebar, review the preset and click **RUN LIVE OPENSEES ANALYSIS**. It builds and solves an independent OpenSees model on the server for each Run. All seven modules are connected: frame analysis, sensor network, damage diagnosis, dynamics, greedy sensor placement, P-Delta/reference buckling and exposure comparison with re-optimization. Structural trial responses use OpenSees; diagnosis, placement and reference buckling are FrameLab algorithms around those responses. The initial active state is Healthy RC. Apply a member-loss table in tab 7 to activate deteriorated RC.

Download the complete Streamlit source and launcher: [FrameLab2D_Streamlit_Full.zip](FrameLab2D_Streamlit_Full.zip). The existing OpenSees QR still opens this updated app.

## Try it in one minute

1. Open the live browser app and choose **Load demo**.
2. Open **Damage diagnosis**, choose **Simulate data**, then **Run diagnosis**.
3. Open **Sensor placement** and choose **Optimize placement**.
4. For the RC example, choose **Deterioration**, then **RC 400 mm preset**, **Assign to group**, and **Apply and compare**.

Use `web/FrameLab2D_Web_QR.png` for direct audience access. `FrameLab2D_QR.png` links to the source repository. The MATLAB files and original MATLAB ZIP remain available below.

## Start the MATLAB app

1. Extract the ZIP. You can keep the extracted `FrameLab2D` folder in OneDrive.
2. Open MATLAB and select that folder as the Current Folder.
3. Type this in the MATLAB Command Window:

```matlab
FrameLab2D
```

The main file and class are both named `FrameLab2D`. MATLAB identifiers cannot contain spaces. The displayed app title is **FrameLab 2D**. This is a self-contained MATLAB UI app; the main file contains its calculation helpers. Use desktop MATLAB with support for AppBase, uifigure, uigridlayout and uitable. No external mesh or input file is needed for the default example. This package does not contain a MATLAB installer or a standalone executable.

## Workflow and features

- **2D Frame Model:** define the planar frame, sections and loads; calculate the linear elastic response. Nodes have horizontal translation, vertical translation and in-plane rotation. Fixed base supports restrain all three DOFs and are shown with base plates and hatch lines.
- **Sensor Network:** view and edit sensors in the front plane. The sensor adequacy check reports sensitivity rank and unresolved parameter combinations. Limited rank can still permit a useful single-member diagnosis.
- **Damage Diagnosis:** simulate measurements and compare candidate single-member stiffness losses. The fitted frame and the trial bar chart have separate plot tabs. Candidate bar values are separate hypotheses, not simultaneous damage in every member.
- **Dynamic Simulation:** inspect the linear elastic dynamic response.
- **Sensor Recommendation:** compare candidate sensor locations and show suggested sensors in front view.
- **P-Delta & Buckling:** solve geometric nonlinear elastic equilibrium using a corotational element and full Newton-Raphson with adaptive load steps and backtracking. Inspect convergence and tangent stiffness. The reported eigenbuckling factor is a reference linearized result, not a nonlinear collapse load.
- **Exposure scenarios:** assign assumed concrete modulus and reinforcement area losses by member or group, compare healthy RC with deteriorated RC, and inspect/export the section calculations and hand-check matrices.

## Exposure demonstration

1. Select **RC 400 mm PRESET** before assigning losses. This sets a valid 400 x 400 mm RC example; the general frame defaults are not that RC section.
2. Choose members or a group. Enter, for example, concrete E loss = 10% and steel area loss = 20%, then choose **ASSIGN TO GROUP**.
3. Choose **APPLY AND COMPARE**. Read **Section calculations** for EA/EI ratios and **Response and hand check** for response comparisons.
4. Use **HAND CHECK** for formulas and element matrices. Use **COMPARE P-DELTA** for the geometric nonlinear comparison. **RESTORE ORIGINAL** returns to the original frame properties.

Exposure labels describe the scenario. They do not automatically predict corrosion or determine loss percentages. With zero loss inputs, the healthy and deteriorated RC models have the same stiffness. Loss percentages are user assumptions, not field measurements or exposure-to-time predictions.

The section model assumes rectangular, uncracked concrete, perfect bond and symmetric top/bottom reinforcement. Removed steel remains a void; the steel elastic modulus is constant. Bar-centre cover is the distance from the concrete surface to the bar centre, not clear cover. Valid reinforcement geometry must fit inside the section.

## Hand calculation and verification

`HandCalculation_Reference.csv` contains an independent reference for a 400 mm square section with Ec = 30 GPa, Es = 200 GPa, two 16 mm bars per face and 50 mm bar-centre cover. With 10% concrete modulus loss and 20% steel area loss, the EA ratio is about 0.896742 and the EI ratio about 0.894595. Its cantilever tip values are a separate hand-check example, not the default frame response.

In MATLAB, run the numerical checks without opening the UI:

```matlab
FrameLab2D.verifyNewtonRaphson
FrameLab2D.verifyEnvironment
```

For the full app smoke check, which opens and closes a temporary app:

```matlab
verifyFrameLab2D
```

The supplied MATLAB source was syntax-checked, but the MATLAB runtime and GUI were not available in the preparation environment. The included MATLAB checks must be run in MATLAB to confirm local execution. The Newton-Raphson formulation was independently implemented in Python and compared with OpenSees 3.8.0 on seven cases; the maximum relative displacement difference was approximately 1.82e-9. These comparisons verify the independent formulation, not execution of the MATLAB file. Details are in `OpenSees_NR_Comparison.json` and `Adaptive_NR_Check.txt`.

Static diagnosis and dynamics remain linear elastic. Newton-Raphson is used for geometric nonlinear equilibrium in the stability module. Material plasticity, concrete cracking/crushing, bond deterioration and corrosion progression are not simulated.

## Files

- `FrameLab2D.m`: main app and solver helpers.
- `verifyFrameLab2D.m`: MATLAB numerical and app smoke checks.
- `HandCalculation_Reference.csv`: section and cantilever hand-check reference.
- `OpenSees_NR_Comparison.json` and `Adaptive_NR_Check.txt`: recorded independent numerical comparisons.


