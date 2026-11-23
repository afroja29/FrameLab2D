# FrameLab 2D: MATLAB desktop package

By Afroja Sultana. Public connected release, updated November 23, 2026.

## Run in under two minutes

1. Download and extract `FrameLab2D_Full_Package.zip`.
2. In MATLAB, select the extracted `FrameLab2D` folder as Current Folder.
3. Type `FrameLab2D` in the Command Window, or open and run `START_FRAMELAB.m`.
4. Click **LOAD RC DEMO**, then **RUN STATIC ANALYSIS**.

`FrameLab2D.m` is the self-contained main file. Its class and filename must remain `FrameLab2D`. Do not rename it to include a space. You need desktop MATLAB with AppBase/uifigure/uigridlayout support. Python, OpenSees and Streamlit are not required to run the MATLAB GUI. This is source code, not a standalone MATLAB Runtime executable.

## One-minute example

Healthy RC composite is the baseline. The demo has two stories, two 6 m bays, 3.5 m story height, Ec = 30 GPa, a 400 × 400 mm gross rectangle and 100 kN total roof force. Reinforcement: Es = 200 GPa, two 16 mm bars per opposite face, 50 mm bar-centre cover. Every column base is fixed.

- Tab 1: Run static analysis. Mean healthy roof X is approximately **8.6465 mm**.
- Tab 2: Existing X channels are nodes 6 and 9. Add/remove sensors and use **SENSOR LIMITATIONS** for rank.
- Tab 3: Simulate member 1 with 30% additional stiffness loss, then run diagnosis. Channel noise and extra simulation noise combine. Set both to zero for an exact grid-search recovery check.
- Tab 4: Run harmonic dynamics. Healthy and deteriorated states give different results.
- Tab 5: Optimize the diagnosis scenarios, or select **Healthy vs deteriorated RC** and compare layouts under identical scenarios, noise and sensor budget.
- Tab 6: The demo gravity is **1000 kN/floor × multiplier 5**, effective EI factor 1, four segments/member. Run P-Delta. Healthy reference roof sway is approximately **9.5084 mm**. The reference elastic critical factor is approximately **50.8317**.
- Tab 7: Assign assumed losses, then **APPLY AND COMPARE**. For all-member 10% Ec loss and 10% steel-area loss, static roof X is approximately **9.6072 mm**. Use **RE-OPTIMIZE SENSORS** for a fresh healthy/active comparison. **RESTORE HEALTHY RC** returns the model to healthy composite stiffness.

Loss percentages are assumptions, not predictions from an exposure label. Every response module uses the same active model. Diagnosis tests an additional single-member loss. Green diamonds mark healthy recommendations, orange stars mark active recommendations. Uniform deterioration can change scores without changing node choices.

## Optional checks in your MATLAB installation

```matlab
FrameLab2D.verifyConnectedWorkflow
FrameLab2D.verifyNewtonRaphson
FrameLab2D.verifyEnvironment
verifyFrameLab2D
```

The last command opens and closes a temporary GUI. `verifyFrameLab2DExport(folder)` recalculates diagnosis residuals/objectives from exact exported observation files. The scripts now use the current button names, including SENSOR LIMITATIONS, CALCULATION DETAILS and RESTORE HEALTHY RC.

This release was syntax-checked. The preparation environment did not have MATLAB, so execution of these MATLAB checks and the native GUI must be confirmed in MATLAB. Independent Python/OpenSees runs check the mathematical implementation under shared assumptions; they are not proof that MATLAB's GUI ran here.

## Public access

- Download hub: https://afroja29.github.io/FrameLab2D/#matlab
- Repository and individual `.m` files: https://github.com/afroja29/FrameLab2D
- Full live web app: https://framelab2d-opensees.streamlit.app/
- Existing JavaScript browser edition: https://afroja29.github.io/FrameLab2D/web/

The MATLAB download QR opens the download hub. MATLAB desktop is needed to execute the `.m` file. The web QR opens the installation-free server app.

## Scope

Planar elastic moment frames with fixed bases. Static diagnosis and dynamics are linear. P-Delta uses a corotational element and Newton-Raphson. Reference elastic buckling is not nonlinear collapse capacity. The RC model assumes uncracked concrete, perfect bond and symmetric bar loss with removed steel treated as void. There is no real sensor hardware connection, plasticity, cracking or time-dependent corrosion model. No structural safety or code-compliance claim is made.
