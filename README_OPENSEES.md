# FrameLab 2D OpenSees Live

Fresh OpenSees analysis with a browser interface. One Python source file contains the OpenSees model, independent matrix solver, Streamlit GUI and Plotly animation. MATLAB is not required.

## Open the hosted app

**https://framelab2d-opensees.streamlit.app/**

Open the sidebar on a phone, review the preset and choose **RUN LIVE OPENSEES ANALYSIS**. The server performs a new OpenSees solve. No installation is needed. Use `FrameLab2D_OpenSees_QR.png` on your audience slide.

## Run locally

Use 64-bit Python 3.12. Extract this folder, open Terminal inside it, then run:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py
```

Open the Local URL printed in Terminal, usually http://localhost:8501. Keep Terminal running. Press Ctrl+C to stop.

On Windows, replace `source .venv/bin/activate` with `.venv\Scripts\activate`.

OpenSees 3.8 wheels require x86_64 Linux/Windows or arm64 Apple Silicon macOS. Use a matching 64-bit interpreter. This package was executed and tested on Linux x86_64 with Python 3.12, not on your Mac. The exact Mac installation should be checked before the presentation. The program contains no MATLAB Engine dependency and no hard-coded personal folder path.

On a compatible Mac, you can also double-click `START_MAC.command` after extracting the folder. It creates a local virtual environment, installs the listed dependencies and launches the app. First launch needs internet.

## One-minute demonstration

1. Review the visible inputs. The screenshot preset is already loaded.
2. Click RUN LIVE OPENSEES ANALYSIS.
3. Read the fresh run ID, OpenSees version and calculation time.
4. Switch between Active RC (deteriorated) and Healthy RC baseline.
5. Open the dynamic animation and press Play. Inspect static deformation, moment, axial force or P-Delta views.
6. Change roof force from 100 to 150 kN. The old result is hidden until Run is clicked again. A fresh OpenSees solve produces a new run ID.
7. Download the run JSON and nodal CSV. The JSON includes actual inputs, output vectors and verification checks.

Animation Play replays the just-computed response; it does not start another solver. Plot magnification is visual only. A changed input requires Run.

## Exact screenshot inputs

- 2 stories, 2 bays, bay width 6 m, story height 3.5 m.
- Fixed bases: Ux = Uy = Rz = 0; no rigid diaphragm.
- Rectangular section 400 × 400 mm; Ec = 30 GPa, Es = 200 GPa.
- Two 16 mm bars per opposite face; bar-centre cover 50 mm.
- Total roof horizontal force 100 kN, divided equally among the three roof nodes.
- Deterioration: 10% concrete E loss and 10% steel area loss on all members. Healthy baseline has zero loss.
- Floor mass 80 tonnes; first-mode damping 5%; harmonic-force frequency 1.2 Hz; duration 6 s; dt = 0.02 s.
- Gravity 1000 kN/floor × factor 5; four segments/member for corotational Newton analysis; EI factor 1.

## Expected actual OpenSees results

| Metric | Healthy RC | Deteriorated RC |
|---|---:|---:|
| Static roof mean X, mm | 8.646462936 | 9.607217788 |
| First natural frequency, Hz | 1.731469718 | 1.642613260 |
| Dynamic roof peak, mm | 22.699651345 | 28.497863724 |
| P-Delta roof mean X, mm | 9.508357682 | 10.683670747 |

These were obtained by executing OpenSees 3.8.0, not by replaying MATLAB numbers. The independent matrix solver is run again for each live analysis. The MATLAB screenshot reference column appears only when the full screenshot preset matches. A screenshot value is a previously recorded reference, not a fresh MATLAB execution.

## Inspect results without installation

Open `OpenSees_Actual_Results.html` in a desktop browser. It contains actual saved results, input records, zoomable plots and animation controls. It does NOT recalculate new inputs. The Python/Streamlit app does.

## Run without the GUI

```bash
python streamlit_app.py --verify
```

This performs fresh analysis and writes `OpenSees_actual_results.json`, `OpenSees_comparison.csv` and `OpenSees_Actual_Results.html` to the current folder. The saved examples included here were produced by this command.

## Scope and assumptions

This companion verifies static analysis, modal frequency, harmonic dynamic response and corotational P-Delta response. Sensor diagnosis/placement are custom inverse-problem routines, not native OpenSees capabilities, and are not independently reproduced in this companion. Reference elastic buckling factors are not computed here. Do not present this app as a seven-tab OpenSees verification of those routines.

Uncracked, perfectly bonded elastic RC; symmetric bar loss with removed steel treated as void. No yielding, cracking, bond slip or time-dependent corrosion model. Dynamic mass and damping intentionally match FrameLab's stated assumptions. Dynamic loading is an applied roof force, not seismic base acceleration. P-Delta gravity uses dead nodal forces. A numerical match checks implementation under shared assumptions; it is not a capacity/design validation.

OpenSees owns a global solver domain; the source uses a process-wide lock to serialize concurrent solver access. New results are stored per Streamlit session, and input hashes prevent stale results being presented as current.

## Tests included

`TEST_RESULTS.json` records executed solver checks and Streamlit AppTest checks: real Run callback, both states, input-change invalidation, fresh run IDs, six diagram views and response changes after a load edit. Tests ran on Linux, not a native Mac GUI. The hosted Streamlit app was also run with the screenshot preset; both healthy and deteriorated results passed the independent matrix and equilibrium checks.

## Official references

- https://opensees.github.io/OpenSeesDocumentation/user/pythonInstall.html
- https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/algorithm/Newton.html
- https://docs.streamlit.io/develop/concepts/architecture/run-your-app
