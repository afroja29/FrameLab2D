# Standalone OpenSeesPy structural solver

This is the actual structural engine extracted from the deployed FrameLab 2D Streamlit source. It builds its own OpenSees frame and calculates results independently. It does not read MATLAB screenshots or GUI result values. No MATLAB or Streamlit installation is required for this standalone file.

## Install and run

Use compatible 64-bit Python 3.12. In this extracted folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements_opensees.txt
python FrameLab2D_OpenSees.py --output OpenSees_results.json
```

On Windows replace the activation command with `.venv\Scripts\activate`.

The default run solves BOTH healthy RC and all-member 10% Ec / 10% steel-area loss. It writes the actual input record, OpenSees version, run ID, displacement/reaction/force vectors, frequencies, complete dynamic histories and verification checks into JSON. Example results included in this ZIP were produced by executing this standalone file, not by copying MATLAB results.

To change inputs, edit `example_inputs.json`, then run:

```bash
python FrameLab2D_OpenSees.py --input example_inputs.json --output my_results.json
```

The code supports 1–4 stories, 1–3 bays, 160 total stability elements and 1500 dynamic time steps. Use a supported OpenSeesPy wheel for your operating system and architecture. This package was executed on Linux/Python 3.12; native Mac execution is not claimed.

## Read the code

- `make_model`: geometry, fixed bases, composite RC properties and member losses.
- `build`: `ops.model`, `ops.node`, `ops.fix`, `ops.geomTransf` and `ops.element('elasticBeamColumn', ...)`.
- `op_static`: real OpenSees static solution with Linear or corotational/Newton analysis.
- `eigen`: OpenSees vibration eigenproblem.
- `op_dynamic`: OpenSees mass, damping, harmonic roof load and Newmark time integration.
- `ref_linear`, `ref_newton`, `ref_dynamics`: independent matrix checks.
- `buckling_reference`: custom elastic eigenbuckling calculation using OpenSees gravity-only axial forces. This is not a native OpenSees buckling command.
- `run`: runs the active/healthy cases and returns the frozen input/output record.

The input `E` is Ec, not steel Es. The composite section supplies equivalent EA and EI with a constant Ec reference. Removed reinforcement remains void. P-Delta uses elastic members with corotational geometry and Newton equilibrium. Dynamics uses a prescribed roof force, not ground acceleration. Reference buckling does not establish nonlinear collapse capacity.

## Expected default results

| Quantity | Healthy RC | Deteriorated RC |
|---|---:|---:|
| Mean static roof X, mm | 8.646462936 | 9.607217788 |
| First vibration frequency, Hz | 1.731469718 | 1.642613260 |
| Dynamic roof peak, mm | 22.699651345 | 28.497863724 |
| P-Delta mean roof X, mm | 9.508357682 | 10.683670747 |

Defaults: 2 stories, 2 bays, 6 m spans, 3.5 m stories, 400 mm square section, Ec 30 GPa, Es 200 GPa, two 16 mm bars per opposite face, 50 mm bar-centre cover, 100 kN total roof force. Mass 80 tonnes/floor, damping 5%, frequency 1.2 Hz, duration 6 s, time step 0.02 s. Gravity 1000 kN/floor × 5, four segments/member.

This standalone engine covers structural analysis and its checks. For sensor-network controls, synthetic data, diagnosis, greedy placement and exposure re-optimization, use the complete `streamlit_app.py` source or `FrameLab2D_Streamlit_Full.zip` from the repository. Those routines are FrameLab algorithms around OpenSees response trials.

- Full web app: https://framelab2d-opensees.streamlit.app/
- Full source: https://github.com/afroja29/FrameLab2D
- Official element: https://openseespydoc.readthedocs.io/en/latest/src/elasticBeamColumn.html
- Corotational transformation: https://openseespydoc.readthedocs.io/en/latest/src/corotational.html
