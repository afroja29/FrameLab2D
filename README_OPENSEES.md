# FrameLab 2D OpenSees Live

Public browser interface for an independent 2D OpenSeesPy frame model. Change the visible inputs and click **RUN LIVE OPENSEES ANALYSIS** to build and solve a fresh model. No MATLAB installation is needed by visitors.

## Deploy on Streamlit Community Cloud

Upload the contents of this folder to the intended GitHub repository. In Streamlit Community Cloud, create an app from that repository, select its actual default branch, set the main file to `streamlit_app.py`, and choose Python **3.12** in Advanced settings. The requirements file and theme configuration are included. Set the app's viewer access to public, then test the actual live URL before generating a QR code.

No API key or application secret is needed. Authentication to the hosting service and permission to access the source repository must be completed by the account owner. A public URL is not reserved or confirmed by this package.

## Default demonstration

2-story, 2-bay frame; 6 m bays; 3.5 m stories; fixed column bases; RC section 400 × 400 mm; Ec 30 GPa; Es 200 GPa; two 16 mm bars per face; 50 mm bar-centre cover. A total 100 kN roof force is equally distributed. The active case applies 10% concrete E loss and 10% steel area loss to all members, and is compared against matching healthy RC.

Dynamics: 80 tonnes/floor, first-mode damping 5%, applied force frequency 1.2 Hz, 6 s duration, 0.02 s step. P-Delta: 1000 kN/floor times gravity factor 5; four segments/member.

| Result | Healthy RC | Deteriorated RC |
|---|---:|---:|
| Static mean roof X, mm | 8.646462936 | 9.607217788 |
| First frequency, Hz | 1.731469718 | 1.642613260 |
| Dynamic roof peak, mm | 22.699651345 | 28.497863724 |
| P-Delta mean roof X, mm | 9.508357682 | 10.683670747 |

These reference values were obtained by executing OpenSees 3.8.0 on Linux, not by displaying MATLAB outputs. New runs include an independent matrix check and a unique run ID. The MATLAB screenshot column is a recorded comparison reference and appears only for the matching preset. Animation replays a newly solved response; input edits require another Run.

## Mobile use

Open the hosted URL or scan its QR. Expand the sidebar using the menu at the top left to edit inputs, then click Run. Results and interactive plots appear in the main panel. The simulation runs on the hosting server, not on the phone.

## Scope

Static, modal, applied-force dynamic and corotational P-Delta analyses are computed. Custom sensor diagnosis, sensor placement and reference buckling factors are not independently reproduced here. The section remains uncracked, elastic and perfectly bonded. Deterioration percentages are assumptions, not exposure-time predictions. This is an educational analysis and verification app.

## Run locally

```bash
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py
```

## Sources

- https://opensees.github.io/OpenSeesDocumentation/user/pythonInstall.html
- https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy
