# Try FrameLab 2D in one minute

With MATLAB already available, open the extracted app folder and run:

```matlab
FrameLab2D
```

1. Keep the default frame inputs and choose **RUN STATIC ANALYSIS** in **2D Frame Model**.
2. Inspect the orange deformed frame and the fixed supports.
3. Open **Damage Diagnosis**, keep the default hidden element and severity, choose **SIMULATE DATA**, then **RUN DIAGNOSIS**.
4. Read the best element/loss estimate and inspect the separate **Frame and supports** and **Single-member trial chart** views.

The measurements are simulated. Noise can change the estimated loss and candidate ranking. Each trial is a separate single-member damage hypothesis.

Optional RC example: in **Exposure scenarios**, choose **RC 400 mm PRESET**, enter concrete E loss = 10% and steel area loss = 20%, choose the intended member group, **ASSIGN TO GROUP**, then **APPLY AND COMPARE**. Read EA/EI ratios in **Section calculations** and displacement comparisons in **Response and hand check**.
