function verifyFrameLab2D
% Optional local MATLAB verification for the public connected release.
% This opens a temporary GUI, checks callbacks, then closes it.
FrameLab2D.verifyNewtonRaphson();
FrameLab2D.verifyEnvironment();
FrameLab2D.verifyConnectedWorkflow();
app=FrameLab2D;
cleanup=onCleanup(@() delete(app)); %#ok<NASGU>
assert(isvalid(app.UIFigure),'App did not open.');
assert(strcmp(app.UIFigure.Name,'FrameLab 2D'),'Unexpected app title.');
click(app,'LOAD RC DEMO');
click(app,'RUN STATIC ANALYSIS');
assert(contains(app.ModelStatus.Text,'Roof mean X'),'Static response missing.');
assert(norm(app.ModelAxes.View-[0 90])<1e-9,'Frame is not in front view.');
assert(~isempty(findall(app.ModelAxes,'Type','line')),'Frame display is missing.');
click(app,'SENSOR LIMITATIONS');
assert(contains(lower(app.ObservabilityStatus.Text),'rank'),'Sensor rank missing.');
click(app,'SIMULATE DATA');
assert(contains(app.DiagnosisStatus.Text,'Synthetic data ready'),'Simulation failed.');
click(app,'RUN DIAGNOSIS');
assert(~isempty(app.DiagnosisTable.Data),'Diagnosis returned no candidates.');
assert(~isequal(app.DiagnosisAxes,app.DiagnosisScanAxes),'Frame and scan axes overlap.');
assert(isempty(findall(app.DiagnosisAxes,'Type','bar')),'Trial bars overlap the frame.');
assert(~isempty(findall(app.DiagnosisScanAxes,'Type','bar')),'Trial chart is missing.');
assert(~isempty(findall(app.DiagnosisAxes,'Tag','FixedSupportPlate')),'Fixed supports are missing.');
click(app,'OPTIMIZE SENSOR LOCATIONS');
assert(~isempty(app.RecommendationTable.Data),'No recommendations returned.');
assert(norm(app.RecommendationAxes.View-[0 90])<1e-9,'Sensor plot is not a front view.');
% Apply deterioration to one member and check unaffected members and repeatability.
click(app,'RC 400 mm PRESET');
mt=findall(app.UIFigure,'Type','uitable','Tag','ExposureMembers');
assert(numel(mt)==1,'Exposure member table missing.');
d=mt.Data;d(1,4:6)={'Combined',10,20};mt.Data=d;
click(app,'APPLY AND COMPARE');
st=findall(app.UIFigure,'Type','uitable','Tag','ExposureSections');
ct=findall(app.UIFigure,'Type','uitable','Tag','ExposureComparison');
assert(height(st.Data)==size(d,1),'Section calculations missing.');
assert(st.Data.EA_ratio(1)<1 && st.Data.EI_ratio(1)<1,'Assigned losses not applied.');
assert(all(st.Data.EA_ratio(2:end)==1),'Unselected members changed.');
assert(height(ct.Data)>=3,'Response comparison missing.');
old=st.Data;click(app,'APPLY AND COMPARE');
assert(isequal(old,st.Data),'Repeated application compounded deterioration.');
click(app,'CALCULATION DETAILS');
au=findall(app.UIFigure,'Tag','ExposureAudit');
assert(~isempty(au) && any(contains(au.Value,'Healthy kLocal')),'Element matrix audit missing.');
click(app,'RE-OPTIMIZE SENSORS');
assert(~isempty(app.RecommendationTable.Data),'Exposure re-optimization returned no rows.');
click(app,'RESTORE HEALTHY RC');click(app,'RUN STATIC ANALYSIS');
assert(contains(app.ModelStatus.Text,'Healthy RC'),'Healthy RC was not restored.');
fprintf('PASS: connected numerical checks and MATLAB GUI smoke checks.\n');
end

function click(app,label)
b=findall(app.UIFigure,'Type','uibutton','Text',label);
assert(~isempty(b),'Missing button: %s',label);
callback=b(1).ButtonPushedFcn;callback(b(1),[]);
end
