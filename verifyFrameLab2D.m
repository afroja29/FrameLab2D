function verifyFrameLab2D
% Run from this folder in MATLAB. Opens and closes a temporary app instance.
FrameLab2D.verifyNewtonRaphson();
FrameLab2D.verifyEnvironment();
app = FrameLab2D;
cleanup = onCleanup(@() delete(app)); %#ok<NASGU>
assert(isvalid(app.UIFigure),'App did not open.');
assert(strcmp(app.UIFigure.Name,'FrameLab 2D'),'App name was not updated.');
assert(strcmp(app.DesignTab4.Title,'1. 2D Frame Model'));
assert(norm(app.ModelAxes.View-[0 90])<1e-9,'Model is not in front view.');
app.TabGroup.SelectedTab = app.DesignTab4;
click(app,'RUN STATIC ANALYSIS');
assert(contains(app.ModelStatus.Text,'Healthy analysis complete'));
% Static analysis must actually draw the requested orange deformation overlay.
lines = findall(app.ModelAxes,'Type','line');
assert(~isempty(lines),'Static deformation overlay is missing.');
click(app,'SIMULATE DATA');
assert(contains(app.DiagnosisStatus.Text,'Synthetic data ready'));
click(app,'RUN DIAGNOSIS');
assert(~isempty(app.DiagnosisTable.Data),'Diagnosis returned no candidates.');
assert(~isequal(app.DiagnosisAxes,app.DiagnosisScanAxes),'Frame and trial chart share axes.');
assert(isempty(findall(app.DiagnosisAxes,'Type','bar')),'Trial bars overlap the frame.');
assert(~isempty(findall(app.DiagnosisScanAxes,'Type','bar')),'Trial chart is missing.');
assert(~isempty(findall(app.DiagnosisAxes,'Tag','FixedSupportPlate')),'Fixed support plates are missing.');
assert(~isempty(findall(app.DiagnosisAxes,'Tag','FixedSupportLabel')),'Fixed support labels are missing.');
click(app,'OPTIMIZE SENSOR LOCATIONS');
assert(~isempty(app.RecommendationTable.Data),'No recommendations returned.');
assert(norm(app.RecommendationAxes.View-[0 90])<1e-9,'Recommended sensors are not in front view.');
assert(isempty(findall(app.UIFigure,'Type','uitab','Title','5. Observability')),'Separate observability tab still exists.');
app.TabGroup.SelectedTab=app.DesignTab5;
click(app,'CHECK SENSOR ADEQUACY');
assert(contains(app.ObservabilityStatus.Text,'rank'));
assert(contains(app.ObservabilityStatus.Text,'Single-member diagnosis') || contains(app.ObservabilityStatus.Text,'passes'));
assert(isequal(ancestor(app.DesignButton24,'uitab'),app.DesignTab5),'Adequacy check is not in Sensor Network.');
% Exposure comparison, changed properties, reset, and reapplication.
% Original app defaults are steel-like: explicitly select valid RC geometry first.
click(app,'RC 400 mm PRESET');
% Regression: the old input container was only 132 px but needed 156 px.
ig=findall(app.UIFigure,'Tag','ExposureInputGrid');
eg=findall(app.UIFigure,'Tag','ExposureLayout');
assert(numel(ig)==1 && numel(eg)==1,'Exposure layout containers missing.');
assert(strcmp(eg.Scrollable,'on') && strcmp(ig.Scrollable,'on'),'Exposure scrolling disabled.');
required=sum(cell2mat(ig.RowHeight))+ig.Padding(2)+ig.Padding(4)+3*ig.RowSpacing;
assert(eg.RowHeight{2}>=required,'Exposure input rows exceed allocated height.');
% Exercise smaller laptop windows; the controls must retain their minimum sizes.
oldPosition=app.UIFigure.Position;
for sizeWH=[900 640;1100 720].'
    position=oldPosition; position(3:4)=sizeWH.';
    app.UIFigure.Position=position; drawnow;
    assert(strcmp(eg.Scrollable,'on') && all(cell2mat(ig.RowHeight)>=28));
end
app.UIFigure.Position=oldPosition; drawnow;
mt=findall(app.UIFigure,'Type','uitable','Tag','ExposureMembers');
assert(numel(mt)==1,'Exposure member table is missing.');
d=mt.Data; d(1,4:6)={'Combined',10,20}; mt.Data=d;
click(app,'APPLY AND COMPARE');
st=findall(app.UIFigure,'Type','uitable','Tag','ExposureSections');
ct=findall(app.UIFigure,'Type','uitable','Tag','ExposureComparison');
assert(height(st.Data)==size(d,1),'No section calculations returned.');
assert(st.Data.EA_ratio(1)<1 && st.Data.EI_ratio(1)<1,'Deterioration was not applied.');
assert(all(st.Data.EA_ratio(2:end)==1),'Unselected members changed.');
assert(height(ct.Data)==3,'Static comparison is missing.');
a=st.Data;
click(app,'APPLY AND COMPARE');
assert(isequal(a,st.Data),'Repeated application compounded deterioration.');
click(app,'HAND CHECK');
au=findall(app.UIFigure,'Tag','ExposureAudit');
assert(any(contains(au.Value,'Healthy kLocal')),'Member stiffness audit missing.');
click(app,'RESTORE ORIGINAL');
click(app,'RUN STATIC ANALYSIS');
assert(contains(app.ModelStatus.Text,'Healthy analysis complete'),'Original model was not restored.');
fprintf('PASS: startup, static display, simulation, diagnosis, recommendation, observability, exposure and reset.\n');
end

function click(app,label)
b = findall(app.UIFigure,'Type','uibutton','Text',label);
assert(~isempty(b),'Missing button: %s',label);
callback = b(1).ButtonPushedFcn;
callback(b(1),[]);
end
