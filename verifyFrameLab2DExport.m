function report=verifyFrameLab2DExport(folder)
% Recalculate objective/RMSE from an exported folder without guessing observations.
if nargin==0,folder=uigetdir(pwd,'Select a FrameLab2D_audit folder');end
assert(~isequal(folder,0),'No audit folder selected.');
measurements=readtable(fullfile(folder,'Measurements_exact.csv'));
p=readtable(fullfile(folder,'Candidate_readings.csv'));
c=readtable(fullfile(folder,'Diagnosis_candidates.csv'));
J=zeros(height(c),1);rmse=J;
for k=1:height(c)
    ix=p.Element==c.Element(k) & abs(p.Loss-c.Loss(k))<1e-12;
    assert(nnz(ix)==height(measurements),'Incomplete candidate readings.');
    for row=find(ix).'
        j=find(measurements.DOF==p.DOF(row));assert(numel(j)==1,'DOF mismatch.');
        assert(abs(measurements.MeasuredReading(j)-p.Measured(row))<1e-12,'Measurement mismatch.');
        assert(abs(measurements.Sigma(j)-p.Sigma(row))<1e-12,'Sigma mismatch.');
    end
    diff=p.Measured(ix)-p.Predicted(ix);
    J(k)=sum((diff./p.Sigma(ix)).^2);rmse(k)=sqrt(mean(diff.^2));
end
assert(all(abs(J-c.Objective)<1e-7*max(1,abs(c.Objective))),'Objective export verification failed.');
assert(all(abs(rmse-c.RMSE)<1e-10+1e-7*abs(c.RMSE)),'RMSE export verification failed.');
report=table(c.Element,c.Loss,c.Objective,J,c.RMSE,rmse, ...
    'VariableNames',{'Element','Loss','StoredObjective','RecomputedObjective','StoredRMSE','RecomputedRMSE'});
disp(report);fprintf('PASS: exact exported measurements reproduce all candidate objectives and RMSE.\n');
end
