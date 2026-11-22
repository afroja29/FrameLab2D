% FrameLab 2D - planar frame analysis, damage studies and sensor placement.
% Run as FrameLab2D from the folder containing FrameLab2D.m.
% Static identification and dynamics remain linear elastic. Stability equilibrium
% uses an analytic corotational element and full Newton-Raphson with load stepping
% and backtracking; no plasticity, nonlinear concrete or arc-length postbuckling.
% Eigenbuckling factors remain the reference linearized gravity-prestress results.
% Kt*du = Fext - fint(u); update u, fint and Kt every Newton iteration.
% Headless check: FrameLab2D.verifyNewtonRaphson
% The formulation's independent Python implementation was compared against
% OpenSees 3.8.0 elasticBeamColumn/Corotational/Newton on seven cases;
% max relative displacement difference 1.82e-9. MATLAB GUI was not run here.
% References: https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/algorithm/Newton.html
% https://openseespydoc.readthedocs.io/en/latest/src/corotational.html
classdef FrameLab2D < matlab.apps.AppBase
    %FRAMELAB2D Structural analysis, deterioration assessment and sensor placement.

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        DesignGridLayout1         matlab.ui.container.GridLayout
        TabGroup                  matlab.ui.container.TabGroup
        DesignTab4                matlab.ui.container.Tab
        DesignGridLayout13        matlab.ui.container.GridLayout
        ModelStatus               matlab.ui.control.Label
        DesignGridLayout20        matlab.ui.container.GridLayout
        DesignButton37            matlab.ui.control.Button
        DesignButton36            matlab.ui.control.Button
        DesignButton35            matlab.ui.control.Button
        DesignGridLayout19        matlab.ui.container.GridLayout
        DesignLabel34             matlab.ui.control.Label
        LoadField                 matlab.ui.control.NumericEditField
        DesignLabel33             matlab.ui.control.Label
        IField                    matlab.ui.control.NumericEditField
        DesignLabel32             matlab.ui.control.Label
        AField                    matlab.ui.control.NumericEditField
        DesignLabel31             matlab.ui.control.Label
        EField                    matlab.ui.control.NumericEditField
        DesignLabel30             matlab.ui.control.Label
        StoryHeightField          matlab.ui.control.NumericEditField
        DesignLabel29             matlab.ui.control.Label
        BayWidthField             matlab.ui.control.NumericEditField
        DesignLabel28             matlab.ui.control.Label
        BaysField                 matlab.ui.control.NumericEditField
        DesignLabel27             matlab.ui.control.Label
        StoriesField              matlab.ui.control.NumericEditField
        DesignLabel26             matlab.ui.control.Label
        ModelAxes                 matlab.ui.control.UIAxes
        DesignTab5                matlab.ui.container.Tab
        DesignGridLayout14        matlab.ui.container.GridLayout
        SensorTable               matlab.ui.control.Table
        DesignGridLayout21        matlab.ui.container.GridLayout
        DesignButton42            matlab.ui.control.Button
        DesignButton41            matlab.ui.control.Button
        SensorNoiseField          matlab.ui.control.NumericEditField
        DesignLabel40             matlab.ui.control.Label
        SensorDirectionDropDown   matlab.ui.control.DropDown
        DesignLabel39             matlab.ui.control.Label
        SensorNodeField           matlab.ui.control.NumericEditField
        DesignLabel38             matlab.ui.control.Label
        SensorAxes                matlab.ui.control.UIAxes
        DesignTab6                matlab.ui.container.Tab
        DesignGridLayout15        matlab.ui.container.GridLayout
        DiagnosisStatus           matlab.ui.control.Label
        DiagnosisTable            matlab.ui.control.Table
        DesignGridLayout22        matlab.ui.container.GridLayout
        DesignButton47            matlab.ui.control.Button
        DesignButton46            matlab.ui.control.Button
        SimulationNoiseField      matlab.ui.control.NumericEditField
        DesignLabel45             matlab.ui.control.Label
        DamageValueLabel          matlab.ui.control.Label
        DamageSlider              matlab.ui.control.Slider
        DesignLabel44             matlab.ui.control.Label
        HiddenElementField        matlab.ui.control.NumericEditField
        DesignLabel43             matlab.ui.control.Label
        DiagnosisAxes             matlab.ui.control.UIAxes
        DiagnosisScanAxes         matlab.ui.control.UIAxes
        DesignTab7                matlab.ui.container.Tab
        DesignGridLayout16        matlab.ui.container.GridLayout
        DynamicStatus             matlab.ui.control.Label
        DesignGridLayout23        matlab.ui.container.GridLayout
        DesignGridLayout53        matlab.ui.container.GridLayout
        DesignButton57            matlab.ui.control.Button
        DesignButton56            matlab.ui.control.Button
        TimeStepField             matlab.ui.control.NumericEditField
        DurationField             matlab.ui.control.NumericEditField
        ExcitationFrequencyField  matlab.ui.control.NumericEditField
        DampingField              matlab.ui.control.NumericEditField
        FloorMassField            matlab.ui.control.NumericEditField
        DesignLabel52             matlab.ui.control.Label
        DesignLabel51             matlab.ui.control.Label
        DesignLabel50             matlab.ui.control.Label
        DesignLabel49             matlab.ui.control.Label
        DesignLabel48             matlab.ui.control.Label
        TimeHistoryAxes           matlab.ui.control.UIAxes
        DynamicAxes               matlab.ui.control.UIAxes
        ObservabilityStatus       matlab.ui.control.Label
        DesignButton24            matlab.ui.control.Button
        DesignTab9                matlab.ui.container.Tab
        DesignGridLayout18        matlab.ui.container.GridLayout
        RecommendationStatus      matlab.ui.control.Label
        RecommendationTable       matlab.ui.control.Table
        DesignGridLayout25        matlab.ui.container.GridLayout
        DesignButton55            matlab.ui.control.Button
        NumberSensorsField        matlab.ui.control.NumericEditField
        DesignLabel54             matlab.ui.control.Label
        RecommendationAxes        matlab.ui.control.UIAxes
        DesignPanel2              matlab.ui.container.Panel
        DesignGridLayout3         matlab.ui.container.GridLayout
        DesignLabel12             matlab.ui.control.Label
        DesignLabel11             matlab.ui.control.Label
        DesignLabel10             matlab.ui.control.Label
    end

    properties (Access = private)
        Model
        AnalysisResult
        ExistingSensors double = zeros(0,3)
        Measurements double = zeros(0,1)
        MeasurementSigma double = zeros(0,1)
        Diagnosis
        Recommendation
        StabilityTab
        StabilityAxes
        StabilityCurveAxes
        StabilityTable
        StabilityStatus
        GravityField
        GravityFactorField
        StiffnessField
        MeshField
        StabilityView
        StabilityDataView
        StabilityResult
        EnvTab
        EnvGroup
        EnvStory
        EnvExposure
        EnvEcLoss
        EnvAsLoss
        EnvEs
        EnvBars
        EnvDiameter
        EnvCover
        EnvMembers
        EnvSections
        EnvComparison
        EnvAudit
        EnvStatus
        EnvOriginal
        EnvHealthy
        EnvReport
        EnvResponse
        EnvStabilitySettings
        EnvApplied logical = false
        EnvDirty logical = true
        EnvSelected double = []
        StopRequested logical = false
    end

    methods (Access = private)
        function createEnvironmentTab(app)
            app.EnvTab=uitab(app.TabGroup,'Title','7. Exposure scenarios');
            g=uigridlayout(app.EnvTab,[5 1]);
            g.Tag='ExposureLayout';
            % Keep a minimum usable content size. Smaller windows scroll instead
            % of clipping controls or collapsing the result table.
            g.Scrollable='on'; g.Padding=[8 8 8 8]; g.RowSpacing=8;
            g.ColumnWidth={'1x'};
            g.RowHeight={54,166,320,50,76};
            note=uilabel(g,'WordWrap','on','Text', ...
                ['RC rectangular sections, symmetric top/bottom bars, uncracked perfect bond. ' ...
                 'Tab 1 E is concrete Ec here. Use RC 400 mm PRESET first for the demo. Loss percentages are user assumptions.']);
            note.Layout.Row=1;
            c=uigridlayout(g,[4 6]); c.Layout.Row=2;
            c.Tag='ExposureInputGrid'; c.Scrollable='on';
            c.Padding=[8 8 8 8]; c.RowSpacing=6; c.ColumnSpacing=6;
            % 4 rows + padding + spacing need 156 px, within the 166 px
            % parent row. Fixed minimum widths allow horizontal scrolling.
            c.RowHeight={28,28,28,38};
            c.ColumnWidth={145,160,145,130,160,190};
            label=uilabel(c,'Text','Member group'); label.Layout.Row=1; label.Layout.Column=1;
            app.EnvGroup=uidropdown(c,'Items',{'All members','Bottom columns','Exterior columns', ...
                'Columns at story','Beams at story','Selected table rows'});
            app.EnvGroup.Layout.Row=1; app.EnvGroup.Layout.Column=2;
            label=uilabel(c,'Text','Story'); label.Layout.Row=1; label.Layout.Column=3;
            app.EnvStory=uieditfield(c,'numeric','Value',1,'Limits',[1 Inf],'RoundFractionalValues','on');
            app.EnvStory.Layout.Row=1; app.EnvStory.Layout.Column=4;
            label=uilabel(c,'Text','Exposure label'); label.Layout.Row=1; label.Layout.Column=5;
            app.EnvExposure=uidropdown(c,'Items',{'Reference','Sulfate','Chloride / seawater','Combined'});
            app.EnvExposure.Layout.Row=1; app.EnvExposure.Layout.Column=6;
            label=uilabel(c,'Text','Concrete E loss (%)'); label.Layout.Row=2; label.Layout.Column=1;
            app.EnvEcLoss=uieditfield(c,'numeric','Value',0,'Limits',[0 80]);
            app.EnvEcLoss.Layout.Row=2; app.EnvEcLoss.Layout.Column=2;
            label=uilabel(c,'Text','Steel area loss (%)'); label.Layout.Row=2; label.Layout.Column=3;
            app.EnvAsLoss=uieditfield(c,'numeric','Value',0,'Limits',[0 80]);
            app.EnvAsLoss.Layout.Row=2; app.EnvAsLoss.Layout.Column=4;
            label=uilabel(c,'Text','Steel Es (GPa)'); label.Layout.Row=2; label.Layout.Column=5;
            app.EnvEs=uieditfield(c,'numeric','Value',200,'Limits',[1 Inf], ...
                'ValueChangedFcn',@(~,~)app.environmentEdited());
            app.EnvEs.Layout.Row=2; app.EnvEs.Layout.Column=6;
            label=uilabel(c,'Text','Bars per face'); label.Layout.Row=3; label.Layout.Column=1;
            app.EnvBars=uieditfield(c,'numeric','Tag','ExposureBars','Value',2,'Limits',[1 20],'RoundFractionalValues','on', ...
                'ValueChangedFcn',@(~,~)app.environmentEdited());
            app.EnvBars.Layout.Row=3; app.EnvBars.Layout.Column=2;
            label=uilabel(c,'Text','Bar diameter (mm)'); label.Layout.Row=3; label.Layout.Column=3;
            app.EnvDiameter=uieditfield(c,'numeric','Tag','ExposureDiameter','Value',16,'Limits',[1 Inf], ...
                'ValueChangedFcn',@(~,~)app.environmentEdited());
            app.EnvDiameter.Layout.Row=3; app.EnvDiameter.Layout.Column=4;
            label=uilabel(c,'Text','Bar-centre cover (mm)'); label.Layout.Row=3; label.Layout.Column=5;
            app.EnvCover=uieditfield(c,'numeric','Tag','ExposureCover','Value',50,'Limits',[1 Inf], ...
                'ValueChangedFcn',@(~,~)app.environmentEdited());
            app.EnvCover.Layout.Row=3; app.EnvCover.Layout.Column=6;
            assign=uibutton(c,'Text','ASSIGN TO GROUP','ButtonPushedFcn',@(~,~)app.assignEnvironment());
            assign.Layout.Row=4; assign.Layout.Column=[1 2];
            preset=uibutton(c,'Text','RC 400 mm PRESET', ...
                'ButtonPushedFcn',@(~,~)app.environmentRCPreset());
            preset.Layout.Row=4; preset.Layout.Column=[3 4];
            n=uilabel(c,'Text','Cover is to BAR CENTRE. Apply preset BEFORE assigning losses.', ...
                'WordWrap','on'); n.Layout.Row=4; n.Layout.Column=[5 6];
            tg=uitabgroup(g); tg.Layout.Row=3;
            t=uitab(tg,'Title','Member scenarios'); app.EnvMembers=uitable(t,'Tag','ExposureMembers','Units','normalized','Position',[0 0 1 1]);
            app.EnvMembers.ColumnName={'Element','Type','Story','Exposure','Ec loss (%)','As loss (%)'};
            app.EnvMembers.ColumnEditable=[false false false true true true];
            app.EnvMembers.ColumnFormat={'numeric','char','numeric', ...
                {'Reference','Sulfate','Chloride / seawater','Combined'},'numeric','numeric'};
            app.EnvMembers.RowName=[];
            app.EnvMembers.CellSelectionCallback=@(~,ev)app.selectEnvironmentRows(ev);
            app.EnvMembers.CellEditCallback=@(~,~)app.environmentEdited();
            t=uitab(tg,'Title','Section calculations');
            app.EnvSections=uitable(t,'Tag','ExposureSections','Units','normalized','Position',[0 0 1 1],'RowName',[]);
            t=uitab(tg,'Title','Response and hand check');
            ag=uigridlayout(t,[2 1]); ag.RowHeight={125,'1x'};
            app.EnvComparison=uitable(ag,'Tag','ExposureComparison','RowName',[]); app.EnvComparison.Layout.Row=1;
            app.EnvAudit=uitextarea(ag,'Tag','ExposureAudit','Editable','off','FontName','Courier New','FontSize',11);
            app.EnvAudit.Layout.Row=2; app.EnvAudit.Value={'Apply a scenario, select a member row, then click HAND CHECK.'};
            b=uigridlayout(g,[1 5]); b.Layout.Row=4; b.Padding=[0 0 0 0];
            b.ColumnWidth={'1x','1x','1x','1x','1x'};
            uibutton(b,'Text','APPLY AND COMPARE','ButtonPushedFcn',@(~,~)app.applyEnvironment());
            uibutton(b,'Text','COMPARE P-DELTA','ButtonPushedFcn',@(~,~)app.compareEnvironmentStability());
            uibutton(b,'Text','HAND CHECK','ButtonPushedFcn',@(~,~)app.auditEnvironment());
            uibutton(b,'Text','EXPORT AUDIT','ButtonPushedFcn',@(~,~)app.exportEnvironment());
            uibutton(b,'Text','RESTORE ORIGINAL','ButtonPushedFcn',@(~,~)app.restoreEnvironment());
            app.EnvStatus=uitextarea(g,'Tag','ExposureStatus','Editable','off'); app.EnvStatus.Layout.Row=5;
            app.EnvStatus.Value={'Generate the frame in Tab 1. Then assign assumptions to member groups.'};
        end

        function environmentRCPreset(app)
            % Explicit opt-in: preserve geometry/load inputs, reset exposure cases.
            app.EField.Value=30;
            app.AField.Value=.4*.4;
            app.IField.Value=.4^4/12;
            app.EnvEs.Value=200; app.EnvBars.Value=2;
            app.EnvDiameter.Value=16; app.EnvCover.Value=50;
            app.generateModel();
            app.invalidateEnvironmentResults();
            app.TabGroup.SelectedTab=app.EnvTab;
            app.EnvStatus.Value={'RC preset: 400 x 400 mm, Ec=30 GPa, Es=200 GPa, 2 bars/face, d=16 mm, bar-centre cover=50 mm.'; ...
                'Member exposure rows reset. Now assign an exposure and losses, then APPLY AND COMPARE.'};
        end

        function resetEnvironment(app)
            app.EnvOriginal=app.Model;
            app.EnvHealthy=[]; app.EnvReport=table(); app.EnvResponse=table(); app.EnvStabilitySettings=[];
            app.EnvApplied=false; app.EnvDirty=true; app.EnvSelected=[];
            ne=size(app.Model.Elements,1); data=cell(ne,6);
            for e=1:ne
                row=app.Model.Elements(e,:);
                y=app.Model.Nodes(row(2:3),3);
                story=max(1,round(max(y)/app.Model.StoryHeight));
                if row(7)==1, kind='Column'; else, kind='Beam'; end
                data(e,:)={row(1),kind,story,'Reference',0,0};
            end
            app.EnvMembers.Data=data;
            app.EnvSections.Data=table(); app.EnvComparison.Data=table();
            app.EnvAudit.Value={'The original gross-section frame is active. APPLY constructs an RC composite baseline and a deteriorated model.'};
            app.EnvStory.Value=1;
            row=app.Model.Elements(1,:); hh=sqrt(12*row(6)/row(5)); bb=row(5)/hh;
            app.EnvStatus.Value={sprintf('Current section: b=%.3f mm, h=%.3f mm; Tab 1 E=%.3g GPa (interpreted as Ec in this module).', ...
                1000*bb,1000*hh,row(4)/1e9); ...
                'For the RC demo, click RC 400 mm PRESET first. This explicitly updates E/A/I and resets exposure rows.'};
        end

        function selectEnvironmentRows(app,ev)
            if isempty(ev.Indices), app.EnvSelected=[];
            else, app.EnvSelected=unique(ev.Indices(:,1)); end
        end

        function environmentEdited(app)
            app.EnvDirty=true;
            app.EnvStatus.Value={'Inputs changed. Click APPLY AND COMPARE to update results.'; ...
                'Existing results in other tabs still describe the last APPLIED model.'};
        end

        function assignEnvironment(app)
            try
                assert(~isempty(app.EnvOriginal),'Generate a frame first.');
                data=app.EnvMembers.Data; ne=size(data,1);
                story=app.EnvStory.Value;
                switch app.EnvGroup.Value
                    case 'All members', mask=true(ne,1);
                    case 'Bottom columns', mask=strcmp(data(:,2),'Column') & cell2mat(data(:,3))==1;
                    case 'Columns at story', mask=strcmp(data(:,2),'Column') & cell2mat(data(:,3))==story;
                    case 'Beams at story', mask=strcmp(data(:,2),'Beam') & cell2mat(data(:,3))==story;
                    case 'Exterior columns'
                        el=app.EnvOriginal.Elements; xy=app.EnvOriginal.Nodes;
                        x=xy(el(:,2),2); tol=1e-9*max(1,max(abs(xy(:,2))));
                        mask=el(:,7)==1 & (abs(x-min(xy(:,2)))<tol | abs(x-max(xy(:,2)))<tol);
                    case 'Selected table rows'
                        mask=false(ne,1); mask(app.EnvSelected)=true;
                end
                assert(any(mask),'No members match. Select table rows or a valid story.');
                ec=app.EnvEcLoss.Value; ar=app.EnvAsLoss.Value;
                % Enforce each simplified mechanism; combined is explicit.
                switch app.EnvExposure.Value
                    case 'Reference', ec=0; ar=0;
                    case 'Sulfate', ar=0;
                    case 'Chloride / seawater', ec=0;
                end
                for e=find(mask).'
                    data(e,4:6)={app.EnvExposure.Value,ec,ar};
                end
                app.EnvMembers.Data=data; app.environmentEdited();
                app.EnvStatus.Value={sprintf('Assigned %s to %d members. Click APPLY AND COMPARE.',app.EnvExposure.Value,nnz(mask)); ...
                    'Assumed deterioration applies over the full selected member length. Exposure labels do not calculate corrosion rate.'};
            catch ME, app.showError(ME); end
        end

        function applyEnvironment(app)
            try
                opts=struct('Es',app.EnvEs.Value*1e9,'nFace',app.EnvBars.Value, ...
                    'd',app.EnvDiameter.Value/1000,'a',app.EnvCover.Value/1000);
                [healthy,scenario,report]=sdt_environmentBuild(app.EnvOriginal,app.EnvMembers.Data,opts);
                rh=sdt_solveFrame(healthy,[]); rs=sdt_solveFrame(scenario,[]);
                comparison=sdt_environmentCompare(healthy,rh,rs);
                % Commit only after validation and both solves have succeeded.
                app.EnvHealthy=healthy; app.EnvReport=report; app.EnvResponse=comparison; app.EnvStabilitySettings=[];
                app.Model=scenario; app.AnalysisResult=rs;
                app.EnvApplied=true; app.EnvDirty=false;
                app.EnvSections.Data=report; app.EnvSections.ColumnName=report.Properties.VariableNames;
                app.EnvComparison.Data=comparison; app.EnvComparison.ColumnName=comparison.Properties.VariableNames;
                app.invalidateEnvironmentResults();
                severity=1-report.EI_ratio;
                app.plotStructure(app.ModelAxes,rs,[],[],severity,true);
                app.plotStructure(app.SensorAxes,rs,sdt_sensorDOFs(app.ExistingSensors(:,1:2)),[],severity);
                app.ModelStatus.Text='Exposure scenario active | static, dynamics and stability use updated EA and EI.';
                app.EnvStatus.Value={'Applied. Comparisons use HEALTHY RC COMPOSITE versus deteriorated RC, not the original gross-concrete model.'; ...
                    'Uncracked perfect bond; symmetric uniform bar loss; removed steel remains void. Diagnosis now tests ADDITIONAL stiffness loss relative to this scenario.'};
                app.auditEnvironment();
            catch ME, app.showError(ME); end
        end

        function invalidateEnvironmentResults(app)
            app.StopRequested=true;
            app.clearInferenceState(); app.StabilityResult=[];
            for ax={app.DiagnosisAxes,app.DiagnosisScanAxes,app.RecommendationAxes, ...
                    app.DynamicAxes,app.TimeHistoryAxes,app.StabilityAxes,app.StabilityCurveAxes}
                cla(ax{1});
            end
            app.StabilityTable.Data=table();
            app.DiagnosisStatus.Text='Model changed. Simulate new measurements before diagnosis.';
            app.RecommendationStatus.Text='Model changed. Run a new diagnosis before recommendations.';
            app.ObservabilityStatus.Text='Model changed. Run CHECK SENSOR ADEQUACY again.';
            app.DynamicStatus.Text='Model changed. Run a new dynamic simulation.';
            app.StabilityStatus.Value={'Model changed. Run P-Delta again.'};
        end

        function requireEnvironmentApplied(app)
            assert(app.EnvApplied,'Apply an exposure scenario first.');
            assert(~app.EnvDirty,'Exposure inputs changed. Click APPLY AND COMPARE first.');
        end

        function auditEnvironment(app)
            try
                app.requireEnvironmentApplied();
                if isempty(app.EnvSelected), e=1; else, e=app.EnvSelected(1); end
                app.EnvAudit.Value=sdt_environmentAudit(app.EnvHealthy,app.Model,app.EnvReport,e);
            catch ME, app.showError(ME); end
        end

        function compareEnvironmentStability(app)
            try
                app.requireEnvironmentApplied();
                args={app.GravityField.Value*1000,app.GravityFactorField.Value, ...
                    app.StiffnessField.Value,round(app.MeshField.Value)};
                h=sdt_stabilitySolve(app.EnvHealthy,args{:});
                s=sdt_stabilitySolve(app.Model,args{:});
                if isempty(h.Factors), hc=NaN; else, hc=h.Factors(1); end
                if isempty(s.Factors), sc=NaN; else, sc=s.Factors(1); end
                name={'PDelta roof sway'}; unit={'mm'};
                extra=table(name,unit,1000*h.RoofSecond,1000*s.RoofSecond, ...
                    sdt_environmentPercent(h.RoofSecond,s.RoofSecond), ...
                    'VariableNames',{'Metric','Unit','HealthyRC','Scenario','Change_percent'});
                extra=[extra;table({'Elastic critical gravity factor'},{'-'},hc,sc, ...
                    sdt_environmentPercent(hc,sc),'VariableNames',extra.Properties.VariableNames)];
                app.EnvResponse=[app.EnvResponse(1:3,:);extra];
                app.EnvComparison.Data=app.EnvResponse; app.EnvStabilitySettings=args;
                app.StabilityResult=s; app.updateStabilityViews();
                app.EnvStatus.Value={sprintf('Stability compared at %.4g kN/floor, lambda %.4g, EI factor %.3g, %d segments/member.', ...
                    args{1}/1000,args{2},args{3},args{4}); ...
                    'Both responses use corotational Newton-Raphson. Buckling factors remain reference elastic eigenvalues.'};
            catch ME
                % Do not retain obsolete stability comparisons after a failed trial.
                if height(app.EnvResponse)>3
                    app.EnvResponse=app.EnvResponse(1:3,:); app.EnvComparison.Data=app.EnvResponse;
                end
                app.EnvStabilitySettings=[];
                app.EnvStatus.Value={['Stability comparison not completed: ' ME.message]};
                app.showError(ME);
            end
        end

        function restoreEnvironment(app)
            try
                app.Model=app.EnvOriginal; app.AnalysisResult=sdt_solveFrame(app.Model,[]);
                app.resetEnvironment(); app.invalidateEnvironmentResults();
                z=zeros(size(app.Model.Elements,1),1);
                app.plotStructure(app.ModelAxes,app.AnalysisResult,[],[],z,true);
                app.plotStructure(app.SensorAxes,app.AnalysisResult,sdt_sensorDOFs(app.ExistingSensors(:,1:2)),[],z);
                app.ModelStatus.Text='Original gross-section properties restored.';
            catch ME, app.showError(ME); end
        end

        function exportEnvironment(app)
            try
                app.requireEnvironmentApplied();
                folder=uigetdir(pwd,'Choose folder for exposure audit');
                if isequal(folder,0), return; end
                writetable(app.EnvReport,fullfile(folder,'Exposure_section_calculations.csv'));
                writetable(app.EnvResponse,fullfile(folder,'Exposure_response_comparison.csv'));
                healthyModel=app.EnvHealthy; scenarioModel=app.Model;
                sectionReport=app.EnvReport; responseReport=app.EnvResponse;
                stabilitySettings=app.EnvStabilitySettings;
                healthyResult=sdt_solveFrame(healthyModel,[]);
                scenarioResult=sdt_solveFrame(scenarioModel,[]);
                save(fullfile(folder,'Exposure_audit.mat'),'healthyModel','scenarioModel', ...
                    'sectionReport','responseReport','healthyResult','scenarioResult','stabilitySettings');
                fid=fopen(fullfile(folder,'Exposure_hand_check.txt'),'w');
                assert(fid>=0,'Could not open hand-check output.');
                closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
                for e=1:height(app.EnvReport)
                    lines=sdt_environmentAudit(healthyModel,scenarioModel,sectionReport,e);
                    fprintf(fid,'%s\n',lines{:}); fprintf(fid,'\n');
                end
                fprintf(fid,'GLOBAL FREE STIFFNESS MATRICES (SI units)\n');
                fprintf(fid,'Healthy Kff:\n%s\n',mat2str(full(healthyResult.K(healthyResult.FreeDOFs,healthyResult.FreeDOFs)),12));
                fprintf(fid,'Scenario Kff:\n%s\n',mat2str(full(scenarioResult.K(scenarioResult.FreeDOFs,scenarioResult.FreeDOFs)),12));
                app.EnvStatus.Value={['Audit exported to ' folder];'CSV properties, response comparison, MAT models/results, and complete member/global matrix hand check.'};
            catch ME, app.showError(ME); end
        end

        function createStabilityTab(app)
            app.StabilityTab=uitab(app.TabGroup,'Title','6. P-Delta & Buckling');
            g=uigridlayout(app.StabilityTab,[4 2]);
            g.RowHeight={114,'1x',190,80}; g.ColumnWidth={'1x','1x'};
            controls=uigridlayout(g,[3 6]); controls.Layout.Row=1;
            controls.Layout.Column=[1 2];
            controls.RowHeight={26,26,32};
            controls.ColumnWidth={150,100,155,100,150,'1x'};
            label=uilabel(controls,'Text','Gravity / floor (kN)'); label.Layout.Row=1; label.Layout.Column=1;
            app.GravityField=uieditfield(controls,'numeric','Value',1000,'Limits',[0 Inf]);
            app.GravityField.Layout.Row=1; app.GravityField.Layout.Column=2;
            label=uilabel(controls,'Text','Gravity factor lambda'); label.Layout.Row=1; label.Layout.Column=3;
            app.GravityFactorField=uieditfield(controls,'numeric','Value',1,'Limits',[0 Inf]);
            app.GravityFactorField.Layout.Row=1; app.GravityFactorField.Layout.Column=4;
            label=uilabel(controls,'Text','Effective EI / gross EI'); label.Layout.Row=1; label.Layout.Column=5;
            app.StiffnessField=uieditfield(controls,'numeric','Value',1,'Limits',[0.01 1]);
            app.StiffnessField.Layout.Row=1; app.StiffnessField.Layout.Column=6;
            label=uilabel(controls,'Text','Segments / member'); label.Layout.Row=2; label.Layout.Column=1;
            app.MeshField=uieditfield(controls,'numeric','Value',4,'Limits',[1 12],'RoundFractionalValues','on');
            app.MeshField.Layout.Row=2; app.MeshField.Layout.Column=2;
            label=uilabel(controls,'Text','Diagram'); label.Layout.Row=2; label.Layout.Column=3;
            app.StabilityView=uidropdown(controls,'Items',{'First vs second order','Buckling mode 1', ...
                'Buckling mode 2','Moment diagram','Axial-force diagram'}, ...
                'ValueChangedFcn',@(~,~)app.updateStabilityViews());
            app.StabilityView.Layout.Row=2; app.StabilityView.Layout.Column=4;
            label=uilabel(controls,'Text','Result table'); label.Layout.Row=2; label.Layout.Column=5;
            app.StabilityDataView=uidropdown(controls,'Items',{'Story drift','Node response','Member end actions', ...
                'Support reactions','Buckling factors','Elastic K (free)','Reference Kg (free)','Newton convergence','Tangent K (free)'}, ...
                'ValueChangedFcn',@(~,~)app.updateStabilityViews());
            app.StabilityDataView.Layout.Row=2; app.StabilityDataView.Layout.Column=6;
            buttons=uigridlayout(controls,[1 5]); buttons.Layout.Row=3; buttons.Layout.Column=[1 6];
            buttons.Padding=[0 0 0 0]; buttons.ColumnWidth={'1x','1x','1x','1x','1x'};
            uibutton(buttons,'Text','RUN P-DELTA','BackgroundColor',[0.16 .52 .74], ...
                'FontColor',[1 1 1],'ButtonPushedFcn',@(~,~)app.runStability());
            uibutton(buttons,'Text','RC FRAME PRESET','ButtonPushedFcn',@(~,~)app.rcStabilityPreset());
            uibutton(buttons,'Text','SHOW 80% CRITICAL','ButtonPushedFcn',@(~,~)app.showNearCritical());
            uibutton(buttons,'Text','VERIFY BENCHMARKS','ButtonPushedFcn',@(~,~)app.verifyStability());
            uibutton(buttons,'Text','EXPORT RESULTS','ButtonPushedFcn',@(~,~)app.exportStability());
            app.StabilityAxes=uiaxes(g); app.StabilityAxes.Layout.Row=2; app.StabilityAxes.Layout.Column=1;
            app.StabilityCurveAxes=uiaxes(g); app.StabilityCurveAxes.Layout.Row=2; app.StabilityCurveAxes.Layout.Column=2;
            app.StabilityTable=uitable(g); app.StabilityTable.Layout.Row=3; app.StabilityTable.Layout.Column=[1 2];
            app.StabilityTable.ColumnName='auto'; app.StabilityTable.RowName=[];
            app.StabilityStatus=uitextarea(g,'Editable','off','FontSize',11);
            app.StabilityStatus.Layout.Row=4; app.StabilityStatus.Layout.Column=[1 2];
            app.StabilityStatus.Value={'Click RUN P-DELTA. Gravity loads are lumped at each floor joint.'; ...
                'Corotational Newton-Raphson, elastic material; reference eigenbuckling is reported separately.'};
        end

        function rcStabilityPreset(app)
            app.StoriesField.Value=3; app.BaysField.Value=2;
            app.BayWidthField.Value=6; app.StoryHeightField.Value=3.5;
            app.EField.Value=30; app.AField.Value=.4*.4; app.IField.Value=.4^4/12;
            app.LoadField.Value=50; app.GravityField.Value=800;
            app.StiffnessField.Value=.5; app.GravityFactorField.Value=1;
            app.MeshField.Value=4; app.generateModel(); app.runStability();
        end

        function runStability(app)
            try
                app.StabilityStatus.Value={'Solving corotational Newton-Raphson equilibrium and reference eigenbuckling...'}; drawnow;
                app.StabilityResult=sdt_stabilitySolve(app.Model,app.GravityField.Value*1000, ...
                    app.GravityFactorField.Value,app.StiffnessField.Value,round(app.MeshField.Value));
                r=app.StabilityResult;
                if isempty(r.Factors), criticalText='No positive buckling factor for this gravity pattern';
                else, criticalText=sprintf('Reference elastic critical factor = %.5g | applied/reference = %.3f', ...
                    r.Factors(1),r.LoadFactor/r.Factors(1)); end
                app.StabilityStatus.Value={criticalText; ...
                    sprintf('Roof sway: first %.4g mm | second %.4g mm | amplification %.4g | residual %.2e', ...
                    1000*r.RoofFirst,1000*r.RoofSecond,r.Amplification,r.EquilibriumResidual); ...
                    sprintf('Gravity = %.4g kN/floor x %.4g; lateral = %.4g kN, held constant; EI factor = %.3f.', ...
                    r.GravityPerFloor/1000,r.LoadFactor,sum(app.Model.Loads(1:3:end))/1000,r.EIFactor); ...
                    sprintf('Full Newton-Raphson: %d history rows, tolerance %.1e. Reference buckling is not nonlinear collapse capacity.', ...
                    height(r.NewtonHistory),r.NewtonTolerance)};
                app.updateStabilityViews();
            catch ME
                app.StabilityResult=[];
                app.StabilityStatus.Value={ME.message}; app.showError(ME);
            end
        end

        function showNearCritical(app)
            try
                r=sdt_stabilitySolve(app.Model,app.GravityField.Value*1000,0, ...
                    app.StiffnessField.Value,round(app.MeshField.Value));
                if isempty(r.Factors), error('There is no positive critical factor. Apply a gravity load.'); end
                app.GravityFactorField.Value=.8*r.Factors(1); app.runStability();
            catch ME, app.showError(ME); end
        end

        function updateStabilityViews(app)
            if isempty(app.StabilityResult), return; end
            r=app.StabilityResult; ax=app.StabilityAxes;
legend(ax,'off');
cla(ax,'reset');
hold(ax,'on');
            span=max([max(r.Mesh.Nodes(:,2))-min(r.Mesh.Nodes(:,2)),max(r.Mesh.Nodes(:,3)),1]);
            viewChoice=app.StabilityView.Value;
            if strcmp(viewChoice,'First vs second order')
                mag=max(hypot(r.Second(1:3:end),r.Second(2:3:end)));
                scale=.12*span/max(mag,eps);
                h0=sdt_drawFrame(ax,r.Mesh,zeros(size(r.Second)),0,[.3 .3 .3],'--');
                h1=sdt_drawFrame(ax,r.Mesh,r.First,scale,[.12 .40 .85],':');
                h2=sdt_drawFrame(ax,r.Mesh,r.Second,scale,[.9 .2 .15],'-');
                legend(ax,[h0 h1 h2],{'Undeformed','First order','Newton-Raphson'},'Location','best');
                title(ax,sprintf('Common visual scale x %.3g | P-Delta comparison',scale));
            elseif startsWith(viewChoice,'Buckling mode')
                k=sscanf(viewChoice,'Buckling mode %d');
                if k>size(r.Modes,2)
                    title(ax,'No mode available for this gravity pattern');
                else
                    mode=r.Modes(:,k); mag=max(hypot(mode(1:3:end),mode(2:3:end)));
                    h0=sdt_drawFrame(ax,r.Mesh,zeros(size(mode)),0,[.3 .3 .3],'--');
                    h1=sdt_drawFrame(ax,r.Mesh,mode,.12*span/max(mag,eps),[.86 .22 .12],'-');
                    legend(ax,[h0 h1],{'Undeformed','Normalized mode'},'Location','best');
                    title(ax,sprintf('Mode %d | gravity factor %.5g | eigen residual %.2e',k,r.Factors(k),r.EigenResiduals(k)));
                end
            else
                sdt_drawFrame(ax,r.Mesh,zeros(size(r.Second)),0,[.35 .35 .35],'--');
                isMoment=strcmp(viewChoice,'Moment diagram');
                vals=r.SegmentActions;
                if isMoment, magnitude=max(abs(vals(:,[3 6])),[],'all'); else, magnitude=max(abs(vals(:,1))); end
                amplitude=.06*span/max(magnitude,eps);
                for e=1:size(r.Mesh.Elements,1)
                    ij=r.Mesh.Elements(e,2:3); xy=r.Mesh.Nodes(ij,2:3); v=xy(2,:)-xy(1,:);
                    normal=[-v(2) v(1)]/norm(v);
                    if isMoment, q=[vals(e,3);-vals(e,6)]; else, q=-vals(e,1)*ones(2,1); end
                    diagram=xy+amplitude*q*normal;
                    plot(ax,diagram(:,1),diagram(:,2),'-','Color',[.7 .1 .7],'LineWidth',1.6);
                    plot(ax,[xy(:,1) diagram(:,1)].',[xy(:,2) diagram(:,2)].','Color',[.75 .6 .75]);
                end
                if isMoment, title(ax,'Corotational Newton-Raphson end moments | N m');
                else, title(ax,'Corotational axial force | tension positive | N'); end
            end
           sdt_drawSupports(ax,r.Mesh,span);
           xlabel(ax,'Width (m)');
ylabel(ax,'Height (m)');
axis(ax,'equal');
axis(ax,'tight');

xBounds=xlim(ax);
yBounds=ylim(ax);

xPad=.03*max(diff(xBounds),1);
yPad=.03*max(diff(yBounds),1);

xlim(ax,xBounds+[-xPad xPad]);
ylim(ax,yBounds+[-yPad yPad]);

ax.YTick=(0:r.Mesh.nStories)*r.Mesh.StoryHeight;

grid(ax,'on');
hold(ax,'off');
            curveAx=app.StabilityCurveAxes; cla(curveAx); hold(curveAx,'on');
            plot(curveAx,1000*r.CurveFirst,r.CurveFactors,'b:','LineWidth',1.8);
            plot(curveAx,1000*r.CurveSecond,r.CurveFactors,'r-','LineWidth',1.8);
            scatter(curveAx,1000*r.RoofSecond,r.LoadFactor,48,[.2 .2 .2],'filled');
            if ~isempty(r.Factors), yline(curveAx,r.Factors(1),'k--','Reference elastic critical factor'); end
            xlabel(curveAx,'Mean roof horizontal displacement (mm)'); ylabel(curveAx,'Gravity factor lambda');
            title(curveAx,'Gravity load-sway path | fixed lateral load');
            legend(curveAx,{'First order','Newton-Raphson','Applied state'},'Location','best');
            grid(curveAx,'on'); hold(curveAx,'off');
            switch app.StabilityDataView.Value
                case 'Story drift', data=r.StoryDrift;
                case 'Node response', data=r.NodeResponse;
                case 'Member end actions', data=r.MemberActions;
                case 'Support reactions', data=r.Reactions;
                case 'Buckling factors'
                    data=table((1:numel(r.Factors)).',r.Factors,r.EigenResiduals, ...
                        'VariableNames',{'Mode','CriticalGravityFactor','EigenResidual'});
                case 'Elastic K (free)', data=array2table(full(r.K(r.Free,r.Free)));
                case 'Reference Kg (free)', data=array2table(full(r.Kg(r.Free,r.Free)));
                case 'Newton convergence', data=r.NewtonHistory;
                case 'Tangent K (free)', data=array2table(full(r.TangentK(r.Free,r.Free)));
            end
            app.StabilityTable.Data=data;
        end

        function verifyStability(app)
            try
                report=sdt_verifyStability();
                uialert(app.UIFigure,report,'Stability benchmark results','Icon','success');
            catch ME, app.showError(ME); end
        end

        function exportStability(app)
            if isempty(app.StabilityResult), app.runStability(); end
            if isempty(app.StabilityResult), return; end
            [name,folder]=uiputfile('*.mat','Save stability model, matrices and results','CE7405_stability.mat');
            if isequal(name,0), return; end
            results=app.StabilityResult; %#ok<NASGU>
            save(fullfile(folder,name),'results');
        end

        function bindRuntimeCallbacks(app)
            % Reconnect runtime actions after App Designer regenerates the
            % protected createComponents section from Design View.
            bindings={ ...
                'INSERT',@generateModel; ...
                'RESET 2D VIEW',@reset2DView; ...
                'RUN STATIC ANALYSIS',@analyzeHealthy; ...
                'ADD SENSOR',@addSensor; ...
                'REMOVE SELECTED',@removeSensor; ...
                'SIMULATE DATA',@simulateMeasurements; ...
                'RUN DIAGNOSIS',@runDiagnosis; ...
                'RUN ANIMATION',@runDynamicSimulation; ...
                'STOP',@stopDynamicSimulation; ...
                'CHECK SENSOR ADEQUACY',@runObservability; ...
                'OPTIMIZE SENSOR LOCATIONS',@runRecommendation};
            for k=1:size(bindings,1)
                button=findall(app.UIFigure,'Type','uibutton','Text',bindings{k,1});
                if ~isempty(button)
                    button(1).ButtonPushedFcn=createCallbackFcn(app,bindings{k,2},true);
                end
            end
            app.SensorTable.CellEditCallback=createCallbackFcn(app,@sensorTableEdited,true);
            app.DamageSlider.ValueChangedFcn=createCallbackFcn(app,@damageSliderMoved,true);
        end

        function configureAxesInteractions(app)
            for ax={app.ModelAxes,app.SensorAxes,app.DynamicAxes,app.RecommendationAxes}
                axtoolbar(ax{1},{'pan','zoomin','zoomout','restoreview'});
                view(ax{1},2); ax{1}.Projection='orthographic';
            end
            rotate3d(app.UIFigure,'off');
        end

        function reset2DView(app,varargin)
            rotate3d(app.UIFigure,'off');
            app.plotStructure(app.ModelAxes,app.AnalysisResult,[],[],zeros(size(app.Model.Elements,1),1),false);
            app.ModelStatus.Text='2D front view restored. Sensors and members share one analytical plane.';
        end

        function generateModel(app,varargin)
            try
                app.Model = sdt_generatePortalFrame(app.StoriesField.Value,app.BaysField.Value, ...
                    app.BayWidthField.Value,app.StoryHeightField.Value, ...
                    app.EField.Value*1e9,app.AField.Value,app.IField.Value);
                nCols = app.Model.nBays+1;
                topNodes = app.Model.nStories*nCols + (1:nCols);
                app.Model.Loads(:) = 0;
                app.Model.Loads(3*topNodes-2) = app.LoadField.Value*1e3/numel(topNodes);
                app.resetEnvironment();
                app.AnalysisResult = sdt_solveFrame(app.Model,[]);
                app.StabilityResult = [];
                node1 = nCols + nCols;
                node2 = app.Model.nStories*nCols+nCols;
                app.ExistingSensors = unique([node1,1,app.SimulationNoiseField.Value; ...
                    node2,1,app.SimulationNoiseField.Value],'rows','stable');
                app.HiddenElementField.Limits = [1,size(app.Model.Elements,1)];
                app.HiddenElementField.Value = 1;
                app.SensorNodeField.Limits = [1,size(app.Model.Nodes,1)];
                app.SensorNodeField.Value = node2;
                app.clearInferenceState();
                app.refreshSensorTable();
                app.plotStructure(app.ModelAxes,app.AnalysisResult,[],[],zeros(size(app.Model.Elements,1),1));
                app.plotStructure(app.SensorAxes,app.AnalysisResult,sdt_sensorDOFs(app.ExistingSensors(:,1:2)),[],zeros(size(app.Model.Elements,1),1));
                app.ModelStatus.Text = sprintf('INPUT ACCEPTED: %s | %d nodes | %d elements | cond(Kff)=%.2e', ...
                    app.Model.Name,size(app.Model.Nodes,1),size(app.Model.Elements,1),app.AnalysisResult.ConditionNumber);
            catch ME
                app.showError(ME);
            end
        end

        function analyzeHealthy(app,varargin)
            try
                app.AnalysisResult = sdt_solveFrame(app.Model,[]);
                app.plotStructure(app.ModelAxes,app.AnalysisResult,[],[],zeros(size(app.Model.Elements,1),1),true);
                maxDisp = max(hypot(app.AnalysisResult.Displacements(1:3:end), ...
                    app.AnalysisResult.Displacements(2:3:end)));
                if app.EnvApplied, label='Exposure scenario analysis complete'; else, label='Healthy analysis complete'; end
                app.ModelStatus.Text = sprintf([label ' | max |u| = %.4g m | equilibrium residual = %.3e N'], ...
                    maxDisp,norm(app.AnalysisResult.Reactions(app.AnalysisResult.FreeDOFs)));
            catch ME
                app.showError(ME);
            end
        end

        function addSensor(app,varargin)
            node = round(app.SensorNodeField.Value);
            direction = find(strcmp(app.SensorDirectionDropDown.Items,app.SensorDirectionDropDown.Value),1);
            newRow = [node,direction,app.SensorNoiseField.Value];
            dof = sdt_sensorDOFs(newRow(1:2));
            if ismember(dof,app.Model.FixedDOFs)
                uialert(app.UIFigure,'A restrained DOF cannot be used as a response sensor.','Invalid sensor');
                return
            end
            oldDOFs = sdt_sensorDOFs(app.ExistingSensors(:,1:2));
            if ismember(dof,oldDOFs)
                uialert(app.UIFigure,'That node-direction channel already exists.','Duplicate sensor');
                return
            end
            app.ExistingSensors(end+1,:) = newRow;
            app.clearInferenceState();
            app.refreshSensorTable();
            app.plotStructure(app.SensorAxes,app.AnalysisResult,sdt_sensorDOFs(app.ExistingSensors(:,1:2)),[],zeros(size(app.Model.Elements,1),1));
        end

        function removeSensor(app,varargin)
            if isempty(app.ExistingSensors), return; end
            idx = size(app.ExistingSensors,1);
            if ~isempty(app.SensorTable.Selection)
                idx = app.SensorTable.Selection(1);
            end
            app.ExistingSensors(idx,:) = [];
            app.clearInferenceState();
            app.refreshSensorTable();
            app.plotStructure(app.SensorAxes,app.AnalysisResult,sdt_sensorDOFs(app.ExistingSensors(:,1:2)),[],zeros(size(app.Model.Elements,1),1));
        end

        function sensorTableEdited(app,varargin)
            try
                data = app.SensorTable.Data;
                if istable(data)
                    data = [data.Node,data.Direction,data.NoisePercent];
                end
                validateattributes(data,{'double'},{'2d','ncols',3,'finite'});
                dofs = sdt_sensorDOFs(data(:,1:2));
                if any(data(:,1)>size(app.Model.Nodes,1)) || ...
                        any(ismember(dofs,app.Model.FixedDOFs))
                    error('shm:InvalidSensor','Choose an existing node and an unrestrained direction.');
                end
                if numel(unique(dofs))~=numel(dofs)
                    error('shm:DuplicateSensor','Each node-direction channel must be unique.');
                end
                if any(data(:,3)<0 | data(:,3)>25)
                    error('shm:InvalidNoise','Sensor noise must be between 0 and 25 percent.');
                end
                app.ExistingSensors = data;
                app.clearInferenceState();
                app.refreshSensorTable();
                app.plotStructure(app.SensorAxes,app.AnalysisResult,dofs,[],zeros(size(app.Model.Elements,1),1));
            catch ME
                app.refreshSensorTable();
                app.showError(ME);
            end
        end

        function damageSliderMoved(app,varargin)
            app.DamageValueLabel.Text = sprintf('%.0f%% stiffness loss',app.DamageSlider.Value);
        end

        function simulateMeasurements(app,varargin)
            try
                app.requireSensors();
                nElements = size(app.Model.Elements,1);
                alpha = zeros(nElements,1);
                hidden = round(app.HiddenElementField.Value);
                alpha(hidden) = app.DamageSlider.Value/100;
                truth = sdt_solveFrame(app.Model,alpha);
                dofs = sdt_sensorDOFs(app.ExistingSensors(:,1:2));
                yTrue = truth.Displacements(dofs);
                noisePercent = app.SimulationNoiseField.Value;
                % Combine channel noise with the simulation noise in quadrature.
                % Each channel uses its own response scale (m or rad).
                healthy = app.AnalysisResult.Displacements(dofs);
                reference = max(max(abs(yTrue),abs(healthy)),1e-6);
                percent = hypot(app.ExistingSensors(:,3),noisePercent);
                sigma = max(percent/100.*reference,1e-9);
                stream = RandStream('mt19937ar','Seed',7405);
                app.Measurements = yTrue + sigma.*randn(stream,size(yTrue));
                app.MeasurementSigma = sigma;
                app.Diagnosis = [];
                app.Recommendation = [];
                app.DiagnosisTable.Data = table();
                app.RecommendationTable.Data = table();
                cla(app.DiagnosisScanAxes,'reset');
                app.plotStructure(app.DiagnosisAxes,truth,dofs,[],alpha,true);
                app.DiagnosisStatus.Text = sprintf('Synthetic data ready | true element %d | true loss %.0f%% | noise uses each channel response scale', ...
                    hidden,100*alpha(hidden));
            catch ME
                app.showError(ME);
            end
        end

        function runDiagnosis(app,varargin)
            try
                if isempty(app.Measurements)
                    app.simulateMeasurements();
                end
                if isempty(app.Measurements), return; end
                dofs = sdt_sensorDOFs(app.ExistingSensors(:,1:2));
                app.DiagnosisStatus.Text = 'Running single-damage inverse search...'; drawnow;
                app.Diagnosis = sdt_estimateSingleDamage(app.Model,dofs,app.Measurements,app.MeasurementSigma,linspace(0.01,0.80,80));
                top = app.Diagnosis.Candidates(1:min(8,end),:);
                app.DiagnosisTable.Data = array2table(top,'VariableNames',{'Element','Alpha','Objective','RMSE_m'});
                severityByElement = zeros(size(app.Model.Elements,1),1);
                for i=1:size(app.Diagnosis.Candidates,1)
                    severityByElement(app.Diagnosis.Candidates(i,1)) = app.Diagnosis.Candidates(i,2);
                end
                ax=app.DiagnosisScanAxes;
                cla(ax,'reset'); hold(ax,'off'); legend(ax,'off');
                bar(ax,1:numel(severityByElement),100*severityByElement,'FaceColor',[.18 .42 .62]);
                ax.DataAspectRatioMode='auto'; ax.PlotBoxAspectRatioMode='auto';
                ylim(ax,[0 100]); xlim(ax,[.5 numel(severityByElement)+.5]);
                xlabel(ax,'Assumed damaged element ID'); ylabel(ax,'Best-fit loss for this trial (%)');
                title(ax,{'Separate single-member trials','Bars do not represent simultaneous damage'});
                grid(ax,'on');
                bestAlpha=zeros(size(app.Model.Elements,1),1);
                bestAlpha(app.Diagnosis.BestElement)=app.Diagnosis.BestAlpha;
                fitted=sdt_solveFrame(app.Model,bestAlpha);
                app.plotStructure(app.DiagnosisAxes,fitted,dofs,[],bestAlpha,true);
                title(app.DiagnosisAxes,{'Best-fit damaged frame | fixed bases', ...
                    sprintf('Estimated element %d: %.1f%% stiffness loss', ...
                    app.Diagnosis.BestElement,100*app.Diagnosis.BestAlpha)});
                app.DiagnosisStatus.Text = sprintf('Best estimate: element %d, %.1f%% loss | weighted objective %.3g | RMSE %.3g m', ...
                    app.Diagnosis.BestElement,100*app.Diagnosis.BestAlpha,app.Diagnosis.Objective,app.Diagnosis.RMSE);
            catch ME
                app.showError(ME);
            end
        end

        function runObservability(app,varargin)
            try
                app.requireSensors();
                dofs = sdt_sensorDOFs(app.ExistingSensors(:,1:2));
                [H,~] = sdt_computeSensitivity(app.Model,dofs,[],1e-3);
                columnScale = max(vecnorm(H,2,1),eps);
                Hs = H./columnScale;
                sv = svd(Hs,'econ');
                tol = max(size(Hs))*eps(max(sv));
                numericalRank = sum(sv>tol);
                if numericalRank<size(H,2)
                    app.ObservabilityStatus.Text=sprintf( ...
                        ['%d measurement channels for %d member stiffness parameters. ' ...
                         'Simultaneous damage in all members cannot be independently resolved ' ...
                         '(local rank %d/%d). Single-member diagnosis can still be attempted.'], ...
                         size(H,1),size(H,2),numericalRank,size(H,2));
                else
                    app.ObservabilityStatus.Text=sprintf( ...
                        ['Current layout passes the local rank check (%d/%d). ' ...
                         'This does not guarantee reliable diagnosis under measurement noise or model error.'], ...
                         numericalRank,size(H,2));
                end
            catch ME
                app.showError(ME);
            end
        end

        function runRecommendation(app,varargin)
            try
                if isempty(app.Diagnosis)
                    app.runDiagnosis();
                end
                if isempty(app.Diagnosis), return; end
                dofs = sdt_sensorDOFs(app.ExistingSensors(:,1:2));
                sigma = max(mean(app.MeasurementSigma),1e-9);
                app.Recommendation = sdt_recommendSensors(app.Model,dofs,app.Diagnosis,round(app.NumberSensorsField.Value),sigma);
                if isempty(app.Recommendation.DOFs)
                    app.RecommendationTable.Data = table();
                    app.RecommendationStatus.Text = 'No additional horizontal sensor channels are available.';
                    return
                end
                tableData = [app.Recommendation.Nodes,app.Recommendation.Directions,app.Recommendation.DOFs,app.Recommendation.Scores];
                app.RecommendationTable.Data = array2table(tableData,'VariableNames',{'Node','Direction','DOF','SeparationScore'});
                bestAlpha = zeros(size(app.Model.Elements,1),1);
                bestAlpha(app.Diagnosis.BestElement)=app.Diagnosis.BestAlpha;
                bestResult = sdt_solveFrame(app.Model,bestAlpha);
                app.plotStructure(app.RecommendationAxes,bestResult,dofs,app.Recommendation.DOFs,bestAlpha,true);
                app.RecommendationStatus.Text = sprintf('Recommended node(s): %s | minimum cumulative scenario separation %.3g', ...
                    strjoin(string(app.Recommendation.Nodes.'),', '),app.Recommendation.Scores(end));
            catch ME
                app.showError(ME);
            end
        end

        function modelInputChanged(app,varargin)
            % Structural inputs regenerate immediately, including floor count.
            app.generateModel();
        end

        function stopDynamicSimulation(app,varargin)
            app.StopRequested = true;
            app.DynamicStatus.Text = 'Stopping animation...';
        end

        function runDynamicSimulation(app,varargin)
            try
                app.StopRequested = false;
                app.DynamicStatus.Text = 'Assembling mass, damping, and stiffness matrices...';
                drawnow;

                nd = 3*size(app.Model.Nodes,1);
                fixed = unique(app.Model.FixedDOFs(:));
                free = setdiff((1:nd).',fixed);
                K = app.AnalysisResult.K(free,free);

                massDiagonal = zeros(nd,1);
                nCols = app.Model.nBays+1;
                nodeMass = app.FloorMassField.Value*1000/nCols; % tonne/floor -> kg/node
                for level = 1:app.Model.nStories
                    floorNodes = level*nCols + (1:nCols);
                    massDiagonal(3*floorNodes-2) = nodeMass;
                    massDiagonal(3*floorNodes-1) = nodeMass;
                    massDiagonal(3*floorNodes) = nodeMass*app.StoryHeightField.Value^2/12;
                end
                M = spdiags(massDiagonal(free),0,numel(free),numel(free));

                lambda = eig(full(K),full(M));
                lambda = real(lambda(isfinite(lambda) & real(lambda)>0));
                omega1 = sqrt(min(lambda));
                zeta = app.DampingField.Value/100;
                C = 2*zeta*omega1*M;

                dt = app.TimeStepField.Value;
                time = 0:dt:app.DurationField.Value;
                nt = numel(time);
                U = zeros(numel(free),nt); V = U; Acc = U;
                loadVector = app.Model.Loads(free);
                frequency = app.ExcitationFrequencyField.Value;
                P = loadVector*sin(2*pi*frequency*time);
                Acc(:,1) = M\(P(:,1)-C*V(:,1)-K*U(:,1));

                beta = 1/4; gamma = 1/2;
                a0 = 1/(beta*dt^2); a1 = gamma/(beta*dt);
                a2 = 1/(beta*dt); a3 = 1/(2*beta)-1;
                a4 = gamma/beta-1; a5 = dt*(gamma/(2*beta)-1);
                effectiveK = K+a0*M+a1*C;
                for k = 2:nt
                    effectiveP = P(:,k) + M*(a0*U(:,k-1)+a2*V(:,k-1)+a3*Acc(:,k-1)) + ...
                        C*(a1*U(:,k-1)+a4*V(:,k-1)+a5*Acc(:,k-1));
                    U(:,k) = effectiveK\effectiveP;
                    Acc(:,k) = a0*(U(:,k)-U(:,k-1))-a2*V(:,k-1)-a3*Acc(:,k-1);
                    V(:,k) = V(:,k-1)+dt*((1-gamma)*Acc(:,k-1)+gamma*Acc(:,k));
                end

                fullU = zeros(nd,nt); fullU(free,:) = U;
                floorResponse=zeros(app.Model.nStories,nt);
                for level=1:app.Model.nStories
                    floorNodes=level*nCols+(1:nCols);
                    floorDOFs=3*floorNodes-2;
                    floorResponse(level,:)=mean(fullU(floorDOFs,:),1);
                end
                roofResponse=floorResponse(end,:);
                floorColors=turbo(app.Model.nStories);
                colororder(app.TimeHistoryAxes,floorColors);
                plot(app.TimeHistoryAxes,time,1000*floorResponse.','LineWidth',1.45);
                grid(app.TimeHistoryAxes,'on'); xlabel(app.TimeHistoryAxes,'Time (s)');
                ylabel(app.TimeHistoryAxes,'Mean floor displacement (mm)');
                title(app.TimeHistoryAxes,'All-floor displacement time histories');
                legend(app.TimeHistoryAxes,compose('Floor %d',1:app.Model.nStories), ...
                    'Location','best','NumColumns',min(2,app.Model.nStories));
                timeCursor=xline(app.TimeHistoryAxes,0,'--','Animation time', ...
                    'Color',[0.12 0.16 0.20],'LineWidth',1.2);

                span = max(max(app.Model.Nodes(:,2))-min(app.Model.Nodes(:,2)), ...
                    max(app.Model.Nodes(:,3))-min(app.Model.Nodes(:,3)));
                maxMotion = max(hypot(fullU(1:3:end,:),fullU(2:3:end,:)),[],'all');
                plotScale = 0;
                if maxMotion>0, plotScale=0.12*span/maxMotion; end
                app.plotStructure(app.DynamicAxes,app.AnalysisResult,[],[],zeros(size(app.Model.Elements,1),1),false);
                hold(app.DynamicAxes,'on');
                edges=app.Model.Elements(:,2:3);
                edgeX=reshape([app.Model.Nodes(edges(:,1),2),app.Model.Nodes(edges(:,2),2),nan(size(edges,1),1)].',[],1);
                edgeY=reshape([app.Model.Nodes(edges(:,1),3),app.Model.Nodes(edges(:,2),3),nan(size(edges,1),1)].',[],1);
                animatedFront=plot(app.DynamicAxes,edgeX,edgeY,'-', ...
                    'Color',[1.00 .25 .08],'LineWidth',2.8,'HandleVisibility','off');
                movingNodes=scatter(app.DynamicAxes,app.Model.Nodes(:,2),app.Model.Nodes(:,3), ...
                    35,'filled','MarkerFaceColor',[1 .65 .10],'HandleVisibility','off');
                xlim(app.DynamicAxes,[min(app.Model.Nodes(:,2))-.18*span,max(app.Model.Nodes(:,2))+.18*span]);
                ylim(app.DynamicAxes,[min(app.Model.Nodes(:,3))-.18*span,max(app.Model.Nodes(:,3))+.18*span]);
                hold(app.DynamicAxes,'off');
                animationFrames = unique(round(linspace(1,nt,min(90,nt))));
                for k = animationFrames
                    if app.StopRequested, break; end
                    ux=fullU(1:3:end,k); uz=fullU(2:3:end,k);
                    deformed=[app.Model.Nodes(:,2)+plotScale*ux,app.Model.Nodes(:,3)+plotScale*uz];
                    edgeX=reshape([deformed(edges(:,1),1),deformed(edges(:,2),1),nan(size(edges,1),1)].',[],1);
                    edgeZ=reshape([deformed(edges(:,1),2),deformed(edges(:,2),2),nan(size(edges,1),1)].',[],1);
                    animatedFront.XData=edgeX; animatedFront.YData=edgeZ;
                    movingNodes.XData=deformed(:,1); movingNodes.YData=deformed(:,2);
                    timeCursor.Value=time(k);
                    title(app.DynamicAxes,sprintf('2D dynamic response | t = %.2f s',time(k)));
                    app.DynamicStatus.Text = sprintf('Animating %d-storey frame | t %.2f / %.2f s | f_1 %.2f Hz', ...
                        app.Model.nStories,time(k),time(end),omega1/(2*pi));
                    drawnow limitrate;
                    pause(0.012);
                end
                if ~app.StopRequested
                    app.DynamicStatus.Text = sprintf('Complete | f_1 = %.2f Hz | peak roof response = %.3f mm', ...
                        omega1/(2*pi),1000*max(abs(roofResponse)));
                end
            catch ME
                app.showError(ME);
            end
        end

        function requireSensors(app)
            if isempty(app.ExistingSensors)
                error('shm:NoSensors','Add at least one existing sensor before continuing.');
            end
        end

        function clearInferenceState(app)
            app.Measurements = zeros(0,1);
            app.MeasurementSigma = zeros(0,1);
            app.Diagnosis = [];
            app.Recommendation = [];
            if ~isempty(app.DiagnosisTable), app.DiagnosisTable.Data = table(); end
            if ~isempty(app.RecommendationTable), app.RecommendationTable.Data = table(); end
            if ~isempty(app.ObservabilityStatus)
                app.ObservabilityStatus.Text='Model or sensors changed. Run CHECK SENSOR ADEQUACY again.';
            end
        end

        function refreshSensorTable(app)
            if isempty(app.ExistingSensors)
                app.SensorTable.Data = table([],[],[],'VariableNames',{'Node','Direction','NoisePercent'});
            else
                app.SensorTable.Data = array2table(app.ExistingSensors,'VariableNames',{'Node','Direction','NoisePercent'});
            end
        end

        function plotStructure(app,ax,result,existingDOFs,proposedDOFs,alpha,showDeformed)
            if nargin<7, showDeformed=false; end
            nodes=app.Model.Nodes; elements=app.Model.Elements;
            if isempty(alpha), alpha=zeros(size(elements,1),1); end
            cla(ax,'reset'); hold(ax,'on'); view(ax,2);
            span=max([max(nodes(:,2))-min(nodes(:,2)),max(nodes(:,3))-min(nodes(:,3)),1]);
            handles=gobjects(0); names={}; hasHealthy=false; hasDamage=false;
            for e=1:size(elements,1)
                xy=nodes(elements(e,2:3),2:3);
                if alpha(e)>0, color=[.85 .18 .12]; else, color=[.20 .27 .34]; end
                h=plot(ax,xy(:,1),xy(:,2),'-','Color',color,'LineWidth',2.8,'HandleVisibility','off');
                if alpha(e)>0 && ~hasDamage
                    h.HandleVisibility='on'; handles(end+1)=h; names{end+1}='Affected member'; hasDamage=true;
                elseif alpha(e)<=0 && ~hasHealthy
                    h.HandleVisibility='on'; handles(end+1)=h; names{end+1}='Frame'; hasHealthy=true;
                end
                mid=mean(xy,1);
                text(ax,mid(1)+.008*span,mid(2)+.008*span,sprintf('e%d',elements(e,1)), ...
                    'FontSize',8,'Color',[.4 .4 .4],'BackgroundColor','white','Margin',1,'Clipping','on');
            end
            scatter(ax,nodes(:,2),nodes(:,3),26,[.2 .25 .3],'filled','HandleVisibility','off');
            for n=1:size(nodes,1)
                if any(ismember(3*n-2:3*n,app.Model.FixedDOFs)), offset=.024*span; else, offset=-.018*span; end
                text(ax,nodes(n,2)-.012*span,nodes(n,3)+offset,sprintf('N%d',nodes(n,1)), ...
                    'HorizontalAlignment','right','VerticalAlignment','top','FontSize',8,'Clipping','on');
            end
            fixedNodes=unique(ceil(app.Model.FixedDOFs(:)/3));
            sdt_drawSupports(ax,app.Model,span);
            bounds=nodes(:,2:3); plotTitle=[app.Model.Name ' | 2D front view'];
            if showDeformed && ~isempty(result)
                ux=result.Displacements(1:3:end); uy=result.Displacements(2:3:end);
                mag=max(hypot(ux,uy));
                if isfield(result,'PlotScale'), scale=result.PlotScale;
                elseif mag>0, scale=.08*span/mag; else, scale=0; end
                h=sdt_drawFrame(ax,app.Model,result.Displacements,scale,[1 .35 .08],'-');
                if ~isempty(h), handles(end+1)=h; names{end+1}='Deformed (scaled)'; end
                bounds=[bounds;nodes(:,2)+scale*ux nodes(:,3)+scale*uy];
                plotTitle=sprintf('%s | 2D | deformation x %.3g',app.Model.Name,scale);
            end
            % All sensors are drawn at analytical nodes in the SAME XY plane.
            % Multiple directions at a node share one marker and a combined label.
            for group=1:2
                if group==1
                    dofs=existingDOFs; color=[0 .42 .85]; marker='o'; label='S'; offset=.024*span;
                else
                    dofs=proposedDOFs; color=[.85 .48 0]; marker='p'; label='New'; offset=-.040*span;
                end
                if isempty(dofs), continue; end
                sensorNodes=ceil(dofs(:)/3); dirs=mod(dofs(:)-1,3)+1;
                uniqueNodes=unique(sensorNodes,'stable');
                h=scatter(ax,nodes(uniqueNodes,2),nodes(uniqueNodes,3),140,marker, ...
                    'MarkerEdgeColor',color,'MarkerFaceColor','none','LineWidth',2);
                handles(end+1)=h;
                if group==1, names{end+1}='Existing sensors'; else, names{end+1}='Recommended sensors'; end
                dirNames={'X','Y','Rz'};
                for n=uniqueNodes.'
                    d=dirs(sensorNodes==n);
                    text(ax,nodes(n,2)+.020*span,nodes(n,3)+offset, ...
                        sprintf('%s N%d: %s',label,n,strjoin(dirNames(d),',')), ...
                        'Color',color,'FontWeight','bold','FontSize',9,'BackgroundColor','white', ...
                        'Margin',1,'Clipping','on');
                end
            end
            axis(ax,'equal'); view(ax,2); ax.Projection='orthographic';
            ax.PlotBoxAspectRatioMode='auto'; ax.Color=[.98 .99 1];
            xlim(ax,[min(bounds(:,1))-.10*span,max(bounds(:,1))+.14*span]);
            ylim(ax,[min(bounds(:,2))-.09*span,max(bounds(:,2))+.09*span]);
            xlabel(ax,'X - horizontal (m)'); ylabel(ax,'Y - height (m)');
            title(ax,{plotTitle,'Fixed column bases: Ux = Uy = Rz = 0'}); grid(ax,'on'); box(ax,'on');
            if ~isempty(handles), legend(ax,handles,names,'Location','best','FontSize',9); end
            hold(ax,'off');
        end

        function showError(app,ME)
            fprintf(2,'APP ERROR: %s\n',getReport(ME,'extended','hyperlinks','off'));
            uialert(app.UIFigure,ME.message,'Analysis error','Icon','error');
        end

    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.bindRuntimeCallbacks();
            app.configureAxesInteractions();
            app.createStabilityTab();
            app.createEnvironmentTab();
            app.generateModel();
            app.runStability();
            app.TabGroup.SelectedTab = app.DesignTab4;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.94 0.96 0.98];
            app.UIFigure.Position = [45 45 1380 840];
            app.UIFigure.Name = 'FrameLab 2D';

            % Create DesignGridLayout1
            app.DesignGridLayout1 = uigridlayout(app.UIFigure);
            app.DesignGridLayout1.ColumnWidth = {'1x'};
            app.DesignGridLayout1.RowHeight = {96, '1x'};
            app.DesignGridLayout1.RowSpacing = 8;

            % Create DesignPanel2
            app.DesignPanel2 = uipanel(app.DesignGridLayout1);
            app.DesignPanel2.BorderType = 'none';
            app.DesignPanel2.BackgroundColor = [0.035 0.12 0.22];
            app.DesignPanel2.Layout.Row = 1;
            app.DesignPanel2.Layout.Column = 1;

            % Create DesignGridLayout3
            app.DesignGridLayout3 = uigridlayout(app.DesignPanel2);
            app.DesignGridLayout3.ColumnWidth = {'1x', 230};
            app.DesignGridLayout3.RowHeight = {40, 32};
            app.DesignGridLayout3.RowSpacing = 4;
            app.DesignGridLayout3.Padding = [18 5 18 4];
            app.DesignGridLayout3.BackgroundColor = [0.035 0.12 0.22];

            % Create DesignLabel10
            app.DesignLabel10 = uilabel(app.DesignGridLayout3);
            app.DesignLabel10.FontSize = 23;
            app.DesignLabel10.FontWeight = 'bold';
            app.DesignLabel10.FontColor = [1 1 1];
            app.DesignLabel10.Layout.Row = 1;
            app.DesignLabel10.Layout.Column = 1;
            app.DesignLabel10.Text = 'FrameLab 2D';

            % Create DesignLabel11
            app.DesignLabel11 = uilabel(app.DesignGridLayout3);
            app.DesignLabel11.FontColor = [0.45 0.86 0.92];
            app.DesignLabel11.Layout.Row = 2;
            app.DesignLabel11.Layout.Column = [1 2];
            app.DesignLabel11.WordWrap = 'on';
            app.DesignLabel11.Text = 'Structural Analysis, Deterioration Assessment and Sensor Placement';

            % Create DesignLabel12
            app.DesignLabel12 = uilabel(app.DesignGridLayout3);
            app.DesignLabel12.HorizontalAlignment = 'center';
            app.DesignLabel12.FontWeight = 'bold';
            app.DesignLabel12.FontColor = [1 1 1];
            app.DesignLabel12.Layout.Row = 1;
            app.DesignLabel12.Layout.Column = 2;
            app.DesignLabel12.Text = 'By Afroja Sultana';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.DesignGridLayout1);
            app.TabGroup.Layout.Row = 2;
            app.TabGroup.Layout.Column = 1;

            % Create DesignTab4
            app.DesignTab4 = uitab(app.TabGroup);
            app.DesignTab4.Title = '1. 2D Frame Model';
            app.DesignTab4.BackgroundColor = [0.96 0.98 1];

            % Create DesignGridLayout13
            app.DesignGridLayout13 = uigridlayout(app.DesignTab4);
            app.DesignGridLayout13.ColumnWidth = {380, '1x'};
            app.DesignGridLayout13.RowHeight = {190, 46, '1x'};
            app.DesignGridLayout13.Padding = [12 12 12 12];

            % Create ModelAxes
            app.ModelAxes = uiaxes(app.DesignGridLayout13);
            xlabel(app.ModelAxes, 'Width (m)')
            ylabel(app.ModelAxes, 'Height (m)')
            zlabel(app.ModelAxes, 'Height (m)')
            app.ModelAxes.CameraViewAngle = 9.00116867899543;
            app.ModelAxes.View = [0 90];
            app.ModelAxes.Projection = 'orthographic';
            app.ModelAxes.DataAspectRatioMode = 'manual';
            app.ModelAxes.PlotBoxAspectRatio = [1.4 1 1.15];
            app.ModelAxes.XLim = [-1 13];
            app.ModelAxes.YLim = [-1 9];
            app.ModelAxes.ZLim = [-0.2 11.3];
            app.ModelAxes.Color = [0.92 0.96 0.99];
            app.ModelAxes.Layout.Row = 3;
            app.ModelAxes.Layout.Column = 2;

            % Create DesignGridLayout19
            app.DesignGridLayout19 = uigridlayout(app.DesignGridLayout13);
            app.DesignGridLayout19.ColumnWidth = {110, 110, 110, 110};
            app.DesignGridLayout19.RowHeight = {28,28,28,28,24};
            app.DesignGridLayout19.RowSpacing=6; app.DesignGridLayout19.Padding=[8 8 8 8];
            app.DesignGridLayout19.Layout.Row = 1;
            app.DesignGridLayout19.Layout.Column = [1 2];

            % Create DesignLabel26
            app.DesignLabel26 = uilabel(app.DesignGridLayout19);
            app.DesignLabel26.Layout.Row = 1;
            app.DesignLabel26.Layout.Column = 1;
            app.DesignLabel26.Text = 'Stories';

            % Create StoriesField
            app.StoriesField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.StoriesField.Limits = [1 12];
            app.StoriesField.RoundFractionalValues = 'on';
            app.StoriesField.Layout.Row = 1;
            app.StoriesField.Layout.Column = 2;
            app.StoriesField.Value = 3;

            % Create DesignLabel27
            app.DesignLabel27 = uilabel(app.DesignGridLayout19);
            app.DesignLabel27.Layout.Row = 1;
            app.DesignLabel27.Layout.Column = 3;
            app.DesignLabel27.Text = 'Bays';

            % Create BaysField
            app.BaysField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.BaysField.Limits = [1 12];
            app.BaysField.RoundFractionalValues = 'on';
            app.BaysField.Layout.Row = 1;
            app.BaysField.Layout.Column = 4;
            app.BaysField.Value = 2;

            % Create DesignLabel28
            app.DesignLabel28 = uilabel(app.DesignGridLayout19);
            app.DesignLabel28.Layout.Row = 2;
            app.DesignLabel28.Layout.Column = 1;
            app.DesignLabel28.Text = 'Bay width (m)';

            % Create BayWidthField
            app.BayWidthField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.BayWidthField.Limits = [0.1 Inf];
            app.BayWidthField.Layout.Row = 2;
            app.BayWidthField.Layout.Column = 2;
            app.BayWidthField.Value = 6;

            % Create DesignLabel29
            app.DesignLabel29 = uilabel(app.DesignGridLayout19);
            app.DesignLabel29.Layout.Row = 2;
            app.DesignLabel29.Layout.Column = 3;
            app.DesignLabel29.Text = 'Story height (m)';

            % Create StoryHeightField
            app.StoryHeightField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.StoryHeightField.Limits = [0.1 Inf];
            app.StoryHeightField.Layout.Row = 2;
            app.StoryHeightField.Layout.Column = 4;
            app.StoryHeightField.Value = 3.5;

            % Create DesignLabel30
            app.DesignLabel30 = uilabel(app.DesignGridLayout19);
            app.DesignLabel30.Layout.Row = 3;
            app.DesignLabel30.Layout.Column = 1;
            app.DesignLabel30.Text = 'E (GPa)';

            % Create EField
            app.EField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.EField.Limits = [1 Inf];
            app.EField.Layout.Row = 3;
            app.EField.Layout.Column = 2;
            app.EField.Value = 200;

            % Create DesignLabel31
            app.DesignLabel31 = uilabel(app.DesignGridLayout19);
            app.DesignLabel31.Layout.Row = 3;
            app.DesignLabel31.Layout.Column = 3;
            app.DesignLabel31.Text = 'A (m^2)';

            % Create AField
            app.AField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.AField.Limits = [1e-06 Inf];
            app.AField.Layout.Row = 3;
            app.AField.Layout.Column = 4;
            app.AField.Value = 0.02;

            % Create DesignLabel32
            app.DesignLabel32 = uilabel(app.DesignGridLayout19);
            app.DesignLabel32.Layout.Row = 4;
            app.DesignLabel32.Layout.Column = 1;
            app.DesignLabel32.Text = 'I (m^4)';

            % Create IField
            app.IField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.IField.Limits = [1e-10 Inf];
            app.IField.Layout.Row = 4;
            app.IField.Layout.Column = 2;
            app.IField.Value = 8e-05;

            % Create DesignLabel33
            app.DesignLabel33 = uilabel(app.DesignGridLayout19);
            app.DesignLabel33.Layout.Row = 4;
            app.DesignLabel33.Layout.Column = 3;
            app.DesignLabel33.Text = 'Lateral load (kN)';

            % Create LoadField
            app.LoadField = uieditfield(app.DesignGridLayout19, 'numeric');
            app.LoadField.Limits = [0 Inf];
            app.LoadField.Layout.Row = 4;
            app.LoadField.Layout.Column = 4;
            app.LoadField.Value = 100;

            % Create DesignLabel34
            app.DesignLabel34 = uilabel(app.DesignGridLayout19);
            app.DesignLabel34.Layout.Row = 5;
            app.DesignLabel34.Layout.Column = 1;
            app.DesignLabel34.Text = '2D frame: X, Y translation and Z rotation';
            app.DesignLabel34.Layout.Column=[1 4];

            % Create DesignGridLayout20
            app.DesignGridLayout20 = uigridlayout(app.DesignGridLayout13);
            app.DesignGridLayout20.ColumnWidth = {'1x', '1x', '1x'};
            app.DesignGridLayout20.RowHeight = {'1x'};
            app.DesignGridLayout20.ColumnSpacing = 12;
            app.DesignGridLayout20.Padding = [0 4 0 4];
            app.DesignGridLayout20.Layout.Row = 2;
            app.DesignGridLayout20.Layout.Column = [1 2];

            % Create DesignButton35
            app.DesignButton35 = uibutton(app.DesignGridLayout20, 'push');
            app.DesignButton35.BackgroundColor = [0.04 0.58 0.44];
            app.DesignButton35.FontSize = 15;
            app.DesignButton35.FontWeight = 'bold';
            app.DesignButton35.FontColor = [1 1 1];
            app.DesignButton35.Tooltip = 'Apply all edited inputs and rebuild the building';
            app.DesignButton35.Layout.Row = 1;
            app.DesignButton35.Layout.Column = 1;
            app.DesignButton35.Text = 'INSERT';

            % Create DesignButton36
            app.DesignButton36 = uibutton(app.DesignGridLayout20, 'push');
            app.DesignButton36.BackgroundColor = [0.4 0.28 0.68];
            app.DesignButton36.FontWeight = 'bold';
            app.DesignButton36.FontColor = [1 1 1];
            app.DesignButton36.Tooltip = 'Enable mouse drag rotation of the building';
            app.DesignButton36.Layout.Row = 1;
            app.DesignButton36.Layout.Column = 2;
            app.DesignButton36.Text = 'RESET 2D VIEW';

            % Create DesignButton37
            app.DesignButton37 = uibutton(app.DesignGridLayout20, 'push');
            app.DesignButton37.BackgroundColor = [0.1 0.36 0.65];
            app.DesignButton37.FontWeight = 'bold';
            app.DesignButton37.FontColor = [1 1 1];
            app.DesignButton37.Layout.Row = 1;
            app.DesignButton37.Layout.Column = 3;
            app.DesignButton37.Text = 'RUN STATIC ANALYSIS';

            % Create ModelStatus
            app.ModelStatus = uilabel(app.DesignGridLayout13);
            app.ModelStatus.VerticalAlignment = 'top';
            app.ModelStatus.FontWeight = 'bold';
            app.ModelStatus.Layout.Row = 3;
            app.ModelStatus.Layout.Column = 1;
            app.ModelStatus.Text = 'INPUT ACCEPTED: 3-story, 2-bay frame | 12 nodes | 15 elements';

            % Create DesignTab5
            app.DesignTab5 = uitab(app.TabGroup);
            app.DesignTab5.Title = '2. Sensor Network';
            app.DesignTab5.BackgroundColor = [0.96 0.98 1];

            % Create DesignGridLayout14
            app.DesignGridLayout14 = uigridlayout(app.DesignTab5);
            app.DesignGridLayout14.ColumnWidth = {420, '1x'};
            app.DesignGridLayout14.RowHeight = {90, '1x', 100};
            app.DesignGridLayout14.Padding = [12 12 12 12];

            % Create SensorAxes
            app.SensorAxes = uiaxes(app.DesignGridLayout14);
            title(app.SensorAxes, '3-story, 2-bay frame | 2D frame front view')
            xlabel(app.SensorAxes, 'Width (m)')
            ylabel(app.SensorAxes, 'Height (m)')
            zlabel(app.SensorAxes, 'Height (m)')
            app.SensorAxes.CameraViewAngle = 9.00116867899543;
            app.SensorAxes.View = [0 90];
            app.SensorAxes.Projection = 'orthographic';
            app.SensorAxes.DataAspectRatioMode = 'manual';
            app.SensorAxes.PlotBoxAspectRatio = [1.4 1 1.15];
            app.SensorAxes.XLim = [-1 13];
            app.SensorAxes.YLim = [-1 9];
            app.SensorAxes.ZLim = [-0.2 11.3];
            app.SensorAxes.Color = [0.92 0.96 0.99];
            app.SensorAxes.Layout.Row = 2;
            app.SensorAxes.Layout.Column = 2;

            % Create DesignGridLayout21
            app.DesignGridLayout21 = uigridlayout(app.DesignGridLayout14);
            app.DesignGridLayout21.ColumnWidth = {45, 80, 70, 100, 65, 90};
            app.DesignGridLayout21.RowHeight = {'1x', '1x', '1x'};
            app.DesignGridLayout21.Layout.Row = 1;
            app.DesignGridLayout21.Layout.Column = [1 2];

            % Create DesignLabel38
            app.DesignLabel38 = uilabel(app.DesignGridLayout21);
            app.DesignLabel38.Layout.Row = 1;
            app.DesignLabel38.Layout.Column = 1;
            app.DesignLabel38.Text = 'Node';

            % Create SensorNodeField
            app.SensorNodeField = uieditfield(app.DesignGridLayout21, 'numeric');
            app.SensorNodeField.Limits = [1 12];
            app.SensorNodeField.RoundFractionalValues = 'on';
            app.SensorNodeField.Layout.Row = 1;
            app.SensorNodeField.Layout.Column = 2;
            app.SensorNodeField.Value = 12;

            % Create DesignLabel39
            app.DesignLabel39 = uilabel(app.DesignGridLayout21);
            app.DesignLabel39.Layout.Row = 1;
            app.DesignLabel39.Layout.Column = 3;
            app.DesignLabel39.Text = 'Direction';

            % Create SensorDirectionDropDown
            app.SensorDirectionDropDown = uidropdown(app.DesignGridLayout21);
            app.SensorDirectionDropDown.Items = {'X displacement', 'Y displacement', 'Rotation'};
            app.SensorDirectionDropDown.Layout.Row = 1;
            app.SensorDirectionDropDown.Layout.Column = 4;
            app.SensorDirectionDropDown.Value = 'X displacement';

            % Create DesignLabel40
            app.DesignLabel40 = uilabel(app.DesignGridLayout21);
            app.DesignLabel40.Layout.Row = 1;
            app.DesignLabel40.Layout.Column = 5;
            app.DesignLabel40.Text = 'Noise %';

            % Create SensorNoiseField
            app.SensorNoiseField = uieditfield(app.DesignGridLayout21, 'numeric');
            app.SensorNoiseField.Limits = [0 25];
            app.SensorNoiseField.Layout.Row = 1;
            app.SensorNoiseField.Layout.Column = 6;
            app.SensorNoiseField.Value = 1;

            % Create DesignButton41
            app.DesignButton41 = uibutton(app.DesignGridLayout21, 'push');
            app.DesignButton41.BackgroundColor = [0.04 0.58 0.44];
            app.DesignButton41.FontWeight = 'bold';
            app.DesignButton41.FontColor = [1 1 1];
            app.DesignButton41.Layout.Row = [2 3];
            app.DesignButton41.Layout.Column = [1 3];
            app.DesignButton41.Text = 'ADD SENSOR';

            % Create DesignButton42
            app.DesignButton42 = uibutton(app.DesignGridLayout21, 'push');
            app.DesignButton42.BackgroundColor = [0.82 0.27 0.25];
            app.DesignButton42.FontWeight = 'bold';
            app.DesignButton42.FontColor = [1 1 1];
            app.DesignButton42.Layout.Row = [2 3];
            app.DesignButton42.Layout.Column = [4 6];
            app.DesignButton42.Text = 'REMOVE SELECTED';

            % Create SensorTable
            app.SensorTable = uitable(app.DesignGridLayout14);
            app.SensorTable.ColumnEditable = [true true true];
            app.SensorTable.Layout.Row = 2;
            app.SensorTable.Layout.Column = 1;

            % Create DesignTab6
            app.DesignTab6 = uitab(app.TabGroup);
            app.DesignTab6.Title = '3. Damage Diagnosis';
            app.DesignTab6.BackgroundColor = [0.96 0.98 1];

            % Create DesignGridLayout15
            app.DesignGridLayout15 = uigridlayout(app.DesignTab6);
            app.DesignGridLayout15.ColumnWidth = {'1x', 410};
            app.DesignGridLayout15.RowHeight = {90, '1x', 35};
            app.DesignGridLayout15.Padding = [12 12 12 12];

            % Separate spatial diagnosis from the single-member trial chart.
            diagnosisPlots=uitabgroup(app.DesignGridLayout15);
            diagnosisPlots.Layout.Row=2; diagnosisPlots.Layout.Column=1;
            frameTab=uitab(diagnosisPlots,'Title','Frame and supports');
            frameGrid=uigridlayout(frameTab,[1 1]); frameGrid.Padding=[6 6 6 6];
            app.DiagnosisAxes=uiaxes(frameGrid); app.DiagnosisAxes.Tag='DiagnosisFrame';
            scanTab=uitab(diagnosisPlots,'Title','Single-member trial chart');
            scanGrid=uigridlayout(scanTab,[1 1]); scanGrid.Padding=[6 6 6 6];
            app.DiagnosisScanAxes=uiaxes(scanGrid); app.DiagnosisScanAxes.Tag='DiagnosisTrials';

            % Create DesignGridLayout22
            app.DesignGridLayout22 = uigridlayout(app.DesignGridLayout15);
            app.DesignGridLayout22.ColumnWidth = {90, 90, 120, '1x', 100, 120, 120};
            app.DesignGridLayout22.Layout.Row = 1;
            app.DesignGridLayout22.Layout.Column = [1 2];

            % Create DesignLabel43
            app.DesignLabel43 = uilabel(app.DesignGridLayout22);
            app.DesignLabel43.Layout.Row = 1;
            app.DesignLabel43.Layout.Column = 1;
            app.DesignLabel43.Text = 'Hidden element';

            % Create HiddenElementField
            app.HiddenElementField = uieditfield(app.DesignGridLayout22, 'numeric');
            app.HiddenElementField.Limits = [1 15];
            app.HiddenElementField.RoundFractionalValues = 'on';
            app.HiddenElementField.Layout.Row = 1;
            app.HiddenElementField.Layout.Column = 2;
            app.HiddenElementField.Value = 1;

            % Create DesignLabel44
            app.DesignLabel44 = uilabel(app.DesignGridLayout22);
            app.DesignLabel44.Layout.Row = 1;
            app.DesignLabel44.Layout.Column = 3;
            app.DesignLabel44.Text = 'Damage severity';

            % Create DamageSlider
            app.DamageSlider = uislider(app.DesignGridLayout22);
            app.DamageSlider.Limits = [5 70];
            app.DamageSlider.MajorTicks=5:5:70; app.DamageSlider.MinorTicks=[];
            app.DamageSlider.FontSize = 8;
            app.DamageSlider.Layout.Row = 1;
            app.DamageSlider.Layout.Column = 4;
            app.DamageSlider.Value = 30;

            % Create DamageValueLabel
            app.DamageValueLabel = uilabel(app.DesignGridLayout22);
            app.DamageValueLabel.Layout.Row = 1;
            app.DamageValueLabel.Layout.Column = 5;
            app.DamageValueLabel.Text = '30% stiffness loss';

            % Create DesignLabel45
            app.DesignLabel45 = uilabel(app.DesignGridLayout22);
            app.DesignLabel45.Layout.Row = 2;
            app.DesignLabel45.Layout.Column = 1;
            app.DesignLabel45.Text = 'Noise %';

            % Create SimulationNoiseField
            app.SimulationNoiseField = uieditfield(app.DesignGridLayout22, 'numeric');
            app.SimulationNoiseField.Limits = [0 25];
            app.SimulationNoiseField.Layout.Row = 2;
            app.SimulationNoiseField.Layout.Column = 2;
            app.SimulationNoiseField.Value = 1;

            % Create DesignButton46
            app.DesignButton46 = uibutton(app.DesignGridLayout22, 'push');
            app.DesignButton46.BackgroundColor = [0.55 0.3 0.72];
            app.DesignButton46.FontWeight = 'bold';
            app.DesignButton46.FontColor = [1 1 1];
            app.DesignButton46.Layout.Row = 2;
            app.DesignButton46.Layout.Column = [5 6];
            app.DesignButton46.Text = 'SIMULATE DATA';

            % Create DesignButton47
            app.DesignButton47 = uibutton(app.DesignGridLayout22, 'push');
            app.DesignButton47.BackgroundColor = [0.9 0.45 0.1];
            app.DesignButton47.FontWeight = 'bold';
            app.DesignButton47.FontColor = [1 1 1];
            app.DesignButton47.Layout.Row = 2;
            app.DesignButton47.Layout.Column = 7;
            app.DesignButton47.Text = 'RUN DIAGNOSIS';

            % Create DiagnosisTable
            app.DiagnosisTable = uitable(app.DesignGridLayout15);
            app.DiagnosisTable.Layout.Row = 2;
            app.DiagnosisTable.Layout.Column = 2;

            % Create DiagnosisStatus
            app.DiagnosisStatus = uilabel(app.DesignGridLayout15);
            app.DiagnosisStatus.FontWeight = 'bold';
            app.DiagnosisStatus.Layout.Row = 3;
            app.DiagnosisStatus.Layout.Column = [1 2];
            app.DiagnosisStatus.Text = 'Simulate measurements to begin.';

            % Create DesignTab7
            app.DesignTab7 = uitab(app.TabGroup);
            app.DesignTab7.Title = '4. Dynamic Simulation';
            app.DesignTab7.BackgroundColor = [0.96 0.98 1];

            % Create DesignGridLayout16
            app.DesignGridLayout16 = uigridlayout(app.DesignTab7);
            app.DesignGridLayout16.RowHeight = {105, '1x', 38};
            app.DesignGridLayout16.Padding = [12 12 12 12];

            % Create DynamicAxes
            app.DynamicAxes = uiaxes(app.DesignGridLayout16);
            app.DynamicAxes.Layout.Row = 2;
            app.DynamicAxes.Layout.Column = 1;

            % Create TimeHistoryAxes
            app.TimeHistoryAxes = uiaxes(app.DesignGridLayout16);
            app.TimeHistoryAxes.Layout.Row = 2;
            app.TimeHistoryAxes.Layout.Column = 2;

            % Create DesignGridLayout23
            app.DesignGridLayout23 = uigridlayout(app.DesignGridLayout16);
            app.DesignGridLayout23.ColumnWidth = {125, 110, 125, 110, 110, '1x'};
            app.DesignGridLayout23.Layout.Row = 1;
            app.DesignGridLayout23.Layout.Column = [1 2];

            % Create DesignLabel48
            app.DesignLabel48 = uilabel(app.DesignGridLayout23);
            app.DesignLabel48.FontWeight = 'bold';
            app.DesignLabel48.Layout.Row = 1;
            app.DesignLabel48.Layout.Column = 1;
            app.DesignLabel48.Text = 'Floor mass (tonne)';

            % Create DesignLabel49
            app.DesignLabel49 = uilabel(app.DesignGridLayout23);
            app.DesignLabel49.FontWeight = 'bold';
            app.DesignLabel49.Layout.Row = 1;
            app.DesignLabel49.Layout.Column = 2;
            app.DesignLabel49.Text = 'Damping (%)';

            % Create DesignLabel50
            app.DesignLabel50 = uilabel(app.DesignGridLayout23);
            app.DesignLabel50.FontWeight = 'bold';
            app.DesignLabel50.Layout.Row = 1;
            app.DesignLabel50.Layout.Column = 3;
            app.DesignLabel50.Text = 'Excitation (Hz)';

            % Create DesignLabel51
            app.DesignLabel51 = uilabel(app.DesignGridLayout23);
            app.DesignLabel51.FontWeight = 'bold';
            app.DesignLabel51.Layout.Row = 1;
            app.DesignLabel51.Layout.Column = 4;
            app.DesignLabel51.Text = 'Duration (s)';

            % Create DesignLabel52
            app.DesignLabel52 = uilabel(app.DesignGridLayout23);
            app.DesignLabel52.FontWeight = 'bold';
            app.DesignLabel52.Layout.Row = 1;
            app.DesignLabel52.Layout.Column = 5;
            app.DesignLabel52.Text = 'Time step (s)';

            % Create FloorMassField
            app.FloorMassField = uieditfield(app.DesignGridLayout23, 'numeric');
            app.FloorMassField.Limits = [0.1 Inf];
            app.FloorMassField.Layout.Row = 2;
            app.FloorMassField.Layout.Column = 1;
            app.FloorMassField.Value = 80;

            % Create DampingField
            app.DampingField = uieditfield(app.DesignGridLayout23, 'numeric');
            app.DampingField.Limits = [0 30];
            app.DampingField.Layout.Row = 2;
            app.DampingField.Layout.Column = 2;
            app.DampingField.Value = 5;

            % Create ExcitationFrequencyField
            app.ExcitationFrequencyField = uieditfield(app.DesignGridLayout23, 'numeric');
            app.ExcitationFrequencyField.Limits = [0.01 20];
            app.ExcitationFrequencyField.Layout.Row = 2;
            app.ExcitationFrequencyField.Layout.Column = 3;
            app.ExcitationFrequencyField.Value = 1.2;

            % Create DurationField
            app.DurationField = uieditfield(app.DesignGridLayout23, 'numeric');
            app.DurationField.Limits = [0.2 60];
            app.DurationField.Layout.Row = 2;
            app.DurationField.Layout.Column = 4;
            app.DurationField.Value = 6;

            % Create TimeStepField
            app.TimeStepField = uieditfield(app.DesignGridLayout23, 'numeric');
            app.TimeStepField.Limits = [0.001 0.2];
            app.TimeStepField.Layout.Row = 2;
            app.TimeStepField.Layout.Column = 5;
            app.TimeStepField.Value = 0.02;

            % Create DesignGridLayout53
            app.DesignGridLayout53 = uigridlayout(app.DesignGridLayout23);
            app.DesignGridLayout53.ColumnWidth = {'1x'};
            app.DesignGridLayout53.Padding = [4 0 4 0];
            app.DesignGridLayout53.Layout.Row = [1 2];
            app.DesignGridLayout53.Layout.Column = 6;

            % Create DesignButton56
            app.DesignButton56 = uibutton(app.DesignGridLayout53, 'push');
            app.DesignButton56.BackgroundColor = [0.04 0.48 0.58];
            app.DesignButton56.FontWeight = 'bold';
            app.DesignButton56.FontColor = [1 1 1];
            app.DesignButton56.Layout.Row = 1;
            app.DesignButton56.Layout.Column = 1;
            app.DesignButton56.Text = 'RUN ANIMATION';

            % Create DesignButton57
            app.DesignButton57 = uibutton(app.DesignGridLayout53, 'push');
            app.DesignButton57.BackgroundColor = [0.82 0.27 0.25];
            app.DesignButton57.FontWeight = 'bold';
            app.DesignButton57.FontColor = [1 1 1];
            app.DesignButton57.Layout.Row = 2;
            app.DesignButton57.Layout.Column = 1;
            app.DesignButton57.Text = 'STOP';

            % Create DynamicStatus
            app.DynamicStatus = uilabel(app.DesignGridLayout16);
            app.DynamicStatus.FontWeight = 'bold';
            app.DynamicStatus.Layout.Row = 3;
            app.DynamicStatus.Layout.Column = [1 2];
            app.DynamicStatus.Text = 'Set the dynamic parameters and run the Newmark-beta time-history simulation.';

            % Compact sensor adequacy check inside Sensor Network.
            adequacy=uigridlayout(app.DesignGridLayout14,[1 2]);
            adequacy.Tag='SensorAdequacyPanel';
            adequacy.Layout.Row=3; adequacy.Layout.Column=[1 2];
            adequacy.ColumnWidth={190,'1x'};
            adequacy.RowHeight={'1x'}; adequacy.Padding=[6 6 6 6];
            app.DesignButton24=uibutton(adequacy,'push');
            app.DesignButton24.BackgroundColor=[0.1 0.36 0.65];
            app.DesignButton24.FontWeight='bold'; app.DesignButton24.FontColor=[1 1 1];
            app.DesignButton24.Layout.Row=1; app.DesignButton24.Layout.Column=1;
            app.DesignButton24.Text='CHECK SENSOR ADEQUACY';
            app.ObservabilityStatus=uilabel(adequacy);
            app.ObservabilityStatus.Layout.Row=1; app.ObservabilityStatus.Layout.Column=2;
            app.ObservabilityStatus.WordWrap='on';
            app.ObservabilityStatus.Text='Check whether the current measurements can distinguish simultaneous member stiffness changes.';

            % Create DesignTab9
            app.DesignTab9 = uitab(app.TabGroup);
            app.DesignTab9.Title = '5. Sensor Recommendation';
            app.DesignTab9.BackgroundColor = [0.96 0.98 1];

            % Create DesignGridLayout18
            app.DesignGridLayout18 = uigridlayout(app.DesignTab9);
            app.DesignGridLayout18.ColumnWidth = {'1x', 420};
            app.DesignGridLayout18.RowHeight = {70, '1x', 40};
            app.DesignGridLayout18.Padding = [12 12 12 12];

            % Create RecommendationAxes
            app.RecommendationAxes = uiaxes(app.DesignGridLayout18);
            app.RecommendationAxes.Layout.Row = 2;
            app.RecommendationAxes.Layout.Column = 1;

            % Create DesignGridLayout25
            app.DesignGridLayout25 = uigridlayout(app.DesignGridLayout18);
            app.DesignGridLayout25.ColumnWidth = {150, 80, 180, '1x'};
            app.DesignGridLayout25.RowHeight = {'1x'};
            app.DesignGridLayout25.Layout.Row = 1;
            app.DesignGridLayout25.Layout.Column = [1 2];

            % Create DesignLabel54
            app.DesignLabel54 = uilabel(app.DesignGridLayout25);
            app.DesignLabel54.Layout.Row = 1;
            app.DesignLabel54.Layout.Column = 1;
            app.DesignLabel54.Text = 'Number of sensors to add';

            % Create NumberSensorsField
            app.NumberSensorsField = uieditfield(app.DesignGridLayout25, 'numeric');
            app.NumberSensorsField.Limits = [1 5];
            app.NumberSensorsField.RoundFractionalValues = 'on';
            app.NumberSensorsField.Layout.Row = 1;
            app.NumberSensorsField.Layout.Column = 2;
            app.NumberSensorsField.Value = 1;

            % Create DesignButton55
            app.DesignButton55 = uibutton(app.DesignGridLayout25, 'push');
            app.DesignButton55.BackgroundColor = [0 0 1];
            app.DesignButton55.FontSize = 10;
            app.DesignButton55.FontWeight = 'bold';
            app.DesignButton55.FontColor = [1 1 1];
            app.DesignButton55.Layout.Row = 1;
            app.DesignButton55.Layout.Column = 3;
            app.DesignButton55.Text = 'OPTIMIZE SENSOR LOCATIONS';

            % Create RecommendationTable
            app.RecommendationTable = uitable(app.DesignGridLayout18);
            app.RecommendationTable.Layout.Row = 2;
            app.RecommendationTable.Layout.Column = 2;

            % Create RecommendationStatus
            app.RecommendationStatus = uilabel(app.DesignGridLayout18);
            app.RecommendationStatus.FontWeight = 'bold';
            app.RecommendationStatus.Layout.Row = 3;
            app.RecommendationStatus.Layout.Column = [1 2];
            app.RecommendationStatus.Text = 'Run diagnosis before requesting recommendations.';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    methods (Static)
        function report = verifyNewtonRaphson()
            report=sdt_verifyNewtonRaphson();
        end

        function report = verifyEnvironment()
            % Run in MATLAB: FrameLab2D.verifyEnvironment
            report=sdt_verifyEnvironment();
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = FrameLab2D

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

% Embedded helper: computeSensitivity.m
function [H,y0] = sdt_computeSensitivity(model, sensorDOFs, alpha0, step)
%COMPUTESENSITIVITY Forward-difference sensor sensitivity to element damage.
if nargin < 3 || isempty(alpha0)
    alpha0 = zeros(size(model.Elements,1),1);
end
if nargin < 4 || isempty(step)
    step = 1e-3;
end
base = sdt_solveFrame(model,alpha0);
y0 = base.Displacements(sensorDOFs);
nParam = size(model.Elements,1);
H = zeros(numel(sensorDOFs),nParam);
for j = 1:nParam
    trial = alpha0;
    trial(j) = min(alpha0(j)+step,0.95);
    actualStep = trial(j)-alpha0(j);
    r = sdt_solveFrame(model,trial);
    H(:,j) = (r.Displacements(sensorDOFs)-y0)/actualStep;
end
end


% Embedded helper: estimateSingleDamage.m
function diagnosis = sdt_estimateSingleDamage(model, sensorDOFs, measurements, sigma, alphaGrid)
%ESTIMATESINGLEDAMAGE Enumerate one damaged element and its severity.
if nargin < 5 || isempty(alphaGrid)
    alphaGrid = linspace(0.02,0.80,80);
end
measurements = measurements(:);
if isscalar(sigma)
    sigma = repmat(sigma,size(measurements));
end
sigma = max(sigma(:),eps);
nElements = size(model.Elements,1);
candidates = zeros(nElements,4); % element, alpha, objective, RMSE

for e = 1:nElements
    best = [e,0,Inf,Inf];
    for a = alphaGrid
        alpha = zeros(nElements,1); alpha(e)=a;
        r = sdt_solveFrame(model,alpha);
        pred = r.Displacements(sensorDOFs);
        residual = (measurements-pred)./sigma;
        objective = residual.'*residual;
        rmse = sqrt(mean((measurements-pred).^2));
        if objective < best(3)
            best = [e,a,objective,rmse];
        end
    end
    candidates(e,:) = best;
end
[~,order] = sort(candidates(:,3),'ascend');
candidates = candidates(order,:);
diagnosis = struct();
diagnosis.BestElement = candidates(1,1);
diagnosis.BestAlpha = candidates(1,2);
diagnosis.Objective = candidates(1,3);
diagnosis.RMSE = candidates(1,4);
diagnosis.Candidates = candidates;
end


% Embedded helper: frameElement2D.m
function [kGlobal, kLocal, T, dofs] = sdt_frameElement2D(model, elementRow, retainedFactor)
%FRAMEELEMENT2D Return a 2D Euler-Bernoulli frame element stiffness.

if nargin < 3
    retainedFactor = 1.0;
end

ni = elementRow(2);
nj = elementRow(3);
E = elementRow(4);
A = elementRow(5);
I = elementRow(6);

xi = model.Nodes(ni,2); yi = model.Nodes(ni,3);
xj = model.Nodes(nj,2); yj = model.Nodes(nj,3);
dx = xj-xi; dy = yj-yi;
L = hypot(dx,dy);
if L <= eps
    error('shm:ZeroLengthElement','Element %d has zero length.',elementRow(1));
end
c = dx/L; s = dy/L;

kLocal = [ ...
    E*A/L,          0,           0, -E*A/L,           0,           0;
    0,       12*E*I/L^3,  6*E*I/L^2,      0, -12*E*I/L^3,  6*E*I/L^2;
    0,        6*E*I/L^2,    4*E*I/L,      0,  -6*E*I/L^2,    2*E*I/L;
   -E*A/L,          0,           0,  E*A/L,           0,           0;
    0,      -12*E*I/L^3, -6*E*I/L^2,      0,  12*E*I/L^3, -6*E*I/L^2;
    0,        6*E*I/L^2,    2*E*I/L,      0,  -6*E*I/L^2,    4*E*I/L];

% T maps global element displacements to local displacements.
T = [c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0; ...
     0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1];
kLocal = retainedFactor*kLocal;
kGlobal = T.'*kLocal*T;
dofs = [3*ni-2, 3*ni-1, 3*ni, 3*nj-2, 3*nj-1, 3*nj];
end


% Embedded helper: generatePortalFrame.m
function model = sdt_generatePortalFrame(nStories, nBays, bayWidth, storyHeight, E, A, I)
%GENERATEPORTALFRAME Create a regular 2D moment-frame benchmark model.
% Units are SI: m, N, Pa, m^2, and m^4.

arguments
    nStories (1,1) double {mustBeInteger,mustBePositive}
    nBays (1,1) double {mustBeInteger,mustBePositive}
    bayWidth (1,1) double {mustBePositive}
    storyHeight (1,1) double {mustBePositive}
    E (1,1) double {mustBePositive}
    A (1,1) double {mustBePositive}
    I (1,1) double {mustBePositive}
end

nCols = nBays + 1;
nLevels = nStories + 1;
nNodes = nCols * nLevels;
nodes = zeros(nNodes, 3);

nodeId = 0;
for level = 0:nStories
    for col = 0:nBays
        nodeId = nodeId + 1;
        nodes(nodeId,:) = [nodeId, col*bayWidth, level*storyHeight];
    end
end

% Elements: [id, nodeI, nodeJ, E, A, I, kind], kind 1=column, 2=beam.
elements = zeros(nStories*nCols + nStories*nBays, 7);
elementId = 0;
for story = 1:nStories
    lower = (story-1)*nCols;
    upper = story*nCols;
    for col = 1:nCols
        elementId = elementId + 1;
        elements(elementId,:) = [elementId, lower+col, upper+col, E, A, I, 1];
    end
    for bay = 1:nBays
        elementId = elementId + 1;
        elements(elementId,:) = [elementId, upper+bay, upper+bay+1, E, A, I, 2];
    end
end

fixedDOFs = reshape([3*(1:nCols)-2; 3*(1:nCols)-1; 3*(1:nCols)], 1, []);
loads = zeros(3*nNodes,1);
topNodes = nStories*nCols + (1:nCols);
totalLateralLoad = 100e3;
loads(3*topNodes-2) = totalLateralLoad/numel(topNodes);

model = struct();
model.Nodes = nodes;
model.Elements = elements;
model.FixedDOFs = sort(fixedDOFs);
model.Loads = loads;
model.nStories = nStories;
model.nBays = nBays;
model.BayWidth = bayWidth;
model.StoryHeight = storyHeight;
model.Units = "SI";
model.Name = sprintf('%d-story, %d-bay frame', nStories, nBays);
end


% Embedded helper: recommendSensors.m
function recommendation = sdt_recommendSensors(model, existingDOFs, diagnosis, numberToAdd, sigma)
%RECOMMENDSENSORS Greedy scenario-separation placement for x-displacement sensors.
arguments
    model struct
    existingDOFs (:,1) double
    diagnosis struct
    numberToAdd (1,1) double {mustBeInteger,mustBePositive}
    sigma (1,1) double {mustBePositive}
end

% Candidate channels are horizontal DOFs at non-base nodes.
baseNodes = unique(ceil(model.FixedDOFs(:)/3));
allNodes = model.Nodes(:,1);
candidateNodes = setdiff(allNodes,baseNodes);
candidateDOFs = 3*candidateNodes-2;
candidateDOFs = setdiff(candidateDOFs,existingDOFs,'stable');

nScenarios = min(5,size(diagnosis.Candidates,1));
scenarioRows = diagnosis.Candidates(1:nScenarios,:);
predictions = zeros(numel(candidateDOFs),nScenarios);
for k = 1:nScenarios
    alpha = zeros(size(model.Elements,1),1);
    alpha(scenarioRows(k,1)) = scenarioRows(k,2);
    r = sdt_solveFrame(model,alpha);
    predictions(:,k) = r.Displacements(candidateDOFs);
end

selected = zeros(0,1);
scores = zeros(0,1);
accumulated = zeros(nScenarios,nScenarios);
for pick = 1:min(numberToAdd,numel(candidateDOFs))
    bestScore = -Inf; bestIndex = NaN; bestAccumulated = [];
    for i = 1:numel(candidateDOFs)
        if ismember(candidateDOFs(i),selected), continue; end
        d = predictions(i,:);
        trial = accumulated;
        for p = 1:nScenarios
            for q = p+1:nScenarios
                trial(p,q) = trial(p,q) + ((d(p)-d(q))/sigma)^2;
            end
        end
        vals = sqrt(trial(triu(true(nScenarios),1)));
        score = min(vals);
        if score > bestScore
            bestScore = score; bestIndex = i; bestAccumulated = trial;
        end
    end
    selected(end+1,1) = candidateDOFs(bestIndex); %#ok<AGROW>
    scores(end+1,1) = bestScore; %#ok<AGROW>
    accumulated = bestAccumulated;
end

recommendation = struct();
recommendation.DOFs = selected;
recommendation.Nodes = ceil(selected/3);
recommendation.Directions = mod(selected-1,3)+1;
recommendation.Scores = scores;
recommendation.ScenarioRows = scenarioRows;
recommendation.Predictions = predictions;
end


% Embedded helper: sensorDOFs.m
function dofs = sdt_sensorDOFs(sensorTable)
%SENSORDOFS Convert [node,direction] rows to global DOF indices.
% Direction: 1=x, 2=y, 3=rotation.
if isempty(sensorTable)
    dofs = zeros(0,1);
    return
end
nodes = sensorTable(:,1);
dirs = sensorTable(:,2);
if any(~isfinite(nodes)) || any(~ismember(dirs,[1 2 3])) || any(nodes < 1) || any(mod(nodes,1)~=0)
    error('shm:InvalidSensor','Sensors require integer node IDs and direction 1, 2, or 3.');
end
dofs = 3*nodes - (3-dirs);
dofs = dofs(:);
end


% Embedded helper: solveFrame.m
function result = sdt_solveFrame(model, alpha)
%SOLVEFRAME Assemble and solve K(alpha)u=F for a 2D frame.

nElements = size(model.Elements,1);
if nargin < 2 || isempty(alpha)
    alpha = zeros(nElements,1);
end
alpha = alpha(:);
if numel(alpha) ~= nElements || any(~isfinite(alpha)) || any(alpha < 0) || any(alpha >= 1)
    error('shm:InvalidDamage','alpha must have one value per element in [0,1).');
end

nDOF = 3*size(model.Nodes,1);
K = sparse(nDOF,nDOF);
elementData = cell(nElements,1);
for e = 1:nElements
    [kg,kl,T,dofs] = sdt_frameElement2D(model,model.Elements(e,:),1-alpha(e));
    K(dofs,dofs) = K(dofs,dofs) + kg;
    elementData{e} = struct('kLocal',kl,'T',T,'dofs',dofs);
end

fixed = unique(model.FixedDOFs(:));
free = setdiff((1:nDOF).',fixed);
Kff = K(free,free);
if isempty(free) || rcond(full(Kff)) < 1e-13
    error('shm:UnstableModel','The free stiffness matrix is singular or nearly singular.');
end

u = zeros(nDOF,1);
u(free) = Kff\model.Loads(free);
reactions = K*u-model.Loads;

localForces = zeros(nElements,6);
for e = 1:nElements
    ed = elementData{e};
    dLocal = ed.T*u(ed.dofs);
    localForces(e,:) = (ed.kLocal*dLocal).';
end

result = struct('Displacements',u,'Reactions',reactions,'K',K, ...
    'FreeDOFs',free,'FixedDOFs',fixed,'LocalForces',localForces, ...
    'ConditionNumber',condest(Kff),'Alpha',alpha);
end

% CE7405 stability helpers. N is positive in tension; compression lowers Kt.
function r = sdt_linearizedStabilitySolve(model, gravityPerFloor, lambda, eiFactor, subdivisions)
validateattributes(gravityPerFloor,{'double'},{'scalar','finite','nonnegative'});
validateattributes(lambda,{'double'},{'scalar','finite','nonnegative'});
validateattributes(eiFactor,{'double'},{'scalar','finite','positive','<=',1});
validateattributes(subdivisions,{'double'},{'scalar','integer','positive'});
nOriginal=size(model.Nodes,1); ne=size(model.Elements,1);
assert(ne*subdivisions<=400,'Reduce subdivisions or frame size to at most 400 stability segments.');
mesh=model; mesh.Nodes=model.Nodes; mesh.Elements=zeros(ne*subdivisions,7);
parents=zeros(ne*subdivisions,1); q=0;
for e=1:ne
    row=model.Elements(e,:); a=row(2); b=row(3);
    ids=[a zeros(1,subdivisions-1) b];
    for j=1:subdivisions-1
        id=size(mesh.Nodes,1)+1;
        xy=model.Nodes(a,2:3)+(j/subdivisions)*(model.Nodes(b,2:3)-model.Nodes(a,2:3));
        mesh.Nodes(end+1,:)=[id xy]; ids(j+1)=id; %#ok<AGROW>
    end
    for j=1:subdivisions
        q=q+1; mesh.Elements(q,:)=[q ids(j) ids(j+1) row(4) row(5) eiFactor*row(6) row(7)];
        parents(q)=e;
    end
end
nd=3*size(mesh.Nodes,1); Fw=zeros(nd,1); Fw(1:numel(model.Loads))=model.Loads;
Fg=zeros(nd,1); nCols=model.nBays+1;
for floor=1:model.nStories
    ids=floor*nCols+(1:nCols); Fg(3*ids-1)=-gravityPerFloor/nCols;
end
K=sparse(nd,nd); meta=cell(size(mesh.Elements,1),1);
for e=1:size(mesh.Elements,1)
    [kg,kl,T,dofs]=sdt_frameElement2D(mesh,mesh.Elements(e,:),1);
    K(dofs,dofs)=K(dofs,dofs)+kg;
    xy=mesh.Nodes(mesh.Elements(e,2:3),2:3); L=norm(xy(2,:)-xy(1,:));
    meta{e}=struct('kl',kl,'T',T,'idx',dofs,'L',L);
end
fixed=unique(mesh.FixedDOFs(:)); free=setdiff((1:nd).',fixed);
Kff=full(K(free,free));
assert(~isempty(free) && rcond(Kff)>1e-13,'Elastic stiffness is singular or poorly conditioned.');
ug=zeros(nd,1); ug(free)=Kff\Fg(free); Kg=sparse(nd,nd);
for e=1:numel(meta)
    p=meta{e}; d=p.T*ug(p.idx); row=mesh.Elements(e,:);
    N=row(4)*row(5)/p.L*(d(4)-d(1)); L=p.L;
    g=zeros(6); idx=[2 3 5 6];
    g(idx,idx)=N/(30*L)*[36 3*L -36 3*L;3*L 4*L^2 -3*L -L^2; ...
        -36 -3*L 36 -3*L;3*L -L^2 -3*L 4*L^2];
    p.kg=g; meta{e}=p; Kg(p.idx,p.idx)=Kg(p.idx,p.idx)+p.T.'*g*p.T;
end
% K phi = -lambda Kg phi, transformed by elastic Cholesky factorization.
Ck=chol(Kff,'lower'); Gff=full(Kg(free,free));
B=-(Ck\Gff)/Ck.'; B=(B+B.')/2;
[V,D]=eig(B); mu=diag(D);
keep=mu>max(1e-14,max(abs(mu))*1e-10);
factors=1./mu(keep); V=V(:,keep); [factors,order]=sort(factors); V=V(:,order);
count=min(5,numel(factors)); factors=factors(1:count);
modes=zeros(nd,count); eigenResiduals=zeros(count,1);
for k=1:count
    v=Ck.'\V(:,k); v=v/max(abs(v)); modes(free,k)=v;
    eigenResiduals(k)=norm((Kff+factors(k)*Gff)*v)/ ...
        max(eps,(norm(Kff,'fro')+factors(k)*norm(Gff,'fro'))*norm(v));
end
if ~isempty(factors) && lambda>=factors(1)*(1-1e-8)
    error('Applied gravity factor %.4g reaches/exceeds the first elastic critical factor %.4g. Reduce lambda or use SHOW 80%% CRITICAL.',lambda,factors(1));
end
Kt=K+lambda*Kg; assert(rcond(full(Kt(free,free)))>1e-12,'Second-order stiffness is near singular.');
F=lambda*Fg+Fw; u1=zeros(nd,1); u2=u1;
u1(free)=Kff\F(free); u2(free)=Kt(free,free)\F(free);
residual=Kt*u2-F;
firstActions=zeros(numel(meta),6); actions=firstActions;
for e=1:numel(meta)
    p=meta{e}; firstActions(e,:)=(p.kl*p.T*u1(p.idx)).';
    actions(e,:)=((p.kl+lambda*p.kg)*p.T*u2(p.idx)).';
end
member=zeros(ne,13);
for e=1:ne
    a=find(parents==e,1,'first'); b=find(parents==e,1,'last');
    member(e,:)=[model.Elements(e,1) firstActions(a,1:3) firstActions(b,4:6) ...
        actions(a,1:3) actions(b,4:6)];
end
roof=model.nStories*nCols+(1:nCols); roofDOFs=3*roof-2;
roof1=mean(u1(roofDOFs)); roof2=mean(u2(roofDOFs));
if abs(roof1)>1e-12, amplification=abs(roof2/roof1); else, amplification=NaN; end
if isempty(factors), limit=max(1,lambda)*2; else, limit=.92*factors(1); end
curveFactors=unique([linspace(0,limit,32) lambda]).';
curve1=zeros(size(curveFactors)); curve2=curve1;
for k=1:numel(curveFactors)
    a=curveFactors(k); ff=a*Fg(free)+Fw(free);
    uu=zeros(nd,1); vv=uu; uu(free)=Kff\ff; vv(free)=(Kff+a*Gff)\ff;
    curve1(k)=mean(uu(roofDOFs)); curve2(k)=mean(vv(roofDOFs));
end
node=zeros(nOriginal,7);
for i=1:nOriginal
    idx=3*i-2:3*i;
    node(i,:)=[model.Nodes(i,1) 1000*u1(idx(1:2)).' u1(idx(3)) ...
        1000*u2(idx(1:2)).' u2(idx(3))];
end
story=zeros(model.nStories,6);
for floor=1:model.nStories
    up=floor*nCols+(1:nCols); lo=(floor-1)*nCols+(1:nCols);
    drift1=u1(3*up-2)-u1(3*lo-2); drift2=u2(3*up-2)-u2(3*lo-2);
    story(floor,:)=[floor 1000*mean(drift1) 1000*mean(drift2) ...
        100*mean(drift1)/model.StoryHeight 100*mean(drift2)/model.StoryHeight ...
        100*max(abs(drift2))/model.StoryHeight];
end
supported=unique(ceil(fixed/3)); reactions=zeros(numel(supported),4);
for i=1:numel(supported)
    n=supported(i); reactions(i,:)=[n residual(3*n-2:3*n).'];
end
r=struct('Mesh',mesh,'K',K,'Kg',Kg,'Free',free,'Fixed',fixed,'First',u1,'Second',u2, ...
    'Factors',factors,'Modes',modes,'EigenResiduals',eigenResiduals,'LoadFactor',lambda, ...
    'GravityPerFloor',gravityPerFloor,'EIFactor',eiFactor,'Subdivisions',subdivisions, ...
    'ReferenceGravity',Fg,'LateralLoad',Fw,'RoofFirst',roof1,'RoofSecond',roof2, ...
    'Amplification',amplification,'EquilibriumResidual',norm(residual(free))/max(1,norm(F(free))), ...
    'SegmentActions',actions,'CurveFactors',curveFactors,'CurveFirst',curve1,'CurveSecond',curve2);
r.StoryDrift=array2table(story,'VariableNames',{'Story','MeanDrift1_mm','MeanDrift2_mm', ...
    'MeanDrift1_percent','MeanDrift2_percent','MaxColumnDrift2_percent'});
r.NodeResponse=array2table(node,'VariableNames',{'Node','Ux1_mm','Uy1_mm','Rz1_rad','Ux2_mm','Uy2_mm','Rz2_rad'});
r.MemberActions=array2table(member,'VariableNames',{'Element','Ni1_N','Vi1_N','Mi1_Nm','Nj1_N','Vj1_N','Mj1_Nm', ...
    'Ni2_N','Vi2_N','Mi2_Nm','Nj2_N','Vj2_N','Mj2_Nm'});
r.Reactions=array2table(reactions,'VariableNames',{'Node','Rx_N','Ry_N','Mz_Nm'});
end

function h = sdt_drawFrame(ax,mesh,u,scale,color,lineStyle)
% Front-view Hermite curve in the displayed current chord frame.
% Removes rigid chord rotation before interpolating relative end rotations.
h=[];
for e=1:size(mesh.Elements,1)
    row=mesh.Elements(e,:); xy=mesh.Nodes(row(2:3),2:3);
    ix=[3*row(2)-2:3*row(2) 3*row(3)-2:3*row(3)]; d=scale*u(ix);
    a=xy(1,:)+d(1:2).'; b=xy(2,:)+d(4:5).'; chord=b-a; ell=norm(chord);
    if ell<=eps, continue; end
    initial=xy(2,:)-xy(1,:);
    beta=atan2(initial(1)*chord(2)-initial(2)*chord(1),dot(initial,chord));
    thetaI=d(3)-beta; thetaJ=d(6)-beta;
    t=linspace(0,1,16).';
    transverse=ell*((t-2*t.^2+t.^3)*thetaI+(-t.^2+t.^3)*thetaJ);
    normal=[-chord(2) chord(1)]/ell;
    points=a+t*chord+transverse*normal;
    line=plot(ax,points(:,1),points(:,2),'Color',color,'LineStyle',lineStyle,'LineWidth',1.9);
    if isempty(h), h=line; else, line.HandleVisibility='off'; end
end
end

function report = sdt_verifyStability()
E=200e9; I=1e-4; A=.02; L=4;
m=struct('Nodes',[1 0 0;2 0 L],'Elements',[1 1 2 E A I 1], ...
    'FixedDOFs',[1 2 4],'Loads',zeros(6,1),'nStories',1,'nBays',0,'StoryHeight',L);
r=sdt_stabilitySolve(m,1000,1,1,8);
pinned=pi^2*E*I/(L^2*1000); errPinned=abs(r.Factors(1)/pinned-1);
assert(errPinned<.001,'Pinned-column Euler benchmark failed.');
m.FixedDOFs=[1 2 3]; Pcr=pi^2*E*I/(4*L^2); P=.3*Pcr; H=1000;
m.Loads(4)=H; r=sdt_stabilitySolve(m,P,1,1,16);
errFixedFree=abs(r.Factors(1)/(Pcr/P)-1);
k=sqrt(P/(E*I)); expected=H/P*(tan(k*L)/k-L);
errDeflection=abs(r.Second(4)/expected-1);
assert(errFixedFree<.001,'Cantilever Euler benchmark failed.');
assert(errDeflection<.005,'Corotational beam-column lateral deflection benchmark failed.');
assert(r.EquilibriumResidual<1e-8 && max(r.EigenResiduals)<1e-8,'Matrix residual verification failed.');
report=sprintf(['PASSED: pinned-column Euler error %.5g%%\n' ...
    'PASSED: cantilever Euler error %.5g%%\n' ...
    'PASSED: beam-column tip-deflection error %.5g%%\n' ...
    'Free-DOF residual %.3e\nReference eigenbuckling and corotational Newton equilibrium checked; elastic material only.'], ...
    100*errPinned,100*errFixedFree,100*errDeflection,r.EquilibriumResidual);
end


function s = sdt_environmentSection(Ec,Ag,Ig,o,ecLoss,steelLoss)
% Uncracked symmetric transformed RC section. SI units throughout.
% Corroded steel is removed as void; concrete Ac,Ic are kept fixed.
% Uniform loss on every bar; no bond/slip, spalling, cracking or strength model.
validateattributes(Ec,{'double'},{'scalar','real','finite','positive'});
validateattributes(Ag,{'double'},{'scalar','real','finite','positive'});
validateattributes(Ig,{'double'},{'scalar','real','finite','positive'});
validateattributes(o.Es,{'double'},{'scalar','real','finite','positive'});
validateattributes(o.nFace,{'double'},{'scalar','real','finite','integer','positive'});
validateattributes(o.d,{'double'},{'scalar','real','finite','positive'});
validateattributes(o.a,{'double'},{'scalar','real','finite','positive'});
validateattributes(ecLoss,{'double'},{'scalar','real','finite','>=',0,'<=',.8});
validateattributes(steelLoss,{'double'},{'scalar','real','finite','>=',0,'<=',.8});
h=sqrt(12*Ig/Ag); b=Ag/h;
assert(o.a>=o.d/2 && h-2*o.a>=o.d,'Bars do not fit vertically. Check bar-centre cover, diameter, A and I.');
if b<2*o.a
    error('sdt:EnvironmentSectionWidth', ...
        ['A=%.8g m^2 and I=%.8g m^4 imply rectangular b=%.3f mm, h=%.3f mm. ' ...
         'Bar-centre cover %.3f mm requires b >= %.3f mm. ' ...
         'For the RC demo, click RC 400 mm PRESET in the Exposure tab BEFORE assigning losses; ' ...
         'or supply consistent RC E/A/I and reinforcement dimensions.'], ...
         Ag,Ig,1000*b,1000*h,1000*o.a,2000*o.a);
end
if o.nFace>1
    assert((b-2*o.a)/(o.nFace-1)>=o.d,'Bars overlap within a face. Reduce count/diameter or increase width.');
end
Ab=pi*o.d^2/4; Ib=pi*o.d^4/64; z=h/2-o.a; n=2*o.nFace;
As0=n*Ab; Is0=n*(Ib+Ab*z^2);
Ac=Ag-As0; Ic=Ig-Is0;
assert(Ac>0 && Ic>0,'Invalid reinforcement: remaining concrete area/inertia must be positive.');
ret=1-steelLoss;
As=ret*As0;
Is=n*(ret^2*Ib+ret*Ab*z^2);
EcD=(1-ecLoss)*Ec;
EA0=Ec*Ac+o.Es*As0; EI0=Ec*Ic+o.Es*Is0;
EA=EcD*Ac+o.Es*As; EI=EcD*Ic+o.Es*Is;
s=struct('Ec0_GPa',Ec/1e9,'Es_GPa',o.Es/1e9,'b_mm',1000*b,'h_mm',1000*h, ...
    'BarDia_mm',1000*o.d,'BarsPerFace',o.nFace,'BarCentroidCover_mm',1000*o.a, ...
    'Ec_loss_percent',100*ecLoss,'As_loss_percent',100*steelLoss, ...
    'Ag_m2',Ag,'Ig_m4',Ig,'As0_m2',As0,'As_m2',As, ...
    'Is0_m4',Is0,'Is_m4',Is,'Ac_m2',Ac,'Ic_m4',Ic, ...
    'EA0_N',EA0,'EA_N',EA,'EI0_Nm2',EI0,'EI_Nm2',EI, ...
    'EA_ratio',EA/EA0,'EI_ratio',EI/EI0, ...
    'Esolver_Pa',Ec,'Aeq_m2',EA/Ec,'Ieq_m4',EI/Ec);
end

function [healthy,scenario,report] = sdt_environmentBuild(base,data,o)
assert(~isempty(base),'Generate a frame first.');
ne=size(base.Elements,1); assert(isequal(size(data),[ne 6]),'Invalid exposure table size.');
healthy=base; scenario=base; records=cell(ne,1);
labels={'Reference','Sulfate','Chloride / seawater','Combined'};
for e=1:ne
    row=base.Elements(e,:);
    assert(isequal(data{e,1},row(1)),'Exposure row IDs do not match frame members.');
    exposure=char(data{e,4}); assert(any(strcmp(exposure,labels)),'Unknown exposure label.');
    ce=data{e,5}; ar=data{e,6};
    validateattributes(ce,{'double'},{'scalar','real','finite','>=',0,'<=',80});
    validateattributes(ar,{'double'},{'scalar','real','finite','>=',0,'<=',80});
    if strcmp(exposure,'Reference'), assert(ce==0 && ar==0,'Reference members must have zero losses.'); end
    if strcmp(exposure,'Sulfate'), assert(ar==0,'Use Combined to include steel loss with sulfate E loss.'); end
    if strcmp(exposure,'Chloride / seawater'), assert(ce==0,'Use Combined to include concrete E loss with chloride steel loss.'); end
    s=sdt_environmentSection(row(4),row(5),row(6),o,ce/100,ar/100);
    healthy.Elements(e,5:6)=[s.EA0_N/row(4) s.EI0_Nm2/row(4)];
    scenario.Elements(e,5:6)=[s.Aeq_m2 s.Ieq_m4];
    r=struct('Element',row(1),'Type',string(data{e,2}),'Story',data{e,3},'Exposure',string(exposure));
    names=fieldnames(s);
    for k=1:numel(names), r.(names{k})=s.(names{k}); end
    records{e}=r;
end
report=struct2table(vertcat(records{:}));
healthy.Environment='Healthy uncracked RC composite';
scenario.Environment='Exposure-related assumed deterioration';
end

function value = sdt_environmentPercent(reference,scenario)
if ~isfinite(reference) || ~isfinite(scenario) || abs(reference)<1e-12, value=NaN;
else, value=100*(scenario-reference)/abs(reference); end
end

function comparison = sdt_environmentCompare(model,h,s)
nCols=model.nBays+1; roof=model.nStories*nCols+(1:nCols);
roofH=mean(h.Displacements(3*roof-2)); roofS=mean(s.Displacements(3*roof-2));
dh=0; ds=0;
for k=1:model.nStories
    up=k*nCols+(1:nCols); lo=(k-1)*nCols+(1:nCols);
    dh=max(dh,max(abs(h.Displacements(3*up-2)-h.Displacements(3*lo-2)))/model.StoryHeight);
    ds=max(ds,max(abs(s.Displacements(3*up-2)-s.Displacements(3*lo-2)))/model.StoryHeight);
end
names={'Roof mean horizontal displacement';'Maximum column drift ratio';'Load compliance F^T u'};
units={'mm';'%';'N m'};
x=[1000*roofH;100*dh;model.Loads.'*h.Displacements];
y=[1000*roofS;100*ds;model.Loads.'*s.Displacements];
pct=arrayfun(@sdt_environmentPercent,x,y);
comparison=table(names,units,x,y,pct,'VariableNames', ...
    {'Metric','Unit','HealthyRC','Scenario','Change_percent'});
end

function lines = sdt_environmentAudit(healthy,scenario,report,e)
validateattributes(e,{'double'},{'scalar','integer','>=',1,'<=',height(report)});
s=report(e,:); er=healthy.Elements(e,:); xy=healthy.Nodes(er(2:3),2:3); L=norm(xy(2,:)-xy(1,:));
[kg0,k0,T,dofs]=sdt_frameElement2D(healthy,er,1);
[kg,k]=sdt_frameElement2D(scenario,scenario.Elements(e,:),1);
lines={sprintf('ELEMENT %d | %s | story %d | %s | L = %.12g m',s.Element,char(s.Type),s.Story,char(s.Exposure)); ...
    'MODEL: uncracked rectangular RC, symmetric top/bottom bars, perfect bond.'; ...
    'Area loss uniformly shrinks each bar at fixed bar-centre location. Removed steel remains void.'; ...
    sprintf('Ag = %.12g m^2; Ig = %.12g m^4; b = %.12g m; h = %.12g m',s.Ag_m2,s.Ig_m4,s.b_mm/1000,s.h_mm/1000); ...
    sprintf('Ec0 = %.12g Pa; Es = %.12g Pa; bars = %d per face; d = %.12g m; a = %.12g m',s.Ec0_GPa*1e9,s.Es_GPa*1e9,s.BarsPerFace,s.BarDia_mm/1000,s.BarCentroidCover_mm/1000); ...
    sprintf('ASSUMED concrete E loss = %.8g%%; steel area loss = %.8g%%',s.Ec_loss_percent,s.As_loss_percent); ...
    'Ab=pi*d^2/4; Ib=pi*d^4/64; z=h/2-a; n=2*barsPerFace'; ...
    'As0=n*Ab; Is0=n*(Ib+Ab*z^2); Ac=Ag-As0; Ic=Ig-Is0'; ...
    sprintf('As0 = %.12g; Is0 = %.12g; Ac = %.12g; Ic = %.12g (SI)',s.As0_m2,s.Is0_m4,s.Ac_m2,s.Ic_m4); ...
    'r=AsLoss/100; As=(1-r)*As0; Is=n*((1-r)^2*Ib+(1-r)*Ab*z^2)'; ...
    sprintf('As = %.12g m^2; Is = %.12g m^4',s.As_m2,s.Is_m4); ...
    'EcD=(1-EcLoss/100)*Ec0; EA=EcD*Ac+Es*As; EI=EcD*Ic+Es*Is'; ...
    sprintf('EA0 = %.12g N; EA = %.12g N; ratio = %.12g',s.EA0_N,s.EA_N,s.EA_ratio); ...
    sprintf('EI0 = %.12g N*m^2; EI = %.12g N*m^2; ratio = %.12g',s.EI0_Nm2,s.EI_Nm2,s.EI_ratio); ...
    'Solver storage: E=Ec0 (constant reference), Aeq=EA/Ec0, Ieq=EI/Ec0. Equivalent properties are not remaining geometric concrete A,I.'; ...
    sprintf('Esolver = %.12g Pa; Aeq = %.12g m^2; Ieq = %.12g m^4',s.Esolver_Pa,s.Aeq_m2,s.Ieq_m4); ...
    sprintf('DOF map = %s; local order [axial_i transverse_i rotation_i axial_j transverse_j rotation_j]',mat2str(dofs)); ...
    'kLocal terms: EA/L; 12EI/L^3; 6EI/L^2; 4EI/L; 2EI/L.'; ...
    'Matrices are listed below, one row per line.'; ...
    'Assemble K(dofs,dofs)+=kGlobal; solve Kff*u_f=F_f, fixed displacements=0.'; ...
    'Stability tab additionally multiplies Ieq by its user-selected effective-EI factor.'};
matrices={k0,k,T,kg0,kg};
names={'Healthy kLocal','Scenario kLocal','Transformation T','Healthy kGlobal','Scenario kGlobal'};
for j=1:numel(matrices)
    lines{end+1,1}=names{j}; %#ok<AGROW>
    for row=1:size(matrices{j},1)
        lines{end+1,1}=sprintf('% .8e  ',matrices{j}(row,:)); %#ok<AGROW>
    end
end
end

function report = sdt_verifyEnvironment()
% Pure numerical checks. No GUI or App Designer session is required.
o=struct('Es',200e9,'nFace',2,'d',.016,'a',.05);
Ec=30e9; A=.16; I=.4^4/12;
s0=sdt_environmentSection(Ec,A,I,o,0,0);
s=sdt_environmentSection(Ec,A,I,o,.1,.2);
assert(abs(s0.EA_ratio-1)<1e-14 && abs(s0.EI_ratio-1)<1e-14,'Zero-loss identity failed.');
% Independent direct rectangle-with-four-circular-bars calculation.
As=4*pi*.016^2/4; Is=4*(pi*.016^4/64+pi*.016^2/4*.15^2);
EA0=30e9*(.16-As)+200e9*As;
EI0=30e9*(.4^4/12-Is)+200e9*Is;
EAd=.9*30e9*(.16-As)+200e9*.8*As;
Id=4*(pi*(.016*sqrt(.8))^4/64+pi*(.016*sqrt(.8))^2/4*.15^2);
EId=.9*30e9*(.4^4/12-Is)+200e9*Id;
assert(abs(s.EA_N-EAd)/EAd<1e-12 && abs(s.EI_Nm2-EId)/EId<1e-12,'Independent section check failed.');
assert(abs(s0.EA0_N-EA0)/EA0<1e-12 && abs(s0.EI0_Nm2-EI0)/EI0<1e-12);
assert(s.EA_ratio<1 && s.EI_ratio<1 && abs(s.EA_ratio-s.EI_ratio)>1e-4,'Separate EA/EI update failed.');
% Cantilever closed-form axial + transverse tip displacement.
L=3; m=struct('Nodes',[1 0 0;2 L 0], ...
    'Elements',[1 1 2 Ec s.Aeq_m2 s.Ieq_m4 2], ...
    'FixedDOFs',[1 2 3],'Loads',[0;0;0;1000;2000;0]);
r=sdt_solveFrame(m,[]);
assert(abs(r.Displacements(4)-1000*L/EAd)<1e-12);
assert(abs(r.Displacements(5)-2000*L^3/(3*EId))<1e-12);
assert(norm(r.Reactions(r.FreeDOFs))<1e-7);
% Unselected members, repeatability, model recovery, and stiffness assembly.
b=sdt_generatePortalFrame(2,1,6,3.5,Ec,A,I); ne=size(b.Elements,1);
data=cell(ne,6);
for e=1:ne
    row=b.Elements(e,:); yy=b.Nodes(row(2:3),3);
    if row(7)==1, kind='Column'; else, kind='Beam'; end
    data(e,:)={e,kind,round(max(yy)/3.5),'Reference',0,0};
end
data(1,4:6)={'Combined',10,20};
[h,d,t]=sdt_environmentBuild(b,data,o);
[h2,d2,t2]=sdt_environmentBuild(b,data,o);
assert(isequal(d.Elements,d2.Elements) && isequal(t,t2),'Reapplying losses compounded them.');
assert(isequal(h.Elements(2:end,:),d.Elements(2:end,:)),'Unselected member changed.');
assert(abs(d.Elements(1,4)*d.Elements(1,5)-EAd)/EAd<1e-12);
assert(abs(d.Elements(1,4)*d.Elements(1,6)-EId)/EId<1e-12);
rh=sdt_solveFrame(h,[]); rd=sdt_solveFrame(d,[]);
assert(norm(rd.K-rd.K.','fro')/norm(rd.K,'fro')<1e-13,'K symmetry failed.');
assert(b.Loads.'*rd.Displacements>=b.Loads.'*rh.Displacements-1e-10,'Compliance unexpectedly decreased.');
sh=sdt_stabilitySolve(h,500e3,.5,1,4); sd=sdt_stabilitySolve(d,500e3,.5,1,4);
assert(sd.Factors(1)<sh.Factors(1),'Deteriorated critical factor did not decrease in this test.');
assert(sd.EquilibriumResidual<1e-9);
% Reject invalid inputs rather than silently using meaningless damage.
bad=o; bad.d=.5; rejected=false;
try, sdt_environmentSection(Ec,A,I,bad,0,0); catch, rejected=true; end
assert(rejected,'Invalid reinforcement geometry was not rejected.');
% Regression: the original app's steel-like default section is too narrow
% for 50 mm bar-centre cover and must yield an actionable error.
rejected=false;
try
    sdt_environmentSection(200e9,.02,8e-5,o,0,0);
catch ME
    rejected=strcmp(ME.identifier,'sdt:EnvironmentSectionWidth') && contains(ME.message,'RC 400 mm PRESET');
end
assert(rejected,'Original default-section guard/preset guidance failed.');
rejected=false; data(1,4:6)={'Reference',10,20};
try, sdt_environmentBuild(b,data,o); catch, rejected=true; end
assert(rejected,'Nonzero losses under Reference were not rejected.');
report=table({'Healthy EA';'Scenario EA';'Healthy EI';'Scenario EI'; ...
    'Cantilever axial displacement';'Cantilever transverse displacement'; ...
    'Healthy critical gravity factor';'Scenario critical gravity factor'}, ...
    [EA0;EAd;EI0;EId;r.Displacements(4);r.Displacements(5);sh.Factors(1);sd.Factors(1)], ...
    'VariableNames',{'Quantity','Value_SI'});
disp(report);
fprintf('PASS: section identity, hand formulas, separate EA/EI, cantilever, repeatability, unaffected members, equilibrium, stability and invalid inputs.\n');
end


function [f,k,actions] = sdt_corotElement(xy,row,u)
% Small-strain corotational Euler-Bernoulli beam with analytic consistent tangent.
% Basic deformations: chord extension and two rotations relative to the chord.
% Engineering axial strain; elastic EI. External loads are dead global nodal loads.
D=[-1 0 0 1 0 0;0 -1 0 0 1 0];
v0=(xy(2,:)-xy(1,:)).'; L=norm(v0); delta=D*u; v=v0+delta; ell=norm(v);
if ~(L>0 && ell>1e-8*L)
    error('sdt:CorotCollapsed','A corotational element has collapsed to zero length.');
end
n=v/ell; c=n(1); s=n(2);
beta=atan2(v0(1)*v(2)-v0(2)*v(1),v0.'*v);
% Stable alternative to ell-L, avoiding subtractive cancellation at small load.
extension=(2*(v0.'*delta)+delta.'*delta)/(ell+L);
qv=[extension;u(3)-beta;u(6)-beta];
E=row(4); A=row(5); I=row(6);
kb=[E*A/L 0 0;0 4*E*I/L 2*E*I/L;0 2*E*I/L 4*E*I/L];
B=zeros(3,6); B(1,:)=n.'*D; g=[-s c]*D/ell;
B(2,:)=-g; B(2,3)=1; B(3,:)=-g; B(3,6)=1;
q=kb*qv; f=B.'*q;
if nargout>=2
    Hl=D.'*(eye(2)-n*n.')*D/ell;
    Hb=D.'*[2*c*s s*s-c*c;s*s-c*c -2*c*s]*D/ell^2;
    k=B.'*kb*B+q(1)*Hl-(q(2)+q(3))*Hb;
    k=(k+k.')/2;
end
if nargout>=3
    actions=[-q(1);(q(2)+q(3))/ell;q(2);q(1);-(q(2)+q(3))/ell;q(3)];
end
end

function [f,K,actions] = sdt_corotAssemble(model,u)
nd=numel(u); ne=size(model.Elements,1); f=zeros(nd,1);
needTangent=nargout>=2;
if needTangent, K=sparse(nd,nd); end
if nargout>=3, actions=zeros(ne,6); end
for e=1:ne
    row=model.Elements(e,:); xy=model.Nodes(row(2:3),2:3);
    ix=[3*row(2)-2:3*row(2) 3*row(3)-2:3*row(3)];
    if nargout>=3
        [fe,ke,ae]=sdt_corotElement(xy,row,u(ix)); actions(e,:)=ae.';
    elseif needTangent
        [fe,ke]=sdt_corotElement(xy,row,u(ix));
    else
        fe=sdt_corotElement(xy,row,u(ix));
    end
    f(ix)=f(ix)+fe;
    if needTangent, K(ix,ix)=K(ix,ix)+ke; end
end
end

function nr = sdt_newtonSolve(model,F,fixed,uStart)
% Full Newton-Raphson: current residual and tangent are rebuilt every iteration.
% Adaptive load increments and backtracking. Only converged states are returned.
nd=numel(F); free=setdiff((1:nd).',fixed(:));
if nargin<4 || isempty(uStart), u=zeros(nd,1); else, u=uStart; end
assert(numel(u)==nd && ~isempty(free),'Invalid nonlinear model or start vector.');
u(fixed)=0; Fstart=sdt_corotAssemble(model,u); Fstart(fixed)=0;
lengthScale=max([1,max(model.Nodes(:,2))-min(model.Nodes(:,2)),max(model.Nodes(:,3))-min(model.Nodes(:,3))]);
weights=ones(nd,1); weights(3:3:end)=1/lengthScale;
tol=1e-8; maxIterations=35; fraction=0; step=.125; attempt=0;
log=zeros(0,5);
while fraction<1-1e-12
    attempt=attempt+1; target=min(1,fraction+step);
    Ft=(1-target)*Fstart+target*F;
    denominator=max(1,norm(weights(free).*Ft(free)));
    trial=u; converged=false;
    for iteration=0:maxIterations
        [f,K]=sdt_corotAssemble(model,trial);
        residual=Ft-f; rel=norm(weights(free).*residual(free))/denominator;
        if rel<=tol
            log(end+1,:)=[attempt target iteration rel 0]; %#ok<AGROW>
            converged=true; break;
        end
        if iteration==maxIterations, break; end
        Kff=K(free,free);
        [~,flag]=chol(Kff);
        if flag~=0, break; end % restrict load control to stable tangent states
        du=zeros(nd,1); du(free)=Kff\residual(free);
        if any(~isfinite(du)), break; end
        accepted=false; a=1;
        for search=1:14
            candidate=trial+a*du;
            try
                fc=sdt_corotAssemble(model,candidate);
                rc=norm(weights(free).*(Ft(free)-fc(free)))/denominator;
            catch ME
                if strcmp(ME.identifier,'sdt:CorotCollapsed'), rc=Inf; else, rethrow(ME); end
            end
            if rc<(1-1e-4*a)*rel || rc<=tol
                accepted=true; break;
            end
            a=a/2;
        end
        log(end+1,:)=[attempt target iteration rel a]; %#ok<AGROW>
        if ~accepted, break; end
        trial=candidate;
    end
    if converged
        u=trial; fraction=target;
        if iteration<=6, step=min(.5,1.5*step); end
    else
        step=step/2;
        if step<1/2048
            error('sdt:NewtonNoConvergence', ...
                ['Newton-Raphson did not converge near load fraction %.5g (residual %.3e). ' ...
                 'Reduce gravity/lateral load or refine the mesh. No unconverged response is reported.'],target,rel);
        end
    end
end
[f,K,actions]=sdt_corotAssemble(model,u);
relative=norm(weights(free).*(F(free)-f(free)))/max(1,norm(weights(free).*F(free)));
assert(relative<=tol,'Final Newton equilibrium residual exceeds tolerance.');
nr=struct('U',u,'Internal',f,'Tangent',K,'Actions',actions,'Residual',relative, ...
    'Tolerance',tol,'History',array2table(log,'VariableNames', ...
    {'Attempt','LoadFraction','Iteration','RelativeResidual','LineSearchFactor'}));
end

function r = sdt_stabilitySolve(model,gravityPerFloor,lambda,eiFactor,subdivisions)
% First order and reference elastic eigenbuckling are retained independently.
% Nonlinear equilibrium uses current-geometry corotational Newton-Raphson.
r=sdt_linearizedStabilitySolve(model,gravityPerFloor,lambda,eiFactor,subdivisions);
F=r.LateralLoad+lambda*r.ReferenceGravity;
nr=sdt_newtonSolve(r.Mesh,F,r.Fixed,[]);
r.Second=nr.U; r.SegmentActions=nr.Actions; r.TangentK=nr.Tangent;
r.NewtonHistory=nr.History; r.NewtonTolerance=nr.Tolerance;
r.EquilibriumResidual=nr.Residual; r.Solver='Corotational full Newton-Raphson with load stepping and line search';
nCols=model.nBays+1; roof=model.nStories*nCols+(1:nCols);
r.RoofSecond=mean(nr.U(3*roof-2));
if abs(r.RoofFirst)>1e-12, r.Amplification=abs(r.RoofSecond/r.RoofFirst); else, r.Amplification=NaN; end
n=size(model.Nodes,1);
r.NodeResponse.Ux2_mm=1000*nr.U(1:3:3*n);
r.NodeResponse.Uy2_mm=1000*nr.U(2:3:3*n);
r.NodeResponse.Rz2_rad=nr.U(3:3:3*n);
for floor=1:model.nStories
    up=floor*nCols+(1:nCols); lo=(floor-1)*nCols+(1:nCols);
    drift=nr.U(3*up-2)-nr.U(3*lo-2);
    r.StoryDrift{floor,[3 5 6]}=[1000*mean(drift),100*mean(drift)/model.StoryHeight,100*max(abs(drift))/model.StoryHeight];
end
for e=1:size(model.Elements,1)
    a=(e-1)*subdivisions+1; b=e*subdivisions;
    r.MemberActions{e,8:13}=[nr.Actions(a,1:3) nr.Actions(b,4:6)];
end
residual=nr.Internal-F; supported=unique(ceil(r.Fixed/3));
for j=1:numel(supported)
    node=supported(j); r.Reactions{j,2:4}=residual(3*node-2:3*node).';
end
if isempty(r.Factors), upper=max(1,lambda)*2; else, upper=.85*r.Factors(1); end
r.CurveFactors=unique([linspace(0,upper,12) lambda]).';
r.CurveFirst=zeros(size(r.CurveFactors)); r.CurveSecond=nan(size(r.CurveFactors));
previous=[]; r.CurveConverged=false(size(r.CurveFactors));
for j=1:numel(r.CurveFactors)
    q=r.CurveFactors(j); ff=r.LateralLoad+q*r.ReferenceGravity;
    first=zeros(size(ff)); first(r.Free)=r.K(r.Free,r.Free)\ff(r.Free);
    r.CurveFirst(j)=mean(first(3*roof-2));
    if abs(q-lambda)<1e-12
        r.CurveSecond(j)=r.RoofSecond; r.CurveConverged(j)=true; previous=nr.U;
    else
        try
            response=sdt_newtonSolve(r.Mesh,ff,r.Fixed,previous);
            r.CurveSecond(j)=mean(response.U(3*roof-2));
            r.CurveConverged(j)=true; previous=response.U;
        catch ME
            if strcmp(ME.identifier,'sdt:NewtonNoConvergence')
                % Unconverged path points remain NaN; no linear substitute.
                previous=[];
            else, rethrow(ME); end
        end
    end
end
end

function report=sdt_verifyNewtonRaphson()
% No GUI required: consistent tangent, rigid motion, linear limit and P-Delta.
E=200e9; A=.02; I=1e-4;
xy=[.2 -.3;3.5 1.7]; row=[1 1 2 E A I 1];
u=[.01;-.02;.015;.03;.005;-.01]; [~,K]=sdt_corotElement(xy,row,u);
fd=zeros(6); step=1e-6;
for j=1:6
    d=zeros(6,1); d(j)=step;
    fd(:,j)=(sdt_corotElement(xy,row,u+d)-sdt_corotElement(xy,row,u-d))/(2*step);
end
tangentError=norm(K-fd,'fro')/norm(K,'fro'); assert(tangentError<1e-7);
angle=.37; R=[cos(angle) -sin(angle);sin(angle) cos(angle)];
rotated=xy*R.'+[1 -2]; rigid=[(rotated(1,:)-xy(1,:)).';angle;(rotated(2,:)-xy(2,:)).';angle];
f=sdt_corotElement(xy,row,rigid); assert(norm(f)<1e-4,'Rigid motion generated internal force.');
m=sdt_generatePortalFrame(1,1,6,3.5,30e9,.16,.4^4/12);
m.Loads=m.Loads*.0005; lin=sdt_solveFrame(m,[]);
nr=sdt_newtonSolve(m,m.Loads,m.FixedDOFs,[]);
linearError=norm(nr.U-lin.Displacements)/norm(lin.Displacements); assert(linearError<1e-5);
L=4; Pcr=pi^2*E*I/(4*L^2); P=.3*Pcr; H=1000; ns=16;
nodes=[(1:ns+1).' zeros(ns+1,1) linspace(0,L,ns+1).'];
elements=[(1:ns).' (1:ns).' (2:ns+1).' repmat([E A I 1],ns,1)];
F=zeros(3*(ns+1),1); F(end-2)=H; F(end-1)=-P;
m=struct('Nodes',nodes,'Elements',elements,'FixedDOFs',[1 2 3],'Loads',F);
nr=sdt_newtonSolve(m,F,m.FixedDOFs,[]);
k=sqrt(P/(E*I)); expected=H/P*(tan(k*L)/k-L);
beamError=abs(nr.U(end-2)/expected-1); assert(beamError<.005);
report=table({'Consistent tangent FD relative error';'Rigid-motion force norm (N, Nm)'; ...
    'Small-load linear limit relative error';'Beam-column reference relative error';'Final NR relative residual'}, ...
    [tangentError;norm(f);linearError;beamError;nr.Residual],'VariableNames',{'Check','Value'});
disp(report); fprintf('PASS: corotational Newton-Raphson numerical checks.\n');
end


function sdt_drawSupports(ax,model,span)
% Fixed: horizontal base plate and hatch, not a pin-support triangle.
supported=unique(ceil(model.FixedDOFs(:)/3));
for n=supported.'
    x=model.Nodes(n,2); y=model.Nodes(n,3);
    fixed=ismember(3*n-2:3*n,model.FixedDOFs);
    halfWidth=.026*span; drop=.018*span;
    if all(fixed)
        plot(ax,x+[-halfWidth halfWidth],[y y],'-','Color',[.05 .05 .05], ...
            'LineWidth',3.5,'HandleVisibility','off','Tag','FixedSupportPlate');
        for q=linspace(-halfWidth,halfWidth,6)
            plot(ax,[x+q x+q-.010*span],[y y-drop],'-','Color',[.15 .15 .15], ...
                'LineWidth',1.2,'HandleVisibility','off','Tag','FixedSupportHatch');
        end
        text(ax,x,y-.029*span,'FIXED','HorizontalAlignment','center', ...
            'VerticalAlignment','top','FontSize',8,'FontWeight','bold','Clipping','on','Tag','FixedSupportLabel');
    else
        labels={'Ux','Uy','Rz'};
        text(ax,x,y-.025*span,['Restrained: ' strjoin(labels(fixed),',')], ...
            'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Clipping','on');
    end
end
end
