function wahi_gui()
% WAH-i GUI — Microphone Array Design & Optimisation (iWAH)
%
% Purpose
% -------
% This application provides an interactive interface to design and optimise 3D
% microphone-array geometries for time-difference-of-arrival (TDOA) localisation.
% It evaluates localisation performance over a user-defined 3D “Field of Accuracy”
% (FoA) and applies an iterative, heuristic optimiser (WAH-i) to improve pass rate
% and error statistics under geometric constraints.
%
% What this GUI does
% ------------------
% 1) Geometry creation and management
%    - Generate random non-coplanar arrays under radius and spacing constraints.
%    - Load/save array geometry from/to CSV (Nx3: [x y z] in metres).
%    - Manually edit individual microphone positions (sliders + edit boxes).
%    - Re-centre arrays to keep the centroid at the origin.
%
% 2) Field-of-Accuracy (FoA) evaluation
%    - Builds a spherical-shell grid (Rin → R) with configurable spacing.
%    - Allows selecting which octants (“Optimisation Octants”) contribute to the grid.
%    - Runs the localisation pipeline at each grid point:
%         simulate -> TDOA estimation (xcorr) -> global solve (lsqnonlin)
%    - Produces performance statistics:
%         pass rate, mean error, 95th percentile, max error, failure counts.
%    - Visualises the FoA as a 3D colour field (position error or angular error).
%
% 3) Optimisation (WAH-i)
%    - Iteratively perturbs microphone positions using a probe-direction rule set.
%    - Accepts updates that improve a scalar score:
%         J = passRate − λ * (meanErr / failPenalty)
%    - Includes stall detection and “shake/reset” steps to escape local optima.
%    - Logs a complete run history (H) for reproducibility and post-analysis.
%
% 4) Reporting and export
%    - Export Run: saves run history (H) and runInfo to a MAT-file (-v7.3).
%    - Report: generates a single-run report figure (before vs after, drift, etc.)
%      with optional PNG export.
%    - Export Workspace: assigns key variables into the base workspace.
%
% Key concepts and definitions
% ----------------------------
% - FoA shell:
%     A hollow spherical region (Rin <= r <= R) sampled on a 3D Cartesian grid.
% - Threshold:
%     threshold_cm = best_cm + tol_cm. “Pass” points satisfy error <= threshold_cm.
% - Failure handling:
%     Points that fail simulation/TDOA/solve are penalised (failPenalty_cm).
%     Grid points too close to microphones can be excluded (minDistToAnyMic) or
%     forced to penalty when within epsDist.
% - Metrics:
%     * Position error (cm): Euclidean distance between estimated and true source.
%     * Angular error (deg): crude azimuth/elevation deviation relative to mic 1.
%
% Requirements
% ------------
% - MATLAB (App Designer UI components: uifigure/uigridlayout/uiaxes/etc.)
% - Optimization Toolbox (lsqnonlin)
% - BatCallLocaliser.m on the MATLAB path (incldued in src/)
%
% File I/O conventions
% --------------------
% - Geometry CSV:
%     Nx3 numeric matrix with columns [x y z] in metres. Values must be finite.
%     Loaded geometries are re-centred to zero mean.
% - Run export MAT:
%     Saves variables:
%       H       : run history struct (H.meta, H.iter records)
%       runInfo : small metadata struct for external scripts/tools
%
% Typical User workflow
% ---------------------------
% 1) Generate or load a geometry:
%       “1 Random Config” or “Load Config”
% 2) Set FoA + thresholds (Array & Field Parameters), choose octants.
% 3) Click “Apply Changes” (after any edits) to synchronise all parameters.
% 4) Evaluate baseline FoA:
%       “2 Run WAH”
% 5) Optimise:
%       “3 Optimise”
% 6) Review results:
%       “Report”, “Export Run”, and/or “Export Workspace”
%
% Computational caution
% ---------------------
% Evaluation cost grows rapidly with:
%   - smaller grid spacing (denser FoA grid),
%   - larger FoA radius (more grid points),
%   - more microphones (more pairwise delays and solve complexity),
%   - longer synthesis tail padding (more samples, slower TDOA estimation).
% Ensure adequate CPU/RAM resources before running high-resolution setups.
% Primary development/testing was performed on Apple M1/M2 systems, with only
% limited testing on an Intel i5–class processor.
%
% Licensing and disclaimers
% -------------------------
% This software is provided “as is”, without warranty of any kind. Use at your
% own risk. The author assumes no liability for results, field use, or downstream
% decisions based on this tool.
%
% Author / Citation
% -----------------
% Ravi Umadi (2026)
% WAH-i —Optimising Array Geometry for Customised Localisation Accuracy
% Related work is referenced in the About dialog within the GUI.
%
% -------------------------------------------------------------------------
% Code: GUI state initialisation, defaults, theme, UI layout,
% callbacks, evaluation routines, optimiser, and report utilities.
% -------------------------------------------------------------------------

%% -------------------- State --------------------
IS_DEPLOYED = isdeployed;
S = struct();
S.M = [];
S.gridPts = [];
S.lastStats = struct('passRate',0,'meanErr',NaN,'p95Err',NaN,'maxErr',NaN,'threshold_cm',NaN);
S.diag = struct();
S.metric = "pos";
S.stopFlag = false;
S.rngSeed = 42;

% Quadrant selection state (all enabled by default)
S.quadrants = struct( ...
    'frontTopLeft',true,'frontTopRight',true,'frontBottomLeft',true,'frontBottomRight',true, ...
    'backTopLeft',true,'backTopRight',true,'backBottomLeft',true,'backBottomRight',true );

% Run history (batch-style)
S.H = [];
S.runInfo = [];
%LabelFontSize = 10;  % overridden by THEME below

%% -------------------- Defaults --------------------
DEF = struct();

% State defaults
DEF.rngSeed = 1;

% Array & Field defaults
DEF.Lmax        = 0.50;
DEF.N           = 4;
DEF.minPairDist = 0.125;
DEF.minDistAny  = 0.08;
DEF.R           = 2.0;
DEF.Rin         = 0.30;
DEF.step        = 0.50;
DEF.epsDist     = 0.03;
DEF.best_cm     = 10;
DEF.tol_cm      = 5;
DEF.lambdaMean  = 0.15;
DEF.failPen_cm  = 1000;

% Optimiser defaults
DEF.maxIters    = 25;
DEF.probes      = 14;
DEF.initStep    = 0.25;
DEF.minStep     = 0.03;
DEF.stall       = 12;
DEF.shake       = 0.08;
DEF.epsJ        = 1e-6;
DEF.targetPass  = 0.95;

% Call synthesis defaults
DEF.fs       = 192000;
DEF.dur_ms   = 5;
DEF.f1       = 30000;
DEF.f0       = 70000;
DEF.snr_db   = 40;
DEF.tail_pct = 200;
DEF.callType = 'FM';

% Quadrants defaults (all enabled)
DEF.quadrants = struct( ...
    'frontTopLeft',true,'frontTopRight',true,'frontBottomLeft',true,'frontBottomRight',true, ...
    'backTopLeft',true,'backTopRight',true,'backBottomLeft',true,'backBottomRight',true );

% -------------------- Modern theme --------------------
THEME = struct();
preferred = ["Lato", "Segoe UI", "Helvetica Neue", "Helvetica", "DejaVu Sans", "Arial"];

try
    if isdeployed
        % Deployed: be conservative (fonts can behave differently)
        available = strings(0,1);
    else
        available = string(listfonts);
    end
catch
    available = strings(0,1);
end

THEME.fontName = 'Helvetica'; % safe fallback
for f = preferred
    if any(strcmpi(available, f))
        THEME.fontName = char(f);
        break;
    end
end
% Fresher, more colourful palette (soft pastels + teal accent)
THEME.baseBg       = [0.98 0.995 0.99]; % very light mint background
THEME.panelBg      = [0.995 1.0 0.99];   % card background with a mint hint
THEME.mutedBg      = [0.98 0.98 1.00];   % soft lavender for muted panels
THEME.panelAlt     = [1.00 0.995 0.98];  % warm cream alternate card
THEME.textColor    = [0.10 0.12 0.15];
THEME.accentColor  = [0.00 0.62 0.60];   % teal-cyan accent
THEME.secondaryBtn = [0.88 0.94 0.97];   % soft blue-gray for secondary buttons
THEME.successColor = [0.15 0.75 0.45];   % brighter green for success/apply
THEME.headerTeal   = [0.62 0.82 0.78];   % neutral teal tint for headers
LabelFontSize = 11;

%% -------------------- UI --------------------
fig = uifigure('Name','WAH-i — Microphone Array Design & Optimisation (Umadi, 2026) V1.0','Position',[60 60 1650 920], 'Resize', 'off');
fig.Color = THEME.baseBg;
% Some MATLAB releases don't support setting FontName/FontSize on uifigure.
% Wrap in try/catch to remain compatible.
try
    fig.FontName = THEME.fontName;
    fig.FontSize = 11;
catch
    % uifigure doesn't support these properties in this MATLAB version — skip.
end

%% -------------------- Top menu --------------------
mHelp = uimenu(fig, 'Text', 'Menu');

uimenu(mHelp, 'Text', 'About', 'MenuSelectedFcn', @(~,~)showAbout());
uimenu(mHelp, 'Text', 'Tips',  'MenuSelectedFcn', @(~,~)showTips());
DOC_URL = 'https://raviumadi.github.io/WAHi-Optimiser';
uimenu(mHelp, 'Text', 'Guide', 'MenuSelectedFcn', @(~,~)openGuide(DOC_URL));

% Set control defaults so created UI components inherit the theme font where supported.
try
    fig.DefaultUilabelFontName    = THEME.fontName;
    fig.DefaultUilabelFontSize    = 11;
    fig.DefaultUibuttonFontName   = THEME.fontName;
    fig.DefaultUibuttonFontSize   = 10;
    fig.DefaultUieditfieldFontName= THEME.fontName;
    fig.DefaultUieditfieldFontSize= 10;
    fig.DefaultUidropdownFontName = THEME.fontName;
    fig.DefaultUidropdownFontSize = 10;
    fig.DefaultUipanelFontName    = THEME.fontName;
    fig.DefaultUipanelFontSize    = 11;
    fig.DefaultUicheckboxFontName = THEME.fontName;
catch
    % Older MATLAB releases may not support default properties on uifigure — ignore.
end

% 1 row, 2 columns: Controls | Plot
root = uigridlayout(fig,[1 2]);
root.ColumnWidth = {800, '1x'};
root.RowHeight   = {'1x'};
root.Padding = [12 12 12 12];
root.ColumnSpacing = 12;
% root.BackgroundColor = 'w';
%% Left panel (scrollable)
pLeft = uipanel(root,'Title','CONFIGURATIONS & CONTROLS','FontWeight','bold','FontSize',12);
pLeft.Layout.Row = 1;
pLeft.Layout.Column = 1;
pLeft.BackgroundColor = THEME.panelAlt;
% pLeft.Scrollable = 'on';

% Two columns inside the controls panel (3 rows × 2 cols)
% Left col:  Array/Field  | Optimiser  | Call synthesis
% Right col: Mic editor   | Display    | Workflow actions
gLeft = uigridlayout(pLeft,[3 2]);
gLeft.RowHeight      = {370, 250, '1x'};     % tune heights as you like
gLeft.ColumnWidth    = {'0.8x','1x'};
gLeft.Padding        = [10 10 10 10];
gLeft.RowSpacing     = 12;
gLeft.ColumnSpacing  = 12;
% gLeft.Scrollable = 'on';
%% Right panel: plot + stats + quadrants
pRight = uipanel(root,'Title','3D FIELD OF ACCURACY','FontWeight','bold','FontSize',12);
pRight.Layout.Row = 1; pRight.Layout.Column = 2;
pRight.BackgroundColor = THEME.panelBg;

gRight = uigridlayout(pRight,[2 1]);
gRight.RowHeight = {'1x', 90};   % plot | stats
gRight.Padding = [10 10 10 10];
gRight.RowSpacing = 10;

% --- Plot area (unchanged) ---
gPlot = uigridlayout(gRight,[2 1]);
gPlot.Layout.Row = 1;
gPlot.RowHeight = {22,'1x'};
gPlot.Padding = [0 0 0 0];
gPlot.RowSpacing = 6;

lblStatus = uilabel(gPlot,'Text','Ready','FontWeight','bold');
lblStatus.Layout.Row = 1;
lblStatus.FontName = THEME.fontName;
lblStatus.FontColor = THEME.textColor;
lblStatus.FontSize = 12;

ax = uiaxes(gPlot);
ax.Layout.Row = 2;
grid(ax,'on'); axis(ax,'equal'); view(ax,3);
xlabel(ax,'X (m)','FontWeight','bold','FontSize',10);
ylabel(ax,'Y (m)','FontWeight','bold','FontSize',10);
zlabel(ax,'Z (m)','FontWeight','bold','FontSize',10);

% --- Bottom row: stats only ---
pStats = uipanel(gRight,'Title','Performance','FontWeight','bold','FontSize',11);
pStats.Layout.Row = 2;
pStats.BackgroundColor = THEME.panelBg;

gStats = uigridlayout(pStats,[1 1]);
gStats.Padding = [10 8 10 10];

lblStats = uilabel(gStats);
% lblStats.FontName  = 'Helvetica';   % nice alignment (optional)
% lblStats.FontWeight = 'bold';
lblStats.FontSize  = 12;
lblStats.FontName  = THEME.fontName;
lblStats.FontColor = THEME.textColor;
lblStats.VerticalAlignment   = 'center';
lblStats.HorizontalAlignment = 'center';
lblStats.WordWrap = 'on';

% Initial placeholder
lblStats.Text = sprintf([ ...
    'Pass: --- |  ' ...
    'Mean: ---  |  ' ...
    'P95:  --- |  ' ...
    'Max:  ---\n\n' ...
    'Fails: ---  |  ' ...
    'Geometry: ---  | ' ...
    'Grid Points: ---' ]);
%% ==================== Array & Field Parameters (geometry + FoA/grid/threshold) ====================
pArray = uipanel(gLeft,'Title','Array & Field Parameters','FontWeight','bold','FontSize',11);
pArray.Layout.Row = 1;
pArray.Layout.Column = 1;
pArray.BackgroundColor = THEME.panelBg;
pArray.Scrollable = 'on';

% Rows: RNG(2 rows) + 10 params = 12 total
gArray = uigridlayout(pArray,[12 2]);
gArray.RowHeight      = repmat({20},1,12);
gArray.ColumnWidth    = {200,'1x'};  % label | control
gArray.Padding        = [10 10 10 10];
gArray.RowSpacing     = 6;
gArray.ColumnSpacing  = 8;
gArray.Scrollable = 'on';
r = 1;

% --- RNG seed ---
lbl = uilabel(gArray,'Text','Random seed:','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edSeed = uieditfield(gArray,'numeric','Value',DEF.rngSeed,'Limits',[0 inf], ...
    'RoundFractionalValues','on');
edSeed.Layout.Row = r; edSeed.Layout.Column = 2;
setTip(edSeed,"Initialises the random number generator used for random geometries and optimiser probes.");
r = r + 1;


edSeed.ValueChangedFcn = @(~,~)markDirty("Seed changed");

% --- Array radius limit (Lmax) ---
lbl = uilabel(gArray,'Text','Array radius limit (m):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edLmax = uieditfield(gArray,'numeric','Value',DEF.Lmax,'Limits',[0.05 5], ...
    'ValueDisplayFormat','%.3f');
edLmax.Layout.Row = r; edLmax.Layout.Column = 2;
setTip(edLmax,"Maximum allowed distance of any microphone from the array centre (arm length constraint).");
r = r + 1;

% --- Number of microphones ---
lbl = uilabel(gArray,'Text','Number of microphones:','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edN = uieditfield(gArray,'numeric','Value',DEF.N,'Limits',[4 24], ...
    'RoundFractionalValues','on');
edN.Layout.Row = r; edN.Layout.Column = 2;
setTip(edN,"Total number of microphones in the array (used for randomisation and validation).");
r = r + 1;

% --- Minimum mic spacing ---
lbl = uilabel(gArray,'Text','Minimum mic spacing (m):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edMinPair = uieditfield(gArray,'numeric','Value',DEF.minPairDist,'Limits',[0 5], ...
    'ValueDisplayFormat','%.3f');
edMinPair.Layout.Row = r; edMinPair.Layout.Column = 2;
setTip(edMinPair,"Minimum allowed distance between any microphone pair (prevents clustering).");
r = r + 1;

% --- Minimum source-to-any-mic distance ---
lbl = uilabel(gArray,'Text','Min source–mic distance (m):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edMinDistAny = uieditfield(gArray,'numeric','Value',DEF.minDistAny,'Limits',[0 5], ...
    'ValueDisplayFormat','%.3f');
edMinDistAny.Layout.Row = r; edMinDistAny.Layout.Column = 2;
setTip(edMinDistAny,"Grid points closer than this to any mic are excluded (avoids near-field/degenerate cases).");
r = r + 1;

% --- Field-of-Accuracy shell (outer radius) ---
lbl = uilabel(gArray,'Text','FoA outer radius (m):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edR = uieditfield(gArray,'numeric','Value',DEF.R,'Limits',[0.2 50], ...
    'ValueDisplayFormat','%.2f');
edR.Layout.Row = r; edR.Layout.Column = 2;
setTip(edR,"Outer radius of the 3D evaluation shell (where accuracy is evaluated).");
r = r + 1;

% --- Field-of-Accuracy inner radius ---
lbl = uilabel(gArray,'Text','FoA inner radius (m):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edRin = uieditfield(gArray,'numeric','Value',DEF.Rin,'Limits',[0 50], ...
    'ValueDisplayFormat','%.2f');
edRin.Layout.Row = r; edRin.Layout.Column = 2;
setTip(edRin,"Inner exclusion radius; evaluation points inside this radius are not used (hollow shell).");
r = r + 1;

% --- Grid spacing ---
lbl = uilabel(gArray,'Text','Grid spacing (m):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edStep = uieditfield(gArray,'numeric','Value',DEF.step ,'Limits',[0.05 5], ...
    'ValueDisplayFormat','%.2f');
edStep.Layout.Row = r; edStep.Layout.Column = 2;
setTip(edStep,"Spacing of the 3D evaluation grid (smaller = denser grid, slower evaluation).");
r = r + 1;

% --- Near-field epsilon distance ---
lbl = uilabel(gArray,'Text','Near-field epsilon (m):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edEps = uieditfield(gArray,'numeric','Value',DEF.epsDist ,'Limits',[0 5], ...
    'ValueDisplayFormat','%.3f');
edEps.Layout.Row = r; edEps.Layout.Column = 2;
setTip(edEps,"If a source point is within epsilon of any mic, a penalty is applied (prevents singularities).");
r = r + 1;

% --- Target error (best) ---
lbl = uilabel(gArray,'Text','Target error (cm):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edBest = uieditfield(gArray,'numeric','Value',DEF.best_cm,'Limits',[0.1 1e5], ...
    'ValueDisplayFormat','%.1f');
edBest.Layout.Row = r; edBest.Layout.Column = 2;
setTip(edBest,"Desired localisation error threshold (base) in centimetres.");
r = r + 1;

% --- Tolerance margin ---
lbl = uilabel(gArray,'Text','Tolerance margin (cm):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edTol = uieditfield(gArray,'numeric','Value',DEF.tol_cm,'Limits',[0 1e5], ...
    'ValueDisplayFormat','%.1f');
edTol.Layout.Row = r; edTol.Layout.Column = 2;
setTip(edTol,"Added margin for pass/fail threshold: threshold = Target error + Tolerance.");
r = r + 1;

% --- Mean-error weight (lambda) ---
lbl = uilabel(gArray,'Text','Mean-error weight (λ):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edLambda = uieditfield(gArray,'numeric','Value',DEF.lambdaMean,'Limits',[0 10], ...
    'ValueDisplayFormat','%.3f');
edLambda.Layout.Row = r; edLambda.Layout.Column = 2;
setTip(edLambda,"Weight of mean error term in the optimiser score J = passRate − λ·(meanErr/penalty).");
r = r + 1;

% --- Failure penalty ---
lbl = uilabel(gArray,'Text','Failure penalty (cm):','FontWeight','bold','FontSize', LabelFontSize);
lbl.Layout.Row = r; lbl.Layout.Column = 1;

edFailPen = uieditfield(gArray,'numeric','Value',DEF.failPen_cm,'Limits',[10 1e7], ...
    'ValueDisplayFormat','%.0f');
edFailPen.Layout.Row = r; edFailPen.Layout.Column = 2;
setTip(edFailPen,"Assigned error (cm) for failed/singular points so they count against the solution.");
% r ends at 12

edLmax.ValueChangedFcn      = @(~,~)onLmaxChanged();
edN.ValueChangedFcn         = @(~,~)markDirty("Mic count changed");
edMinPair.ValueChangedFcn   = @(~,~)markDirty("Array/Field changed");
edMinDistAny.ValueChangedFcn= @(~,~)markDirty("Array/Field changed");
edR.ValueChangedFcn         = @(~,~)markDirty("Array/Field changed");
edRin.ValueChangedFcn       = @(~,~)markDirty("Array/Field changed");
edStep.ValueChangedFcn      = @(~,~)markDirty("Array/Field changed");
edEps.ValueChangedFcn       = @(~,~)markDirty("Array/Field changed");
edBest.ValueChangedFcn      = @(~,~)markDirty("Threshold changed");
edTol.ValueChangedFcn       = @(~,~)markDirty("Threshold changed");
edLambda.ValueChangedFcn    = @(~,~)markDirty("Scoring changed");
edFailPen.ValueChangedFcn   = @(~,~)markDirty("Scoring changed");

%% ==================== Optimiser Parameters (single-column, unchanged variables) ====================
pOpt = uipanel(gLeft,'Title','Optimiser Parameters','FontWeight','bold','FontSize',LabelFontSize);
pOpt.Layout.Row = 2;
pOpt.Layout.Column = 1;
pOpt.BackgroundColor = THEME.mutedBg;


gOpt = uigridlayout(pOpt,[8 2]);
gOpt.RowHeight      = repmat({20},1,8);
gOpt.ColumnWidth    = {200,'1x'};
gOpt.Padding        = [10 10 10 10];
gOpt.RowSpacing     = 6;
gOpt.ColumnSpacing  = 8;
gOpt.Scrollable = 'on';
r2 = 1;

lbl = uilabel(gOpt,'Text','Max iterations:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edMaxIters = uieditfield(gOpt,'numeric','Value',DEF.maxIters,'Limits',[1 5000],'RoundFractionalValues','on');
edMaxIters.Layout.Row = r2; edMaxIters.Layout.Column = 2;
setTip(edMaxIters,"Maximum optimisation iterations (one mic perturbed per iteration).");
r2 = r2 + 1;

lbl = uilabel(gOpt,'Text','Probes per mic:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edProbes = uieditfield(gOpt,'numeric','Value',DEF.probes,'Limits',[1 64],'RoundFractionalValues','on');
edProbes.Layout.Row = r2; edProbes.Layout.Column = 2;
setTip(edProbes,"Number of direction probes tested for the selected mic per iteration.");
r2 = r2 + 1;

lbl = uilabel(gOpt,'Text','Initial step fraction:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edInitStep = uieditfield(gOpt,'numeric','Value',DEF.initStep,'Limits',[0 1],'ValueDisplayFormat','%.3f');
edInitStep.Layout.Row = r2; edInitStep.Layout.Column = 2;
setTip(edInitStep,"Initial probe step size as a fraction of Lmax (step = frac·Lmax).");
r2 = r2 + 1;

lbl = uilabel(gOpt,'Text','Minimum step fraction:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edMinStep = uieditfield(gOpt,'numeric','Value',DEF.minStep,'Limits',[0 1],'ValueDisplayFormat','%.3f');
edMinStep.Layout.Row = r2; edMinStep.Layout.Column = 2;
setTip(edMinStep,"Lower bound on step fraction (prevents steps becoming too small late in optimisation).");
r2 = r2 + 1;

lbl = uilabel(gOpt,'Text','Stall limit:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edStall = uieditfield(gOpt,'numeric','Value',DEF.stall,'Limits',[1 1e4],'RoundFractionalValues','on');
edStall.Layout.Row = r2; edStall.Layout.Column = 2;
setTip(edStall,"Number of consecutive non-improving iterations before a shake/reset is triggered.");
r2 = r2 + 1;

lbl = uilabel(gOpt,'Text','Shake fraction:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edShake = uieditfield(gOpt,'numeric','Value',DEF.shake,'Limits',[0 1],'ValueDisplayFormat','%.3f');
edShake.Layout.Row = r2; edShake.Layout.Column = 2;
setTip(edShake,"Shake magnitude as fraction of Lmax (adds jitter to escape local minima).");
r2 = r2 + 1;

lbl = uilabel(gOpt,'Text','Min improvement (ΔJ):','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edEpsJ = uieditfield(gOpt,'numeric','Value',DEF.epsJ,'Limits',[0 1],'ValueDisplayFormat','%.1e');
edEpsJ.Layout.Row = r2; edEpsJ.Layout.Column = 2;
setTip(edEpsJ,"Minimum score improvement required to accept a candidate geometry.");
r2 = r2 + 1;

lbl = uilabel(gOpt,'Text','Target pass rate:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = r2; lbl.Layout.Column = 1;
edTargetPass = uieditfield(gOpt,'numeric','Value',DEF.targetPass,'Limits',[0 1],'ValueDisplayFormat','%.2f');
edTargetPass.Layout.Row = r2; edTargetPass.Layout.Column = 2;
setTip(edTargetPass,"Stop early once pass rate reaches this threshold.");
% r2 ends at 8

edMaxIters.ValueChangedFcn  = @(~,~)markDirty("Optimiser changed");
edProbes.ValueChangedFcn    = @(~,~)markDirty("Optimiser changed");
edInitStep.ValueChangedFcn  = @(~,~)markDirty("Optimiser changed");
edMinStep.ValueChangedFcn   = @(~,~)markDirty("Optimiser changed");
edStall.ValueChangedFcn     = @(~,~)markDirty("Optimiser changed");
edShake.ValueChangedFcn     = @(~,~)markDirty("Optimiser changed");
edEpsJ.ValueChangedFcn      = @(~,~)markDirty("Optimiser changed");
edTargetPass.ValueChangedFcn= @(~,~)markDirty("Optimiser changed");

%% ==================== Call Synthesis Parameters (moved BELOW optimiser) ====================
pCall = uipanel(gLeft,'Title','Call Synthesis Parameters','FontWeight','bold','FontSize',LabelFontSize);
pCall.Layout.Row = 3;
pCall.Layout.Column = 1;
pCall.BackgroundColor = THEME.mutedBg;
pCall.Scrollable = 'on';

% 7 rows: fs, d, f0, f1, SNR, Tail, Call type
gCall = uigridlayout(pCall,[7 2]);
gCall.RowHeight      = repmat({20},1,7);
gCall.ColumnWidth    = {200,'1x'};
gCall.Padding        = [10 10 10 10];
gCall.RowSpacing     = 6;
gCall.ColumnSpacing  = 8;

rc = 1;

lbl = uilabel(gCall,'Text','Sample rate (Hz):','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = rc; lbl.Layout.Column = 1;
edFs = uieditfield(gCall,'numeric','Value',DEF.fs,'Limits',[8000 1e6],'RoundFractionalValues','on');
edFs.Layout.Row = rc; edFs.Layout.Column = 2;
setTip(edFs,"Sampling rate used to synthesise calls and simulate microphone signals.");
rc = rc + 1;

lbl = uilabel(gCall,'Text','Call duration (ms):','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = rc; lbl.Layout.Column = 1;
edDur = uieditfield(gCall,'numeric','Value',DEF.dur_ms,'Limits',[0.1 200],'ValueDisplayFormat','%.2f');
edDur.Layout.Row = rc; edDur.Layout.Column = 2;
setTip(edDur,"Call length in milliseconds (converted to seconds internally).");
rc = rc + 1;

lbl = uilabel(gCall,'Text','Start frequency f0 (Hz):','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = rc; lbl.Layout.Column = 1;
edF0 = uieditfield(gCall,'numeric','Value',DEF.f0,'Limits',[1 3e5],'RoundFractionalValues','on');
edF0.Layout.Row = rc; edF0.Layout.Column = 2;
setTip(edF0,"Higher frequency of the sweep (FM) or carrier for CF call type.");
rc = rc + 1;

lbl = uilabel(gCall,'Text','End frequency f1 (Hz):','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = rc; lbl.Layout.Column = 1;
edF1 = uieditfield(gCall,'numeric','Value',DEF.f1,'Limits',[1 3e5],'RoundFractionalValues','on');
edF1.Layout.Row = rc; edF1.Layout.Column = 2;
setTip(edF1,"Lower frequency of FM sweep (ignored for CF call type).");
rc = rc + 1;

lbl = uilabel(gCall,'Text','SNR (dB):','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = rc; lbl.Layout.Column = 1;
edSNR = uieditfield(gCall,'numeric','Value',DEF.snr_db,'Limits',[-30 100],'ValueDisplayFormat','%.1f');
edSNR.Layout.Row = rc; edSNR.Layout.Column = 2;
setTip(edSNR,"Signal-to-noise ratio injected into simulated microphone signals.");
rc = rc + 1;

lbl = uilabel(gCall,'Text','Envelope tail (%):','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = rc; lbl.Layout.Column = 1;
edTail = uieditfield(gCall,'numeric','Value',DEF.tail_pct,'Limits',[0 200],'ValueDisplayFormat','%.0f');
edTail.Layout.Row = rc; edTail.Layout.Column = 2;
setTip(edTail,"Percentage of duration used to pad the synthesised signal. Always ensure the maximum delay from a distant source covers simulated shifts. Larger tails increase computational cost.");
rc = rc + 1;

lbl = uilabel(gCall,'Text','Call type:','FontWeight','bold','FontSize',LabelFontSize);
lbl.Layout.Row = rc; lbl.Layout.Column = 1;
ddCallType = uidropdown(gCall,'Items',{'FM','CF'},'Value',DEF.callType);
ddCallType.Layout.Row = rc; ddCallType.Layout.Column = 2;
setTip(ddCallType,"FM = frequency-modulated sweep; CF = constant-frequency tone.");

edFs.ValueChangedFcn        = @(~,~)markDirty("Call synthesis changed");
edDur.ValueChangedFcn       = @(~,~)markDirty("Call synthesis changed");
edF0.ValueChangedFcn        = @(~,~)markDirty("Call synthesis changed");
edF1.ValueChangedFcn        = @(~,~)markDirty("Call synthesis changed");
edSNR.ValueChangedFcn       = @(~,~)markDirty("Call synthesis changed");
edTail.ValueChangedFcn      = @(~,~)markDirty("Call synthesis changed");
ddCallType.ValueChangedFcn  = @(~,~)markDirty("Call synthesis changed");

%% ==================== Workflow Actions ====================
pAct = uipanel(gLeft,'Title','Workflow Actions','FontWeight','bold','FontSize',LabelFontSize);
pAct.Layout.Row = 3;
pAct.BackgroundColor = THEME.panelAlt;
pAct.Layout.Row  = 3; pAct.Layout.Column  = 2;

% Outer grid to vertically centre content
gActOuter = uigridlayout(pAct,[3 1]);
gActOuter.RowHeight = {'fit','fit','fit'};   % top spacer | buttons | bottom spacer
gActOuter.ColumnWidth = {'1x'};
gActOuter.Padding = [10 10 10 10];

% Inner grid = your actual buttons
gAct = uigridlayout(gActOuter,[4 3]);
gAct.Layout.Row = 2;    % MIDDLE ROW → centred vertically
gAct.RowHeight = {34, 34, 34, 34};
gAct.ColumnWidth = {'1x','1x','1x'};
gAct.RowSpacing = 6;
gAct.ColumnSpacing = 8;

btnRand = uibutton(gAct,'Text','1 Random Config ♯','FontWeight','bold','FontSize',LabelFontSize,'ButtonPushedFcn',@(~,~)doRandomise());
btnRand.BackgroundColor = THEME.accentColor;
btnRand.FontColor = [1 1 1];
setTip(btnRand, "Generate a random array configuuration for N mics");
btnRand.Layout.Row = 1; btnRand.Layout.Column = 1;

btnLoad = uibutton(gAct,'Text','Load Config ⬆','FontWeight','bold','FontSize',10,'ButtonPushedFcn',@(~,~)doLoadGeometry());
btnLoad.BackgroundColor = THEME.secondaryBtn;
btnLoad.FontColor = THEME.textColor;
setTip(btnLoad, "Get microphone array configirations from a 3XN CSV file. Overrides Lmax and N");
btnLoad.Layout.Row = 1; btnLoad.Layout.Column = 2;

btnSaveCSV = uibutton(gAct,'Text','Save Config ⬇','FontSize',LabelFontSize,'ButtonPushedFcn',@(~,~)doSaveGeometry());
btnSaveCSV.BackgroundColor = THEME.secondaryBtn;
btnSaveCSV.FontColor = THEME.textColor;
setTip(btnSaveCSV, "Export the current micrphone geometry to a 3XN CSV file.")
btnSaveCSV.Layout.Row = 1; btnSaveCSV.Layout.Column = 3;

btnRun = uibutton(gAct,'Text','2 Run WAH ▶︎','FontWeight','bold','FontSize',LabelFontSize,'ButtonPushedFcn',@(~,~)doRunWAH());
btnRun.BackgroundColor = THEME.accentColor;
btnRun.FontColor = [1 1 1];
setTip(btnRun, "Run the Widefield Acoustics Heuristic algorithm. Estimates the FoA for the current array geometry.");
btnRun.Layout.Row = 2; btnRun.Layout.Column = [1 3];

btnOpt = uibutton(gAct,'Text','3 Optimise ᯓ★','FontWeight','bold','FontSize',LabelFontSize,'ButtonPushedFcn',@(~,~)doOptimiser());
btnOpt.BackgroundColor = THEME.accentColor;
btnOpt.FontColor = [1 1 1];
setTip(btnOpt, "Apply the WAH-i algorithm to current array geometry.");
btnOpt.Layout.Row = 3; btnOpt.Layout.Column = [1 3];

btnStop = uibutton(gAct,'Text','Stop ⏹','FontSize',LabelFontSize,'Enable','off','ButtonPushedFcn',@(~,~)doStop());
btnStop.BackgroundColor = THEME.secondaryBtn;
btnStop.FontColor = THEME.textColor;
setTip(btnStop, "Stop the optimiser. Finishes the run before stopping, takes time.");
btnStop.Layout.Row = 4; btnStop.Layout.Column = 1;

btnExportRun = uibutton(gAct,'Text','Export Run ➜]','FontSize',LabelFontSize,'ButtonPushedFcn',@(~,~)doExportRun(), ...
    'Enable','off');btnExportRun.BackgroundColor = THEME.secondaryBtn;
btnExportRun.FontColor = THEME.textColor;
setTip(btnExportRun, "Save the optimiser run log to a .mat file for further analyses via Matlab.");
btnExportRun.Layout.Row = 4; btnExportRun.Layout.Column = 2;

btnReport = uibutton(gAct,'Text','Report ℹ','FontSize',LabelFontSize,'ButtonPushedFcn',@(~,~)doReport(), ...
    'Enable','off');btnReport.BackgroundColor = THEME.secondaryBtn;
btnReport.FontColor = THEME.textColor;
setTip(btnReport, "See the report for the current optimisation session.");
btnReport.Layout.Row = 4; btnReport.Layout.Column = 3;

%% ==================== Mic Editor ====================
pMic = uipanel(gLeft,'Title','Manual Management','FontWeight','bold','FontSize',LabelFontSize);
pMic.BackgroundColor = THEME.panelBg;
pMic.Layout.Row  = 1;          % keep as you had
pMic.Layout.Column = 2;

% --- OUTER container: placement area + bottom buttons area ---
gMicOuter = uigridlayout(pMic,[2 1]);
gMicOuter.RowHeight = {'1x','fit'};   % top grows, bottom stays snug
gMicOuter.ColumnWidth = {'1x'};
gMicOuter.Padding = [10 10 10 10];
gMicOuter.RowSpacing = 10;

% ==================== (A) Placement area ====================
pMicPlace = uipanel(gMicOuter,'Title','Microphone Placement','FontWeight','bold','FontSize',LabelFontSize);
pMicPlace.Layout.Row = 1;
pMicPlace.BackgroundColor = THEME.panelAlt;

gMic = uigridlayout(pMicPlace,[6 4]);               % <-- ONLY the placement rows
gMic.RowHeight = {20, 22, 22, 22, 22, 18};          % Mic row + X/Y/Z + Tip + spare
gMic.ColumnWidth = {60, 120, '3x', '1x'};
gMic.Padding = [10 10 10 10];
gMic.RowSpacing = 18;
gMic.ColumnSpacing = 8;

rr = 1;

% Selected Mic + recenter
lbl = uilabel(gMic); lbl.Text = 'Mic:'; lbl.FontWeight = 'bold'; lbl.FontSize = 10;
lbl.Layout.Row = rr; lbl.Layout.Column = 1;

ddMic = uidropdown(gMic,'Items',compose('%d',1:edN.Value),'Value','1', ...
    'ValueChangedFcn',@(~,~)syncMicToControls());
ddMic.Layout.Row = rr; ddMic.Layout.Column = 2;

btnRecenter = uibutton(gMic,'Text','Re-centre','FontSize',9,'ButtonPushedFcn',@(~,~)recenterArray());
btnRecenter.Layout.Row = rr; btnRecenter.Layout.Column = [3 4];
btnRecenter.BackgroundColor = THEME.secondaryBtn;
btnRecenter.FontColor = THEME.textColor;
rr = rr + 1;

% X
lbl = uilabel(gMic); lbl.Text = 'X (m):'; lbl.FontWeight = 'bold'; lbl.FontSize = 10;
lbl.Layout.Row = rr; lbl.Layout.Column = 1;

edX = uieditfield(gMic,'numeric','Value',0,'ValueDisplayFormat','%.4f');
edX.Layout.Row = rr; edX.Layout.Column = 2;

slX = uislider(gMic,'Limits',[-edLmax.Value edLmax.Value],'Value',0);
slX.Layout.Row = rr; slX.Layout.Column = [3 4];
slX.FontSize  = 6;
slX.FontColor = THEME.accentColor;
updateSliderLimits();   % sets Limits
% Clear ticks/labels so sliders show no tick marks or labels
try
    slX.MajorTicks = [];
    slX.MajorTickLabels = {};
catch
end
drawnow;                % forces first paint

rr = rr + 1;

% Y
lbl = uilabel(gMic); lbl.Text = 'Y (m):'; lbl.FontWeight = 'bold'; lbl.FontSize = 10;
lbl.Layout.Row = rr; lbl.Layout.Column = 1;

edY = uieditfield(gMic,'numeric','Value',0,'ValueDisplayFormat','%.4f');
edY.Layout.Row = rr; edY.Layout.Column = 2;

slY = uislider(gMic,'Limits',[-edLmax.Value edLmax.Value],'Value',0);
slY.FontSize  = 6;
slY.FontColor = THEME.accentColor;
slY.Layout.Row = rr; slY.Layout.Column = [3 4];
updateSliderLimits();   % sets Limits
try
    slY.MajorTicks = [];
    slY.MajorTickLabels = {};
catch
end
drawnow;                % forces first paint
rr = rr + 1;

% Z
lbl = uilabel(gMic); lbl.Text = 'Z (m):'; lbl.FontWeight = 'bold'; lbl.FontSize = 10;
lbl.Layout.Row = rr; lbl.Layout.Column = 1;

edZ = uieditfield(gMic,'numeric','Value',0,'ValueDisplayFormat','%.4f');
edZ.Layout.Row = rr; edZ.Layout.Column = 2;

slZ = uislider(gMic,'Limits',[-edLmax.Value edLmax.Value],'Value',0);
slZ.Layout.Row = rr; slZ.Layout.Column = [3 4];
slZ.FontSize  = 6;
slZ.FontColor = THEME.accentColor;
updateSliderLimits();   % sets Limits
try
    slZ.MajorTicks = [];
    slZ.MajorTickLabels = {};
catch
end
drawnow;                % forces first paint
rr = rr + 1;


% ==================== (B) Apply / Reset area ====================
pMicApply = uipanel(gMicOuter,'Title','Value Changes','FontWeight','bold','FontSize',LabelFontSize);
pMicApply.Layout.Row = 2;
pMicApply.BackgroundColor = THEME.headerTeal;

gApply = uigridlayout(pMicApply,[2 1]);
gApply.RowHeight = {34, 34};
gApply.ColumnWidth = {'1x'};
gApply.Padding = [10 10 10 10];
gApply.RowSpacing = 8;

btnApplyAll = uibutton(gApply,'Text','Apply Changes', ...
    'FontWeight','bold','FontSize',LabelFontSize, ...
    'Enable','off', ...
    'ButtonPushedFcn',@(~,~)applyAllChanges());
btnApplyAll.Layout.Row = 1;
btnApplyAll.Layout.Column = 1;
btnApplyAll.BackgroundColor = THEME.accentColor;
btnApplyAll.FontColor = [1 1 1];
setTip(btnApplyAll,"Applies all edited settings (seed, N, field/grid params, octants, etc.) and refreshes the geometry display.");

btnReset = uibutton(gApply,'Text','Reset Defaults', ...
    'FontWeight','bold','FontSize',LabelFontSize, ...
    'ButtonPushedFcn',@(~,~)doResetDefaults());
btnReset.Layout.Row = 2;
btnReset.Layout.Column = 1;
btnReset.BackgroundColor = THEME.secondaryBtn;
btnReset.FontColor = THEME.textColor;
setTip(btnReset,"Resets all parameters to defaults, clears results/history, and generates a fresh default geometry.");

% Wire up callbacks (same as before)
slX.ValueChangingFcn = @(~,evt)micMoveFromSlider('x',evt.Value);
slY.ValueChangingFcn = @(~,evt)micMoveFromSlider('y',evt.Value);
slZ.ValueChangingFcn = @(~,evt)micMoveFromSlider('z',evt.Value);

edX.ValueChangedFcn = @(~,~)micMoveFromEdit();
edY.ValueChangedFcn = @(~,~)micMoveFromEdit();
edZ.ValueChangedFcn = @(~,~)micMoveFromEdit();

%% ==================== Display Options ====================
pDisp = uipanel(gLeft,'Title','Display & Handling','FontWeight','bold','FontSize',LabelFontSize);
pDisp.Layout.Row = 2;
pDisp.Layout.Column = 2;
pDisp.BackgroundColor = THEME.mutedBg;

% --- OUTER grid to vertically centre content ---
gDispOuter = uigridlayout(pDisp,[3 1]);
gDispOuter.RowHeight   = {'1x','fit','1x'};   % spacer | content | spacer
gDispOuter.ColumnWidth = {'1x'};
gDispOuter.Padding     = [10 10 10 10];

% ==================== ACTUAL content grid ====================
gDisp = uigridlayout(gDispOuter,[2 2]);
gDisp.Layout.Row = 2;          % <-- centre row
gDisp.Padding = [0 0 0 0];
gDisp.RowHeight = {28,'fit'};
gDisp.ColumnWidth = {'1x','1x'};
gDisp.RowSpacing = 10;
gDisp.ColumnSpacing = 8;

ddMetric = uidropdown(gDisp,'Items',{'Position Error (cm)','Angular Error (deg)'},...
    'Value','Position Error (cm)',...
    'ValueChangedFcn',@(~,~)metricChanged());
ddMetric.Layout.Row = 1; ddMetric.Layout.Column = 1;

btnExportWS = uibutton(gDisp,'Text','Export Workspace','ButtonPushedFcn',@(~,~)doExportWorkspace());

btnExportWS.Layout.Row = 1; btnExportWS.Layout.Column = 2;
btnExportWS.BackgroundColor = THEME.secondaryBtn;
btnExportWS.FontColor = THEME.textColor;
if IS_DEPLOYED
    btnExportWS.Enable = 'off';
    setTip(btnExportWS, "Feature only available via MATLAB launch");
else
    btnExportWS.Enable = 'on';
    setTip(btnExportWS, "Export the program variables to active MATLAB workspace.");
end

% ---- Octant selector panel INSIDE Display Options ----
pQuad = uipanel(gDisp,'Title','Optimisation Octants', 'BorderType','none', 'FontWeight','bold','FontSize',10);
setTip(pQuad, "Select which regions to include in the optimisation.");
pQuad.Layout.Row = 2;
pQuad.Layout.Column = [1 2];
% pQuad.BackgroundColor = [0.99 0.99 1];

gQuad = uigridlayout(pQuad,[4 2]);
gQuad.RowHeight = repmat({24},1,4);
gQuad.ColumnWidth = {'1x','1x'};
gQuad.Padding = [8 6 8 8];
gQuad.RowSpacing = 4;
gQuad.ColumnSpacing = 8;

cbFrontTopLeft     = mkCB(gQuad,1,1,'Front Top Left',true);
cbFrontTopRight    = mkCB(gQuad,1,2,'Front Top Right',true);
cbFrontBottomLeft  = mkCB(gQuad,2,1,'Front Bottom Left',true);
cbFrontBottomRight = mkCB(gQuad,2,2,'Front Bottom Right',true);
cbBackTopLeft      = mkCB(gQuad,3,1,'Back Top Left',true);
cbBackTopRight     = mkCB(gQuad,3,2,'Back Top Right',true);
cbBackBottomLeft   = mkCB(gQuad,4,1,'Back Bottom Left',true);
cbBackBottomRight  = mkCB(gQuad,4,2,'Back Bottom Right',true);

cbFrontTopLeft.ValueChangedFcn     = @(~,~)markDirty("Octants changed");
cbFrontTopRight.ValueChangedFcn    = @(~,~)markDirty("Octants changed");
cbFrontBottomLeft.ValueChangedFcn  = @(~,~)markDirty("Octants changed");
cbFrontBottomRight.ValueChangedFcn = @(~,~)markDirty("Octants changed");
cbBackTopLeft.ValueChangedFcn      = @(~,~)markDirty("Octants changed");
cbBackTopRight.ValueChangedFcn     = @(~,~)markDirty("Octants changed");
cbBackBottomLeft.ValueChangedFcn   = @(~,~)markDirty("Octants changed");
cbBackBottomRight.ValueChangedFcn  = @(~,~)markDirty("Octants changed");

%% ==================== Initialise ====================
% Ensure every UI component that supports FontName uses the theme font.
try
    uic = findall(fig,'-property','FontName');
    for ii = 1:numel(uic)
        try
            uic(ii).FontName = THEME.fontName;
        catch
            % ignore objects that for some reason can't be set
        end
    end
catch
    % older MATLAB versions may not allow findall on uifigure — ignore safely
end

applyAllChanges();  % sets seed + reconciles geometry for current N/Lmax
doRandomise();
%% ==================== Callbacks ====================

    function showAbout()

        html = strjoin([ ...
            "<!doctype html><html><head><meta charset='utf-8'>", ...
            "<style>", ...
            "body{font-family:system-ui,-apple-system,Segoe UI,Helvetica,Arial,sans-serif;", ...
            "font-size:14px;color:#111;line-height:1.35;margin:0;padding:14px;}", ...
            "h2{margin:0 0 10px 0;font-size:16px;}", ...
            ".muted{color:#444;}", ...
            ".code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;}", ...
            "ul{margin:6px 0 10px 18px;padding:0;}", ...
            "</style></head><body>", ...
            " " ...
            "<h2>ABOUT</h2>", ...
            "<div class='muted'>WAH-<i>i</i> — Widefield Acoustics Heuristic (inverse iterative)</div>", ...
            "<div>Version: 1.0</div>", ...
            "<div>Initial release: February 2026</div>", ...
            "<div>Copyright: Ravi Umadi - <span class='code'>https://biosonix.io</span></div>", ...
            "<div>License: GPLv3</div>", ...
            "<div>GitHub: <span class='code'>https://github.com/raviumadi/WAHi-Optimiser</span></div>", ...
            "<div>Documentation and User Guide <span class='code'>https://raviumadi.github.io/WAHi-Optimiser/</span></div>", ...
            "<br>", ...
            " " ...
            "<b>Related research articles</b><br>", ...
            "1. Umadi (2025), <b>Widefield acoustics heuristic: advancing microphone array design for accurate spatial tracking of echolocating bats</b>. <i>BMC Ecology and Evolution</i> — <span class='code'>https://doi.org/10.1186/s12862-025-02441-4</span><br>", ...
            ...
            " ", ...
            "2. Umadi (2026), <b>WAH-<i>i</i> : Optimising Microphone Array Geometry for Customised Localisation Accuracy</b>. <span class='code'>https://doi.org/10.64898/2026.02.07.704547</span><br>", ... % <-- UPDATE THIS
            " ", ...
            "3. Umadi (2025), <b>BATSY4-PRO: An Open-Source Multichannel Ultrasound Recorder for Field Bioacoustics</b>. <span class='code'>https://doi.org/10.1101/2025.08.11.669530</span><br>", ...
            " ", ...
            "4. Umadi (2026), <b>ESPERDYNE: A dual-band heterodyne monitor and ultrasound recorder for bioacoustic field surveys</b>. <i>Methods in Ecology and Evolution</i> - <span class='code'>https://doi.org/10.1111/2041-210x.70241</span><br>", ...
            " ", ...
            "<br>", ...
            "<div>WAH-<i>i</i> is an optimiser for finding microphone-array geometries that meet a desired field of accuracy.</div>", ...
            "<div>The interface exposes the full set of algorithm parameters, enabling unusually fine control over array-geometry trade-offs, ", ...
            "and making it practical to reach very high accuracy when combined with staged optimisation and careful refinement.</div>", ...
            "<br>", ...
            " " ...
            "<div><b>Reports and export options support a clean workflow:</b></div>", ...
            "<ul>", ...
            "<li>Export geometry to CSV for CAD pipelines</li>", ...
            "<li>Export run history to MAT for deeper analysis in MATLAB or Python</li>", ...
            "</ul>", ...
            "<br>", ...
            " " ...
            "<b>SUPPORT</b><br>", ...
            "<div>To help me continue releasing tools like this, and other related projects, freely and to promote science, consider buying me a coffee:</div>", ...
            "<div class='code'>https://buymeacoffee.com/raviumadi</div>", ...
            "<br>", ...
            "<div><i>If you can, help promote education in developing countries.</i></div>", ...
            "</body></html>" ...
            ], "");

        % Add as many as you want here (scales)
        linkLabels = { ...
            'Biosonix', ...
            'User Guide', ...
            'Widefield Acoustics Heuristic (Research Paper)', ...
            'WAH-i: Optimising Microphone Array Geometry (Research Paper)', ...
            'BATSY4-PRO (Research Paper)', ...
            'ESPERDYNE (Research Paper)', ...
            'Github Repo', ...
            'Buy me a coffee' ...
            };

        linkUrls = { ...
            'https://biosonix.io', ...
            'https://raviumadi.github.io/WAHi-Optimiser/', ...
            'https://doi.org/10.1186/s12862-025-02441-4', ...
            'https://doi.org/10.64898/2026.02.07.704547', ... % <-- UPDATE
            'https://doi.org/10.1101/2025.08.11.669530', ...
            'https://doi.org/10.1111/2041-210x.70241', ...
            'https://github.com/raviumadi/WAHi-Optimiser', ...
            'https://buymeacoffee.com/raviumadi' ...
            };

        showHtmlPopupWindow('About WAH-i', html, linkLabels, linkUrls);
    end
    function showHtmlPopupWindow(titleStr, htmlText, linkLabels, linkUrls)

        % Ensure scalar char (uihtml wants a single document string)
        if isstring(htmlText) && ~isscalar(htmlText)
            htmlText = strjoin(htmlText, "");
        end
        htmlText = char(htmlText);

        % Defensive defaults
        if nargin < 3 || isempty(linkLabels), linkLabels = {}; end
        if nargin < 4 || isempty(linkUrls),   linkUrls   = {}; end

        w = 760; h = 860;
        pos = fig.Position;
        x = pos(1) + (pos(3)-w)/2;
        y = pos(2) + (pos(4)-h)/2;

        d = uifigure('Name', titleStr, 'Position', [x y w h], ...
            'Color', THEME.panelBg, 'Resize','off', 'WindowStyle','modal');

        gl = uigridlayout(d,[3 1]);
        gl.RowHeight = {'1x', 40, 48};
        gl.Padding = [12 12 12 12];
        gl.RowSpacing = 10;
        
        htm = uihtml(gl);
        htm.Layout.Row = 1;
        htm.HTMLSource = htmlText;
            
        % --- Single link opener row ---
        row = uigridlayout(gl, [1 3]);
        row.Layout.Row = 2;
        row.ColumnWidth = {'1x', 120, 90};
        row.ColumnSpacing = 8;
        row.Padding = [0 0 0 0];
        row.BackgroundColor = THEME.panelBg;

        dd = uidropdown(row, 'Items', linkLabels, 'ItemsData', linkUrls);
        dd.Tooltip = 'Select a link to open in your browser';

        bOpen = uibutton(row, 'Text', 'Open link', 'ButtonPushedFcn', @(~,~)openUrl(dd.Value));
        bOpen.BackgroundColor = THEME.secondaryBtn;
        bOpen.FontColor = THEME.textColor;

        bCopy = uibutton(row, 'Text', 'Copy', 'ButtonPushedFcn', @(~,~)copyUrl(dd.Value));
        bCopy.BackgroundColor = THEME.secondaryBtn;
        bCopy.FontColor = THEME.textColor;

        if isempty(linkUrls)
            row.Visible = 'off';
        end

        % Close button
        btn = uibutton(gl, 'Text','Close', 'ButtonPushedFcn', @(~,~)delete(d));
        btn.Layout.Row = 3;
        btn.BackgroundColor = THEME.secondaryBtn;
        btn.FontColor = THEME.textColor;

        function openUrl(url)
            if isempty(url), return; end
            try
                web(char(url), '-browser');   % reliable: system browser
            catch ME
                uialert(d, sprintf('Could not open link:\n%s\n\n%s', char(url), ME.message), 'Link error');
            end
        end

        function copyUrl(url)
            if isempty(url), return; end
            try
                clipboard('copy', char(url));
            catch ME
                uialert(d, sprintf('Could not copy link:\n%s\n\n%s', char(url), ME.message), 'Clipboard error');
            end
        end
    end

    function showTips()
        txt = sprintf([ ...
            'Workflow Guidelines\n\n' ...
            '1) Multistep optimisation\n' ...
            'Begin with a coarse (low-resolution) grid to identify a strong initial array configuration. ' ...
            'Once a stable solution is reached, reuse that configuration as the starting point for a higher-resolution grid optimisation. ' ...
            'This staged approach is often more efficient and reliable than optimising at full resolution from the outset.\n\n' ...
            '2) Manual refinement\n' ...
            'Some initial configurations converge to local performance ceilings or become trapped in shallow optima. ' ...
            'In such cases, a full randomisation is not always necessary. ' ...
            'Small, deliberate manual adjustments to microphone positions can reshape the error distribution and enable the optimiser to escape these valleys.\n\n' ...
            '3) Building very high pass rates\n' ...
            'For stringent accuracy requirements, avoid targeting very high pass rates in a single run. ' ...
            'Instead, begin with a modest target (e.g. 60%% pass rate), then progressively tighten the accuracy threshold while reducing step size and adjusting penalty or weighting parameters. ' ...
            'Successive optimisation runs refine the geometry without undoing earlier gains, and typically yield more robust high-accuracy solutions.\n' ...
            ]);
        showPopupWindow('Optimisation Tips', txt);
    end

    function showPopupWindow(titleStr, bodyText)
        % Modal dialog-like uifigure (works reliably across versions)
        d = uifigure('Name', titleStr, ...
            'Position', centreOnParent(fig, [640 420]), ...
            'Color', THEME.panelBg, ...
            'Resize', 'off', ...
            'WindowStyle', 'modal');

        try
            d.FontName = THEME.fontName;
            d.FontSize = 11;
        catch
        end

        g = uigridlayout(d, [3 1]);
        g.RowHeight = {28, '1x', 40};
        g.Padding = [14 14 14 14];
        g.RowSpacing = 10;

        % Header
        lbl = uilabel(g, 'Text', titleStr, 'FontWeight', 'bold');
        lbl.FontColor = THEME.textColor;
        lbl.FontSize = 14;

        % Body (scrollable text)
        ta = uitextarea(g, ...
            'Value', splitlines(string(bodyText)), ...
            'Editable', 'off');
        ta.FontName = THEME.fontName;
        ta.FontSize = 11;

        % Close button
        btn = uibutton(g, 'Text', 'Close', 'ButtonPushedFcn', @(~,~)delete(d));
        btn.BackgroundColor = THEME.secondaryBtn;
        btn.FontColor = THEME.textColor;
    end

    function openGuide(url)
        try
            web(char(url), '-browser');   % system default browser
        catch ME
            uialert(fig, sprintf('Could not open documentation link:\n%s\n\n%s', char(url), ME.message), ...
                'Link error');
        end
    end

    function pos = centreOnParent(parentFig, sz)
        % sz = [w h], returns [x y w h] centred on parentFig
        w = sz(1); h = sz(2);
        try
            p = parentFig.Position; % [x y w h]
            x = p(1) + (p(3)-w)/2;
            y = p(2) + (p(4)-h)/2;
        catch
            % fallback centre-ish
            x = 200; y = 200;
        end
        pos = [round(x) round(y) w h];
    end

    % function updateQuadrants()
    %     S.quadrants.frontTopLeft     = cbFrontTopLeft.Value;
    %     S.quadrants.frontTopRight    = cbFrontTopRight.Value;
    %     S.quadrants.frontBottomLeft  = cbFrontBottomLeft.Value;
    %     S.quadrants.frontBottomRight = cbFrontBottomRight.Value;
    %     S.quadrants.backTopLeft      = cbBackTopLeft.Value;
    %     S.quadrants.backTopRight     = cbBackTopRight.Value;
    %     S.quadrants.backBottomLeft   = cbBackBottomLeft.Value;
    %     S.quadrants.backBottomRight  = cbBackBottomRight.Value;
    % 
    %     nEnabled = sum(structfun(@(x) x, S.quadrants));
    %     updateStatus(sprintf('Quadrants: %d/8 enabled', nEnabled));
    % end
    % 
    % function applySeed()
    %     S.rngSeed = edSeed.Value;
    %     rng(S.rngSeed,'twister');
    %     updateStatus(sprintf('RNG seed: %d', S.rngSeed));
    % end

    function metricChanged()
        if strcmp(ddMetric.Value,'Position Error (cm)')
            S.metric = "pos";
        else
            S.metric = "ang";
        end
        refreshPlot(false);
    end

    function doStop()
        S.stopFlag = true;
        updateStatus('Stopping optimiser...');
    end

    function doExportWorkspace()
    if isdeployed
        uialert(fig, 'Workspace export is only available in MATLAB (not in deployed apps).', ...
            'Not available');
        return;
    end

    assignin('base','iWAH_M',S.M);
    assignin('base','iWAH_gridPts',S.gridPts);
    assignin('base','iWAH_stats',S.lastStats);
    assignin('base','iWAH_diag',S.diag);
    assignin('base','iWAH_H',S.H);
    assignin('base','iWAH_runInfo',S.runInfo);
    updateStatus('Exported to workspace: iWAH_M, iWAH_gridPts, iWAH_stats, iWAH_diag, iWAH_H, iWAH_runInfo');
end

    function doExportRun()
        if isempty(S.H) || isempty(S.runInfo)
            uialert(fig,'No run history found. Run WAH or Optimise first.','Nothing to export');
            return;
        end
        try
            [file, path] = uiputfile({'*.mat','MAT-files (*.mat)'}, 'Save run as .mat', defaultRunFileName());
            if isequal(file,0)
                updateStatus('Export cancelled');
                return;
            end
            H = S.H;
            runInfo = S.runInfo;
            save(fullfile(path,file), 'H', 'runInfo', '-v7.3');
            updateStatus(sprintf('Saved run: %s', file));
        catch ME
            uialert(fig, ME.message, 'Export Error');
            updateStatus('Export failed');
        end
    end

    function doReport()
        if strcmp(btnReport.Enable,'off')
            uialert(fig,'Run the optimiser first to generate a report.','Report unavailable');
            return;
        end
        if isempty(S.H) || numel(S.H.iter) < 1
            uialert(fig,'No run history. Run WAH or Optimise first.','No data');
            return;
        end

        try
            f = iwah_make_single_run_report(S.H, "final", S.metric);
            f.Name = 'iWAH Report';
            movegui(f,'center');

            choice = uiconfirm(fig,'Save report figure as PNG?','Save report?', ...
                'Options',{'Yes','No'},'DefaultOption',2,'CancelOption',2);
            if strcmp(choice,'Yes')
                [file,path] = uiputfile({'*.png','PNG (*.png)'}, 'Save report', 'iwah_report.png');
                if ~isequal(file,0)
                    exportgraphics(f, fullfile(path,file), 'Resolution', 180);
                    updateStatus(sprintf('Report saved: %s', file));
                else
                    updateStatus('Report not saved');
                end
            else
                updateStatus('Report generated');
            end
        catch ME
            uialert(fig, ME.message, 'Report Error');
            updateStatus('Report failed');
        end
    end

    function doLoadGeometry()
        try
            [file, path] = uigetfile({'*.csv','CSV Files (*.csv)';'*.*','All Files'},...
                'Select Microphone Geometry CSV File');
            if isequal(file,0), updateStatus('Load cancelled'); return; end

            fullpath = fullfile(path, file);
            data = readmatrix(fullpath);

            if size(data,2) ~= 3
                uialert(fig, 'CSV must have exactly 3 columns (x, y, z)', 'Invalid Format');
                return;
            end
            if size(data,1) < 4
                uialert(fig, 'Need at least 4 microphones', 'Too Few Mics');
                return;
            end
            if size(data,1) > 24
                uialert(fig, 'Maximum 24 microphones supported', 'Too Many Mics');
                return;
            end
            if any(~isfinite(data(:)))
                uialert(fig, 'CSV contains invalid values (NaN or Inf)', 'Invalid Values');
                return;
            end

            S.M = data;
            edN.Value = size(S.M, 1);
            S.M = S.M - mean(S.M, 1);

            if ~isValidGeometry(S.M, edLmax.Value, edMinPair.Value)
                uialert(fig, sprintf(['Loaded geometry violates constraints.\n' ...
                    'Consider adjusting Lmax (%.3f m) or Min Pair Dist (%.3f m).'], ...
                    edLmax.Value, edMinPair.Value), ...
                    'Geometry Constraint Violation','Icon','warning');
            end

            resetResults();
            ddMic.Items = compose('%d',1:size(S.M,1));
            ddMic.Value = '1';
            updateSliderLimits();
            syncMicToControls();
            refreshPlot(true);
            updateStatsTable();
            updateStatus(sprintf('Loaded %d mics from %s', size(S.M,1), file));
        catch ME
            uialert(fig, ME.message, 'Load Error');
            updateStatus('Load failed');
        end
    end

    function doSaveGeometry()
        if isempty(S.M)
            uialert(fig, 'No geometry to save', 'No Geometry');
            return;
        end
        try
            [file, path] = uiputfile({'*.csv','CSV Files (*.csv)'},...
                'Save Microphone Geometry', 'mic_geometry.csv');
            if isequal(file,0), updateStatus('Save cancelled'); return; end
            writematrix(S.M, fullfile(path, file));
            updateStatus(sprintf('Saved %d mics to %s', size(S.M,1), file));
        catch ME
            uialert(fig, ME.message, 'Save Error');
            updateStatus('Save failed');
        end
    end

    function doRandomise()
        try
            N = round(edN.Value);
            ddMic.Items = compose('%d',1:N);
            ddMic.Value = '1';

            Lmax = edLmax.Value;
            minD = edMinPair.Value;

            S.M = randomArrayNonCoplanar(N, Lmax, minD);
            S.M = S.M - mean(S.M,1);

            resetResults();
            updateSliderLimits();
            updateSliderTicks();
            syncMicToControls();
            refreshPlot(true);
            updateStatsTable();
            updateStatus('Geometry randomised');
        catch ME
            uialert(fig, ME.message, 'Randomise Error');
            updateStatus('Randomise failed');
        end
    end

    function doRunWAH()
        if isempty(S.M)
            uialert(fig,'No geometry. Click Randomise or Load CSV first.','Missing Geometry');
            return;
        end
        updateStatus('Running WAH...');
        drawnow;

        try
            P = readParams();
            [loc, gridPts] = buildLocaliserAndGrid(P);

            [stats, diag] = evaluateFoA_full(loc, gridPts, P.threshold_cm, ...
                P.minDistToAnyMic, P.epsDist, P.failPenalty_cm, P.c);

            % Update state
            S.gridPts = gridPts;
            S.lastStats = stats;
            S.diag = diag;

            % Create history with baseline only
            initHistoryIfNeeded(P);
            logIter(0, S.M, NaN, NaN, [], P, readOPT(), P.quadrants, gridPts, stats, diag, scoreJ(stats, P.lambdaMean, P.failPenalty_cm), [], true, "WAH baseline");

            updateStatsTable();
            refreshPlot(false);
            updateStatus(sprintf('WAH done: %.1f%% pass', 100*stats.passRate));
        catch ME
            uialert(fig, ME.message, 'WAH Error');
            updateStatus('WAH failed');
        end
    end

    function doOptimiser()
        if isempty(S.M)
            uialert(fig,'No geometry. Click Randomise or Load CSV first.','Missing Geometry');
            return;
        end

        % Always start a fresh run history for optimiser runs
        S.H = [];
        S.runInfo = [];
        btnReport.Enable = 'off';      % disable until optimiser completes
        btnExportRun.Enable = 'off';   % disable until optimiser completes

        S.stopFlag = false;
        btnStop.Enable = 'on';
        btnOpt.Enable = 'off';
        btnRun.Enable = 'off';
        btnRand.Enable = 'off';
        btnLoad.Enable = 'off';

        updateStatus('Optimiser running...');
        drawnow;

        try
            P = readParams();
            OPT = readOPT();

            [locTemplate, gridPts] = buildLocaliserAndGrid(P);
            paramTemplate = locTemplate.param;

            initHistoryIfNeeded(P);
            S.H.meta.OPT = OPT; % keep latest OPT
            S.H.meta.Q = P.quadrants;

            % Baseline eval (fresh instance)
            paramBase = paramTemplate;
            paramBase.mic_positions = S.M;
            locBase = BatCallLocaliser(paramBase);

            [bestStats, bestDiag] = evaluateFoA_full(locBase, gridPts, P.threshold_cm, ...
                P.minDistToAnyMic, P.epsDist, P.failPenalty_cm, P.c);
            bestJ = scoreJ(bestStats, P.lambdaMean, P.failPenalty_cm);

            % Log baseline if history empty
            if isempty(S.H.iter)
                logIter(0, S.M, NaN, NaN, [], P, OPT, P.quadrants, gridPts, bestStats, bestDiag, bestJ, [], true, "baseline");
            end
            acceptUpdate(bestStats, bestDiag, gridPts);

            baseDirs = ruleSet2Dirs();

            stall = 0;
            for it = 1:OPT.maxIters
                if S.stopFlag, break; end

                Lmax = P.Lmax;
                minD = P.minPairDist;
                N = size(S.M,1);

                if OPT.maxIters > 1
                    stepFrac = max(OPT.minStepFrac, OPT.initStepFrac * (1 - (it-1)/(OPT.maxIters-1)));
                else
                    stepFrac = OPT.initStepFrac;
                end
                stepSize = stepFrac * Lmax;

                micIdx = mod(it-1, N) + 1;

                dirs = baseDirs(1:min(OPT.probesPerMic,size(baseDirs,1)),:);
                v = S.M(micIdx,:);
                nv = norm(v);
                radial = [1 0 0];
                if nv > 1e-9, radial = v / nv; end
                dirs = [dirs; radial; -radial]; %#ok<AGROW>

                bestLocalJ = -inf;
                bestLocalM = [];
                bestLocalStats = [];
                bestLocalDiag = [];

                for k = 1:size(dirs,1)
                    Mtry = S.M;
                    Mtry(micIdx,:) = Mtry(micIdx,:) + stepSize * dirs(k,:);
                    Mtry = Mtry - mean(Mtry,1);

                    r = vecnorm(Mtry,2,2);
                    if r(micIdx) > Lmax
                        Mtry(micIdx,:) = Mtry(micIdx,:) * (Lmax/r(micIdx)) * 0.999;
                        Mtry = Mtry - mean(Mtry,1);
                    end

                    if ~isValidGeometry(Mtry, Lmax, minD)
                        continue;
                    end

                    paramTry = paramTemplate;
                    paramTry.mic_positions = Mtry;
                    locTry = BatCallLocaliser(paramTry);

                    [stTry, dgTry] = evaluateFoA_full(locTry, gridPts, P.threshold_cm, ...
                        P.minDistToAnyMic, P.epsDist, P.failPenalty_cm, P.c);

                    Jtry = scoreJ(stTry, P.lambdaMean, P.failPenalty_cm);

                    if Jtry > bestLocalJ
                        bestLocalJ = Jtry;
                        bestLocalM = Mtry;
                        bestLocalStats = stTry;
                        bestLocalDiag = dgTry;
                    end
                end

                accepted = false;
                note = sprintf('rejected (stall=%d)', stall+1);

                if ~isempty(bestLocalM) && (bestLocalJ > bestJ + OPT.improveEpsJ)
                    S.M = bestLocalM;
                    bestJ = bestLocalJ;
                    bestStats = bestLocalStats;
                    bestDiag = bestLocalDiag;
                    stall = 0;
                    accepted = true;
                    note = "accepted (improved J)";
                    acceptUpdate(bestStats, bestDiag, gridPts);
                    updateStatus(sprintf('Iter %d/%d: pass=%.1f%% J=%.3f', it, OPT.maxIters, 100*bestStats.passRate, bestJ));
                else
                    stall = stall + 1;
                end

                bestLocal = struct('J',bestLocalJ,'M',bestLocalM,'stats',bestLocalStats,'diag',bestLocalDiag);
                logIter(it, S.M, micIdx, stepSize, dirs, P, OPT, P.quadrants, gridPts, bestStats, bestDiag, bestJ, bestLocal, accepted, note);

                if bestStats.passRate >= OPT.targetPass
                    break;
                end

                if stall >= OPT.stallLimit
                    Mshake = shakeGeometry(S.M, OPT.shakeFrac*Lmax, Lmax, minD);

                    if ~isValidGeometry(Mshake, Lmax, minD)
                        Mshake = randomArrayNonCoplanar(N, Lmax, minD);
                        Mshake = Mshake - mean(Mshake,1);
                    end

                    paramShake = paramTemplate;
                    paramShake.mic_positions = Mshake;
                    locShake = BatCallLocaliser(paramShake);

                    [stShake, dgShake] = evaluateFoA_full(locShake, gridPts, P.threshold_cm, ...
                        P.minDistToAnyMic, P.epsDist, P.failPenalty_cm, P.c);

                    Jshake = scoreJ(stShake, P.lambdaMean, P.failPenalty_cm);

                    S.M = Mshake;
                    bestStats = stShake;
                    bestDiag  = dgShake;
                    bestJ     = Jshake;
                    stall = 0;

                    acceptUpdate(bestStats, bestDiag, gridPts);
                    logIter(it+0.5, S.M, NaN, NaN, [], P, OPT, P.quadrants, gridPts, bestStats, bestDiag, bestJ, [], true, "shake/reset");
                    updateStatus(sprintf('Shake/reset: pass=%.1f%% J=%.3f', 100*bestStats.passRate, bestJ));
                end

                drawnow limitrate;
            end

            if S.stopFlag
                updateStatus('Optimiser stopped');
            else
                updateStatus(sprintf('Optimiser done: %.1f%% pass', 100*S.lastStats.passRate));
            end

            % Finalise history
            S.H.meta.timeEnd = datestr(now);
            S.H.meta.final = struct('M',S.M,'stats',S.lastStats,'diag',S.diag,'J',scoreJ(S.lastStats,P.lambdaMean,P.failPenalty_cm));

        catch ME
            uialert(fig, ME.message, 'Optimiser Error');
            updateStatus('Optimiser failed');
        end
        % Optimiser produced a valid history -> enable reporting/export
        if ~isempty(S.H) && isfield(S.H,'iter') && numel(S.H.iter) >= 1
            btnReport.Enable = 'on';
            btnExportRun.Enable = 'on';
        end

        btnStop.Enable = 'off';
        btnOpt.Enable  = 'on';
        btnRun.Enable  = 'on';
        btnRand.Enable = 'on';
        btnLoad.Enable = 'on';
    end

    function acceptUpdate(stats, diag, gridPts)
        updateSliderLimits();
        syncMicToControls();

        S.gridPts   = gridPts;
        S.lastStats = stats;
        S.diag      = diag;

        updateStatsTable();
        refreshPlot(false);
    end

    function updateSliderLimits()
        Lmax = edLmax.Value;

        slX.Limits = [-Lmax Lmax];
        slY.Limits = [-Lmax Lmax];
        slZ.Limits = [-Lmax Lmax];

        updateSliderTicks();   % keep ticks consistent with new limits
    end

    function onLmaxChanged()
        % When array radius changes we must immediately update slider limits
        % and ticks and ensure current slider values lie within the new limits.
        markDirty("Array/Field changed");
        try
            updateSliderLimits();
            updateSliderTicks();
            drawnow;
            % Clamp slider values to new limits
            try
                slX.Value = min(max(slX.Value, slX.Limits(1)), slX.Limits(2));
            catch
            end
            try
                slY.Value = min(max(slY.Value, slY.Limits(1)), slY.Limits(2));
            catch
            end
            try
                slZ.Value = min(max(slZ.Value, slZ.Limits(1)), slZ.Limits(2));
            catch
            end
        catch
            % ignore any compatibility issues
        end
    end

    function updateSliderTicks()
        Lmax = edLmax.Value;

        % Pick a nice mid tick; handle degenerate Lmax==0 case
        if ~isfinite(Lmax) || Lmax <= 0
            ticks = 0;
        else
            ticks = [-Lmax 0 Lmax];
            ticks = unique(ticks);
        end

        % Build labels (try string array first, then cell array fallback)
        try
            labsStr = compose('%.2f', ticks); % string array
        catch
            labsStr = [];
        end
        try
            labsCell = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
        catch
            labsCell = cellstr(num2str(ticks(:),'%.2f'))';
        end

        % Helper to assign ticks+labels robustly
        assignTicks = @(sld) tryAssign(sld, ticks, labsStr, labsCell);

        % Clear ticks/labels for the mic sliders - the lables were taking
        % space, the edit boxes update value, so not needed.
        try
            slX.MajorTicks = [];
            slX.MajorTickLabels = {};
        catch
        end
        try
            slY.MajorTicks = [];
            slY.MajorTickLabels = {};
        catch
        end
        try
            slZ.MajorTicks = [];
            slZ.MajorTickLabels = {};
        catch
        end
    end

    function syncMicToControls()
        if isempty(S.M), return; end
        i = str2double(ddMic.Value);
        p = S.M(i,:);
        slX.Value = p(1); slY.Value = p(2); slZ.Value = p(3);
        edX.Value = p(1); edY.Value = p(2); edZ.Value = p(3);

        % Keep WAH overlay if it exists
        hasOverlay = ~isempty(S.gridPts) && isfield(S.lastStats,'err_cm') && any(isfinite(S.lastStats.err_cm));
        refreshPlot(~hasOverlay);
    end

    function micMoveFromSlider(axisName, val)
        if isempty(S.M), return; end
        i = str2double(ddMic.Value);
        p = S.M(i,:);
        switch axisName
            case 'x', p(1) = val;
            case 'y', p(2) = val;
            case 'z', p(3) = val;
        end
        [ok,p2] = enforceMicConstraints(i,p);
        if ok
            S.M(i,:) = p2;
            edX.Value = p2(1); edY.Value = p2(2); edZ.Value = p2(3);
            refreshPlot(true);
        end
    end

    function micMoveFromEdit()
        if isempty(S.M), return; end
        i = str2double(ddMic.Value);
        p = [edX.Value, edY.Value, edZ.Value];
        [ok,p2] = enforceMicConstraints(i,p);
        if ok
            S.M(i,:) = p2;
            slX.Value = p2(1); slY.Value = p2(2); slZ.Value = p2(3);
            refreshPlot(true);
        end
    end

    function [ok,pOut] = enforceMicConstraints(i,pIn)
        ok = true;
        pOut = pIn;

        Lmax = edLmax.Value;
        minD = edMinPair.Value;

        r = norm(pOut);
        if r > Lmax
            pOut = pOut * (Lmax/r) * 0.999;
        end

        Mtry = S.M;
        Mtry(i,:) = pOut;
        Mtry = Mtry - mean(Mtry,1);

        if ~isValidGeometry(Mtry, Lmax, minD)
            ok = false;
            return;
        end

        S.M = Mtry;
        pOut = S.M(i,:);
    end

    function recenterArray()
        if isempty(S.M), return; end
        S.M = S.M - mean(S.M,1);
        syncMicToControls();
        updateStatus('Array re-centred');
    end

    function updateStatus(msg)
        lblStatus.Text = msg;
    end

    function refreshPlot(geometryOnly)
        cla(ax);
        hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax,3);

        Lmax = edLmax.Value;
        R = edR.Value;

        [sx,sy,sz] = sphere(30);
        surf(ax, Lmax*sx, Lmax*sy, Lmax*sz, 'FaceAlpha',0.05,'EdgeAlpha',0.12);
        surf(ax, R*sx,    R*sy,    R*sz,    'FaceAlpha',0.03,'EdgeAlpha',0.08);

        if ~isempty(S.M)
            cmap = lines(size(S.M,1));
            for ii = 1:size(S.M,1)
                scatter3(ax, S.M(ii,1),S.M(ii,2),S.M(ii,3), 120, cmap(ii,:), 'filled',...
                    'MarkerEdgeColor','k','LineWidth',1.5);
                text(ax, S.M(ii,1),S.M(ii,2),S.M(ii,3), sprintf('  M%d',ii), ...
                    'FontWeight','bold','FontSize',10,'Color',cmap(ii,:));
            end
        end
        plot3(ax,0,0,0,'pentagram','MarkerSize',12,'LineWidth',2,...
            'MarkerFaceColor','r','MarkerEdgeColor','k');

        if ~geometryOnly && ~isempty(S.gridPts) && isfield(S.lastStats,'err_cm')
            if S.metric == "pos"
                vals = S.lastStats.err_cm;
                ttl = 'Position Error (cm)';
            else
                vals = S.lastStats.ang_deg;
                ttl = 'Angular Error (deg)';
            end

            scatter3(ax, S.gridPts(:,1),S.gridPts(:,2),S.gridPts(:,3), ...
                28, vals, 'filled','MarkerFaceAlpha',0.65);

            colormap(ax, turbo);
            % --- colourbar ---
            cb = findall(ax, 'Type', 'ColorBar');    % look on the axes (not the figure)
            if isempty(cb) || ~isvalid(cb)
                cb = colorbar(ax);
            else
                cb = cb(1);                         % in case there are multiples
            end

            cb.Location      = 'eastoutside';
            cb.TickDirection = 'in';
            cb.Box           = 'off';
            cb.FontSize      = 10;

            % Make it thin
            if isprop(cb,'Thickness')               % newer MATLAB releases
                cb.Thickness = 6;                   % try 6–12
            else                                    % older: shrink Position width
                drawnow;                            % ensure cb.Position is initialised
                p = cb.Position;                    % [x y w h]
                p(3) = p(3) * 0.35;                 % thinner
                cb.Position = p;
            end

            cb.Label.String  = ttl;
            cb.Label.FontWeight = 'normal';
        end

        xlabel(ax,'X (m)','FontWeight','bold','FontSize',10);
        ylabel(ax,'Y (m)','FontWeight','bold','FontSize',10);
        zlabel(ax,'Z (m)','FontWeight','bold','FontSize',10);

        drawnow limitrate;
    end

    function updateStatsTable()
        % Geometry string
        if isempty(S.M)
            geomStr = '---';
        else
            Dmat = squareform(pdist(S.M));
            Dmat(Dmat==0) = inf;
            minPair = min(Dmat(:));
            r = vecnorm(S.M,2,2);
            geomStr = sprintf('%.3fm / %.3fm', minPair, max(r));
        end

        % Stats strings
        if isempty(S.lastStats) || ~isfield(S.lastStats,'meanErr') || isnan(S.lastStats.meanErr)
            passStr = '---'; meanStr = '---'; p95Str = '---'; maxStr = '---';
            usedStr = '---';
        else
            passStr = sprintf('%.1f%%', 100*S.lastStats.passRate);
            meanStr = sprintf('%.2f cm', S.lastStats.meanErr);
            p95Str  = sprintf('%.2f cm', S.lastStats.p95Err);
            maxStr  = sprintf('%.2f cm', S.lastStats.maxErr);
            if isfield(S.diag,'gridUsed')
                usedStr = sprintf('%d', S.diag.gridUsed);
            else
                usedStr = '---';
            end
        end

        if isempty(fieldnames(S.diag))
            failStr = '---';
        else
            failStr = sprintf('%d/%d/%d', S.diag.simFail, S.diag.tdoaFail, S.diag.solveFail);
        end

        lblStats.Text = sprintf([ ...
            'Pass:  %s  |  ' ...
            'Mean:  %s  |  ' ...
            'P95:  %s  |  ' ...
            'Max:  %s\n\n' ...
            'Fails:  %s (sim/tdoa/solve)  |  ' ...
            'Geometry:  %s (minPair/maxArm)  |  ' ...
            'Grid Points:  %s' ], ...
            passStr, meanStr, p95Str, maxStr, failStr, geomStr, usedStr);
    end

    function resetResults()
        S.gridPts = [];
        S.lastStats = struct('passRate',0,'meanErr',NaN,'p95Err',NaN,'maxErr',NaN,'threshold_cm',NaN);
        S.diag = struct();
        S.H = [];
        S.runInfo = [];
    end

    function P = readParams()
        P = struct();

        % ---- Array & Field ----
        P.Lmax = edLmax.Value;
        P.N = round(edN.Value);
        P.minPairDist = edMinPair.Value;
        P.R = edR.Value;
        P.Rin = edRin.Value;
        P.step = edStep.Value;

        P.best_cm = edBest.Value;
        P.tol_cm = edTol.Value;
        P.threshold_cm = P.best_cm + P.tol_cm;

        P.lambdaMean = edLambda.Value;
        P.failPenalty_cm = edFailPen.Value;

        P.minDistToAnyMic = edMinDistAny.Value;
        P.epsDist = edEps.Value;
        P.c = 343;

        % ---- Call synthesis ----
        P.fs = edFs.Value;
        P.d = edDur.Value/1000;      % ms -> s
        P.f0 = edF0.Value;
        P.f1 = edF1.Value;
        if P.f0 > P.f1
            tmp = P.f0; P.f0 = P.f1; P.f1 = tmp;
        end
        P.snr_db = edSNR.Value;
        P.tail = edTail.Value;
        P.callType = ddCallType.Value;

        % ---- Quadrants ----
        P.quadrants = S.quadrants;

        % ---- Validations ----
        if P.Rin >= P.R
            error('FoA inner radius must be smaller than FoA outer radius.');
        end
        if P.step <= 0
            error('Grid spacing must be > 0.');
        end
        if P.d <= 0
            error('Call duration must be > 0.');
        end
        if P.fs <= 0
            error('Sample rate must be > 0.');
        end
        if ~any(structfun(@(x)x, P.quadrants))
            error('No octants enabled.');
        end
    end

    function OPT = readOPT()
        OPT = struct();
        OPT.maxIters      = round(edMaxIters.Value);
        OPT.probesPerMic  = round(edProbes.Value);
        OPT.initStepFrac  = edInitStep.Value;
        OPT.minStepFrac   = edMinStep.Value;
        OPT.stallLimit    = round(edStall.Value);
        OPT.shakeFrac     = edShake.Value;
        OPT.improveEpsJ   = edEpsJ.Value;
        OPT.targetPass    = edTargetPass.Value;
    end

    function [loc, gridPts] = buildLocaliserAndGrid(P)
        if isempty(S.M), error('No mic geometry'); end

        gridPts = buildShellGrid(P.R, P.Rin, P.step, P.quadrants);
        if isempty(gridPts), error('Empty grid — check octant selection'); end

        param = struct();
        param.fs = P.fs;
        param.d = P.d;
        param.f0 = P.f0;
        param.f1 = P.f1;
        param.tail = P.tail;
        param.snr_db = P.snr_db;
        param.callType = P.callType;
        param.velocity = [0 0 0];
        param.mic_positions = S.M;

        loc = BatCallLocaliser(param);
    end

    function initHistoryIfNeeded(P)
        if ~isempty(S.H), return; end

        S.H = struct();
        S.H.meta = struct();
        S.H.meta.rngSeed = S.rngSeed;
        S.H.meta.timeStart = datestr(now);
        S.H.meta.P = P;
        S.H.meta.OPT = readOPT();
        S.H.meta.Q = P.quadrants;
        S.H.meta.start = struct('source',"GUI",'csvPath','', 'M', S.M);
        S.H.iter = struct('it',{},'M',{},'micIdx',{},'stepSize',{},'dirsTried',{}, ...
            'P',{},'OPT',{},'Q',{},'gridPts',{},'stats',{},'diag',{},'J',{}, ...
            'bestLocal',{},'accepted',{},'note',{});

        S.runInfo = struct();
        S.runInfo.seed = S.rngSeed;
        S.runInfo.N = size(S.M,1);
        S.runInfo.startNote = "GUI start";
        S.runInfo.csvPath = "";
        S.runInfo.mode = "gui";
        S.runInfo.configTag = "interactive";
        S.runInfo.rep = 1;
        S.runInfo.csvFile = "";
    end

    function logIter(it, M, micIdx, stepSize, dirs, P, OPT, Q, gridPts, stats, diag, J, bestLocal, accepted, note)
        rec = struct();
        rec.it = it;
        rec.M = M;
        rec.micIdx = micIdx;
        rec.stepSize = stepSize;
        rec.dirsTried = dirs;
        rec.P = P;
        rec.OPT = OPT;
        rec.Q = Q;
        rec.gridPts = gridPts;
        rec.stats = stats;
        rec.diag = diag;
        rec.J = J;
        rec.bestLocal = bestLocal;
        rec.accepted = accepted;
        rec.note = char(note);
        S.H.iter(end+1) = rec; %#ok<AGROW>
    end

    function fn = defaultRunFileName()
        N = size(S.M,1);
        ts = datestr(now,'yyyymmdd_HHMMSS');
        fn = sprintf('iwah_gui_run__%dmics__%s.mat', N, ts);
    end

    function pos = localiseTDOA_global(tdoa, mic_pos, c)
        % Solve: ||mi-x|| - ||m1-x|| = c*tdoa_i
        tdoa = tdoa(:);
        m1 = mic_pos(1,:);
        mi = mic_pos(2:end,:);

        fun = @(x) vecnorm(mi - x,2,2) - norm(m1 - x) - c*tdoa;

        x0 = mean(mic_pos,1) + [0.25 0 0];
        opts = optimoptions('lsqnonlin', ...
            'Display','off', ...
            'FunctionTolerance',1e-10, ...
            'StepTolerance',1e-10, ...
            'MaxIterations',200);

        pos = lsqnonlin(fun, x0, [], [], opts);
    end

% helper to set tooltip
    function setTip(h, tip)
        if isempty(h), return; end
        try
            if isvalid(h)
                h.Tooltip = char(tip);  % char() is safe across string/char inputs
            end
        catch
            % ignore environments where isvalid or Tooltip may error
        end
    end


    function markDirty(reason)
        if nargin < 1, reason = "Changes pending"; end
        S.dirty = true;

        if exist('btnApplyAll','var')
            try
                if isvalid(btnApplyAll)
                    btnApplyAll.Enable = 'on';
                end
            catch
            end
        end

        % Update status label directly (no dependency on updateStatus)
        if exist('lblStatus','var')
            try
                if isvalid(lblStatus)
                    lblStatus.Text = sprintf('%s — click Apply Changes', reason);
                end
            catch
            end
        end
    end

    function applyAllChanges()
        % Applies ALL pending edits in one shot:
        %  - octants -> S.quadrants
        %  - RNG seed
        %  - mic count reconcile (trim/extend)
        %  - refresh geometry-only plot + stats label
        %  - invalidate evaluation results

        % ---- local "safe UI" helpers (no dependency on updateStatus/fig scope) ----
        safeStatus = @(msg) setStatusSafe(msg);
        safeAlert  = @(msg,title) alertSafe(msg,title);

        try
            % 0) Guard: ensure S exists
            if ~exist('S','var')
                safeAlert('Internal state (S) not found. Ensure applyAllChanges is nested inside iwah_gui_full2().', ...
                    'Apply Changes Error');
                return;
            end

            % 1) Pull octant UI -> state
            S.quadrants.frontTopLeft     = cbFrontTopLeft.Value;
            S.quadrants.frontTopRight    = cbFrontTopRight.Value;
            S.quadrants.frontBottomLeft  = cbFrontBottomLeft.Value;
            S.quadrants.frontBottomRight = cbFrontBottomRight.Value;
            S.quadrants.backTopLeft      = cbBackTopLeft.Value;
            S.quadrants.backTopRight     = cbBackTopRight.Value;
            S.quadrants.backBottomLeft   = cbBackBottomLeft.Value;
            S.quadrants.backBottomRight  = cbBackBottomRight.Value;

            if ~any(structfun(@(x)x, S.quadrants))
                safeAlert('At least one octant must be enabled.', 'Invalid octant selection');
                return;
            end

            % 2) Apply RNG seed (no separate "apply seed" button needed)
            S.rngSeed = round(edSeed.Value);
            rng(S.rngSeed,'twister');

            % 3) Apply mic count changes (trim/extend geometry)
            Nnew = round(edN.Value);
            Lmax = edLmax.Value;
            minD = edMinPair.Value;

            if isempty(S.M)
                % If no geometry exists yet, initialise one
                S.M = randomArrayNonCoplanar(Nnew, Lmax, minD);
            else
                S.M = reconcileMicCount(S.M, Nnew, Lmax, minD);
            end

            % Keep centred
            S.M = S.M - mean(S.M,1);

            % 4) Update mic dropdown + sliders to new limits
            ddMic.Items = compose('%d', 1:size(S.M,1));
            if isempty(ddMic.Value) || ~any(strcmp(ddMic.Value, ddMic.Items))
                ddMic.Value = '1';
            end

            updateSliderLimits();
            syncMicToControls();  % includes geometry refresh

            % 5) Invalidate previous evaluation (grid+stats no longer valid)
            clearEvaluationResultsOnly();

            % 6) Update stats + geometry plot
            updateStatsTable();
            refreshPlot(true); % geometry + shells only (no colours)

            % 7) Clear dirty flag & disable button
            S.dirty = false;
            if exist('btnApplyAll','var') && isvalid(btnApplyAll)
                btnApplyAll.Enable = 'off';
            end

            nEnabled = sum(structfun(@(x)x, S.quadrants));
            safeStatus(sprintf('Applied changes. Run WAH to recompute.', ...
                S.rngSeed, size(S.M,1), nEnabled));

        catch ME
            safeAlert(ME.message, 'Apply Changes Error');
            safeStatus('Apply failed');
        end

        % ---------------- nested local helpers ----------------
        function setStatusSafe(msg)
            % Prefer existing updateStatus if available, else write lblStatus directly.
            try
                if exist('updateStatus','var') == 1
                    updateStatus(msg);
                    return;
                end
            catch
                % fall through
            end
            try
                if exist('lblStatus','var') && isvalid(lblStatus)
                    lblStatus.Text = msg;
                end
            catch
                % last resort
                warning('%s', msg);
            end
        end

        function alertSafe(msg, ttl)
            % Prefer fig if valid, else gcf, else just warn.
            try
                if exist('fig','var') && isvalid(fig)
                    uialert(fig, msg, ttl);
                    return;
                end
            catch
                % fall through
            end
            try
                uialert(gcf, msg, ttl);
            catch
                warning('%s: %s', ttl, msg);
            end
        end
    end


    function M = reconcileMicCount(M, Nnew, Lmax, minD)
        if isempty(M)
            M = randomArrayNonCoplanar(Nnew, Lmax, minD);
            return;
        end

        Nold = size(M,1);
        if Nnew == Nold
            return;
        elseif Nnew < Nold
            M = M(1:Nnew,:);
            % ensure still valid (if not, regenerate)
            if ~isValidGeometry(M, Lmax, minD)
                M = randomArrayNonCoplanar(Nnew, Lmax, minD);
            end
        else
            % Extend: keep existing, add new random points until valid
            Madd = randomArrayNonCoplanar(Nnew, Lmax, minD);
            % Keep your existing mics as much as possible (first Nold)
            Madd(1:Nold,:) = M;
            Madd = Madd - mean(Madd,1);
            if ~isValidGeometry(Madd, Lmax, minD)
                Madd = randomArrayNonCoplanar(Nnew, Lmax, minD);
            end
            M = Madd;
        end
    end

    function clearEvaluationResultsOnly()
        % Clear everything that depends on grid/params, but keep geometry
        S.gridPts = [];
        S.lastStats = struct('passRate',0,'meanErr',NaN,'p95Err',NaN,'maxErr',NaN,'threshold_cm',NaN);
        S.diag = struct();

        % history becomes inconsistent with new settings
        S.H = [];
        S.runInfo = [];
    end

    function doResetDefaults()
        choice = uiconfirm(fig, ...
            'Reset all parameters to default values and clear results/history?', ...
            'Reset defaults', ...
            'Options',{'Reset','Cancel'}, ...
            'DefaultOption',2, 'CancelOption',2);

        if ~strcmp(choice,'Reset')
            updateStatus('Reset cancelled');
            return;
        end

        % --- Set UI controls to defaults ---
        edSeed.Value      = DEF.rngSeed;

        edLmax.Value      = DEF.Lmax;
        edN.Value         = DEF.N;
        edMinPair.Value   = DEF.minPairDist;
        edMinDistAny.Value= DEF.minDistAny;
        edR.Value         = DEF.R;
        edRin.Value       = DEF.Rin;
        edStep.Value      = DEF.step;
        edEps.Value       = DEF.epsDist;
        edBest.Value      = DEF.best_cm;
        edTol.Value       = DEF.tol_cm;
        edLambda.Value    = DEF.lambdaMean;
        edFailPen.Value   = DEF.failPen_cm;

        edMaxIters.Value  = DEF.maxIters;
        edProbes.Value    = DEF.probes;
        edInitStep.Value  = DEF.initStep;
        edMinStep.Value   = DEF.minStep;
        edStall.Value     = DEF.stall;
        edShake.Value     = DEF.shake;
        edEpsJ.Value      = DEF.epsJ;
        edTargetPass.Value= DEF.targetPass;

        edFs.Value        = DEF.fs;
        edDur.Value       = DEF.dur_ms;
        edF0.Value        = DEF.f0;
        edF1.Value        = DEF.f1;
        edSNR.Value       = DEF.snr_db;
        edTail.Value      = DEF.tail_pct;
        ddCallType.Value  = DEF.callType;

        % Quadrants (checkboxes)
        cbFrontTopLeft.Value     = DEF.quadrants.frontTopLeft;
        cbFrontTopRight.Value    = DEF.quadrants.frontTopRight;
        cbFrontBottomLeft.Value  = DEF.quadrants.frontBottomLeft;
        cbFrontBottomRight.Value = DEF.quadrants.frontBottomRight;
        cbBackTopLeft.Value      = DEF.quadrants.backTopLeft;
        cbBackTopRight.Value     = DEF.quadrants.backTopRight;
        cbBackBottomLeft.Value   = DEF.quadrants.backBottomLeft;
        cbBackBottomRight.Value  = DEF.quadrants.backBottomRight;

        % --- Reset state ---
        S.rngSeed = DEF.rngSeed;
        rng(S.rngSeed,'twister');

        % Clear evaluation + history
        resetResults();           % wipes grid/stats/diag + H + runInfo
        btnReport.Enable    = 'off';
        btnExportRun.Enable = 'off';

        % Apply geometry-related updates + regenerate default geometry
        btnApplyAll.Enable = 'off';
        applyAllChanges();  % will reconcile N/Lmax etc and clear eval again safely
        doRandomise();      % fresh geometry under default constraints

        updateStatus('Defaults restored');
    end
    function [stats, diag] = evaluateFoA_full(loc, gridPts, threshold_cm, minDistToAnyMic, epsDist, failPenalty_cm, c)
        % Evaluates FoA using your simulate->TDOA->solve pipeline.
        nP = size(gridPts,1);

        err_cm  = nan(nP,1);
        ang_deg = nan(nP,1);

        simFail   = 0;
        tdoaFail  = 0;
        solveFail = 0;
        tooClose  = 0;

        mic_pos = loc.mic_positions;
        fs = loc.param.fs;

        usedMask = false(nP,1);

        for i = 1:nP
            src = gridPts(i,:);
            dToMics = vecnorm(mic_pos - src,2,2);

            if min(dToMics) < minDistToAnyMic
                tooClose = tooClose + 1;
                continue;
            end
            usedMask(i) = true;

            if min(dToMics) < epsDist
                err_cm(i)  = failPenalty_cm;
                ang_deg(i) = 180;
                continue;
            end

            try
                simres = loc.simulate(src);
            catch
                simFail = simFail + 1;
                continue;
            end

            try
                tdoa = estimateTDOA_xcorr(simres.signals, fs);
                if any(~isfinite(tdoa))
                    tdoaFail = tdoaFail + 1;
                    continue;
                end
            catch
                tdoaFail = tdoaFail + 1;
                continue;
            end

            try
                est = localiseTDOA_global(tdoa, mic_pos, c);
                if any(~isfinite(est))
                    solveFail = solveFail + 1;
                    continue;
                end

                e = norm(est - src) * 100; % m -> cm
                if ~isfinite(e)
                    solveFail = solveFail + 1;
                    continue;
                end
                err_cm(i) = e;

                % crude angular error (az/el from mic 1)
                vT = src - mic_pos(1,:);
                vE = est - mic_pos(1,:);
                azT = atan2d(vT(2),vT(1));
                elT = asind(vT(3)/max(norm(vT),eps));
                azE = atan2d(vE(2),vE(1));
                elE = asind(vE(3)/max(norm(vE),eps));

                dAz = mod((azE-azT)+180,360)-180;
                dEl = elE-elT;
                ang_deg(i) = hypot(dAz,dEl);
            catch
                solveFail = solveFail + 1;
                continue;
            end

            if mod(i,25)==0
                drawnow limitrate;
            end
        end

        % Fill NaNs on used points with penalty (so pass mask behaves sensibly)
        err_f = err_cm;
        err_f(isnan(err_f) & usedMask) = failPenalty_cm;

        ang_f = ang_deg;
        ang_f(isnan(ang_f) & usedMask) = 180;

        passMask = (err_f <= threshold_cm) & usedMask;
        passRate = sum(passMask) / max(1,sum(usedMask));

        usedErrs = err_f(usedMask);
        meanErr = mean(usedErrs);
        p95Err  = prctile(usedErrs,95);
        maxErr  = max(usedErrs);

        stats = struct();
        stats.err_cm = err_f;
        stats.ang_deg = ang_f;
        stats.passMask = passMask;
        stats.passRate = passRate;
        stats.meanErr = meanErr;
        stats.p95Err = p95Err;
        stats.maxErr = maxErr;
        stats.threshold_cm = threshold_cm;
        stats.usedMask = usedMask; % IMPORTANT: report/plotters rely on this

        diag = struct();
        diag.nP = nP;
        diag.gridUsed = sum(usedMask);
        diag.tooClose = tooClose;
        diag.simFail = simFail;
        diag.tdoaFail = tdoaFail;
        diag.solveFail = solveFail;
    end
end

%% ===================== REPORT FUNCTION (standalone figure) =====================
function fig = iwah_make_single_run_report(H, AFTER, METRIC)
% iWAH single-run report (enhanced)
% Adds:
%   - Better global BEFORE vs AFTER distributions (valid-only + catastrophic annotation)
%   - Accuracy vs radius with median + 10–90% band (for BOTH pos and ang)
%   - Mic drift plots (dx,dy,dz vs iteration; per mic)
%   - Mic start/end table (uitable)
%
% Inputs:
%   H     : run history struct (H.meta.P, H.iter)
%   AFTER : "final" | "best"
%   METRIC: "pos" | "ang"  (used for the 3D colour tiles; other tiles show both)

if nargin < 2, AFTER  = "final"; end
if nargin < 3, METRIC = "pos";   end

P = H.meta.P;
if ~isfield(P,'threshold_cm')
    P.threshold_cm = P.best_cm + P.tol_cm;
end

% ---- pick baseline and after ----
iBase  = 1;
iFinal = numel(H.iter);

passAll = arrayfun(@(x) x.stats.passRate, H.iter(:))';
p95All  = arrayfun(@(x) x.stats.p95Err,    H.iter(:))';
JAll    = arrayfun(@(x) x.J,              H.iter(:))';
iBest   = maxPassTieBreak(passAll, p95All, JAll);

switch string(AFTER)
    case "best",  iAfter = iBest;
    otherwise,    iAfter = iFinal;
end

base = H.iter(iBase);
aft  = H.iter(iAfter);

% ---- grid ----
gridPts = base.gridPts;
if isempty(gridPts)
    error('Report: base.gridPts is empty. Run WAH/Optimiser with grid logging enabled.');
end

% ---- masks ----
maskBase = true(size(gridPts,1),1);
maskAft  = true(size(gridPts,1),1);
if isfield(base.stats,'usedMask') && ~isempty(base.stats.usedMask), maskBase = base.stats.usedMask(:); end
if isfield(aft.stats,'usedMask')  && ~isempty(aft.stats.usedMask),  maskAft  = aft.stats.usedMask(:);  end

% ---- metrics ----
[posBase, ~] = pickMetric(base.stats, "pos");
[posAft,  ~] = pickMetric(aft.stats,  "pos");
[angBase, ~] = pickMetric(base.stats, "ang");
[angAft,  ~] = pickMetric(aft.stats,  "ang");

thr_cm = P.threshold_cm;            % for position thresholding
thr2   = 4 * thr_cm;                % "valid region" cut-off used in your script
cmax3d = 2 * thr_cm;                % colour cap for 3D tiles

% ---- pass trajectory extras (accepted / shake) if present ----
itVals = arrayfun(@(x) x.it, H.iter(:))';
accepted = false(numel(H.iter),1);
note = strings(numel(H.iter),1);
for i = 1:numel(H.iter)
    if isfield(H.iter(i),'accepted'), accepted(i) = logical(H.iter(i).accepted); end
    if isfield(H.iter(i),'note') && ~isempty(H.iter(i).note), note(i) = string(H.iter(i).note); end
end
isShake = contains(note,"shake/reset");

% ---- drift over iters ----
[M0, dM, dx, dy, dz] = micDriftXYZ(H); % baseline + per-iter deltas

% ---- table figure (mic start/end) ----
makeMicStartEndTableFigure(base.M, aft.M, AFTER);

% ---- main figure layout ----
fig = figure('Name','iWAH Report (enhanced)','Color',[1 1 1], 'Position',[40 40 1760 1040]);
t = tiledlayout(fig, 3, 4, 'Padding','compact', 'TileSpacing','compact');

% ============================================================
% (1) 3D BEFORE (chosen METRIC)
% ============================================================
ax1 = nexttile(t,1);
if string(METRIC) == "ang"
    valsB = angBase; lab = 'Angular error (deg)';
else
    valsB = posBase; lab = 'Position error (cm)';
end
plot3D(ax1, P, base.M, gridPts, valsB, maskBase, lab, thr_cm, cmax3d);
title(ax1, sprintf('BEFORE (it=%.1f): pass=%.1f%%, mean=%.2f cm, p95=%.2f cm', ...
    base.it, 100*base.stats.passRate, base.stats.meanErr, base.stats.p95Err));

% ============================================================
% (2) 3D AFTER (chosen METRIC)
% ============================================================
ax2 = nexttile(t,2);
if string(METRIC) == "ang"
    valsA = angAft; lab = 'Angular error (deg)';
else
    valsA = posAft; lab = 'Position error (cm)';
end
plot3D(ax2, P, aft.M, gridPts, valsA, maskAft, lab, thr_cm, cmax3d);
title(ax2, sprintf('AFTER (%s, it=%.1f): pass=%.1f%%, mean=%.2f cm, p95=%.2f cm', ...
    AFTER, aft.it, 100*aft.stats.passRate, aft.stats.meanErr, aft.stats.p95Err));

% ============================================================
% (3) Pass trajectory
% ============================================================
ax3 = nexttile(t,3);
grid(ax3,'on'); hold(ax3,'on');
plot(ax3, itVals, 100*passAll, 'LineWidth', 1.5);
if any(accepted), scatter(ax3, itVals(accepted), 100*passAll(accepted), 40, 'filled'); end
if any(isShake),  scatter(ax3, itVals(isShake),  100*passAll(isShake),  55, 'd');      end
xlabel(ax3,'Iteration (it)'); ylabel(ax3,'Pass rate (%)');
title(ax3,'Pass-rate progression');
legendItems = {'pass'};
if any(accepted), legendItems{end+1} = 'accepted'; end %#ok<AGROW>
if any(isShake),  legendItems{end+1} = 'shake/reset'; end %#ok<AGROW>
legend(ax3, legendItems, 'Location','best');
hold(ax3,'off');

% ============================================================
% (4) Mic geometry: before vs after
% ============================================================
ax4 = nexttile(t,4);
grid(ax4,'on'); axis(ax4,'equal'); hold(ax4,'on'); view(ax4,3);
xlabel(ax4,'X (m)'); ylabel(ax4,'Y (m)'); zlabel(ax4,'Z (m)');
title(ax4,'Mic geometry: before (dull) vs after (bright)');
beforeCol = repmat([0.55 0.55 0.55], size(base.M,1), 1);
afterCol  = lines(size(aft.M,1));
scatter3(ax4, base.M(:,1), base.M(:,2), base.M(:,3), 120, beforeCol, 'filled', ...
    'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth',1.2);
scatter3(ax4, aft.M(:,1),  aft.M(:,2),  aft.M(:,3),  120, afterCol, 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.4);
legend(ax4, {'before','after'}, 'Location','best');
hold(ax4,'off');

% ============================================================
% (5) Global distribution (POS) — valid region only + catastrophic annotation
% ============================================================
ax5 = nexttile(t,5);
plotGlobalDistributionValidOnly(ax5, posBase(maskBase), posAft(maskAft), thr_cm, 'Position error (cm)');

% ============================================================
% (6) Global distribution (ANG)
% ============================================================

ax6 = nexttile(t,6);
thr_ang = deriveAngularThreshold(angBase(maskBase));
plotGlobalDistributionValidOnly(ax6, angBase(maskBase), angAft(maskAft), thr_ang, 'Angular error (deg)');

% ============================================================
% (7) Accuracy vs distance (POS): points + median + 10–90% band
% ============================================================
ax7 = nexttile(t,7);
plotAccuracyVsRadius(ax7, gridPts, maskBase, maskAft, posBase, posAft, ...
    P, thr_cm, 'Position error (cm)');

% ============================================================
% (8) Accuracy vs distance (ANG): points + median + 10–90% band
% ============================================================
ax8 = nexttile(t,8);
plotAccuracyVsRadius(ax8, gridPts, maskBase, maskAft, angBase, angAft, ...
    P, thr_ang, 'Angular error (deg)');

% ============================================================
% (9)-(11) Mic drift: dx, dy, dz per mic across iteration
% ============================================================
ax9  = nexttile(t,9);
ax10 = nexttile(t,10);
ax11 = nexttile(t,11);
plotMicDriftAxes(ax9, ax10, ax11, itVals, dx, dy, dz, accepted, isShake);

% ============================================================
% (12) Summary text tile
% ============================================================
ax12 = nexttile(t,12);
axis(ax12,'off');
txt = sprintf([ ...
    'Run summary\n' ...
    '-----------\n' ...
    'N mics: %d\n' ...
    'Grid used (base/after): %d / %d\n' ...
    'Pass (base/after): %.1f%% / %.1f%%\n' ...
    'Pos mean (base/after): %.2f / %.2f cm\n' ...
    'Pos p95  (base/after): %.2f / %.2f cm\n' ...
    '\n' ...
    'Mic drift (after - before)\n' ...
    '  median |Δx|: %.4f m\n' ...
    '  median |Δy|: %.4f m\n' ...
    '  median |Δz|: %.4f m\n'], ...
    size(base.M,1), ...
    sum(maskBase), sum(maskAft), ...
    100*base.stats.passRate, 100*aft.stats.passRate, ...
    base.stats.meanErr, aft.stats.meanErr, ...
    base.stats.p95Err,  aft.stats.p95Err, ...
    median(abs(aft.M(:,1)-base.M(:,1))), ...
    median(abs(aft.M(:,2)-base.M(:,2))), ...
    median(abs(aft.M(:,3)-base.M(:,3))) ...
    );
text(ax12, 0, 1, txt, 'Units','normalized', 'VerticalAlignment','top', ...
    'FontName','Consolas', 'FontSize',10);

% Attach useful things for callers
fig.UserData = struct();
fig.UserData.base = base;
fig.UserData.after = aft;
fig.UserData.P = P;
fig.UserData.M0 = M0;
fig.UserData.dM = dM;

end

%% ===================== Plot helpers =====================

function plotGlobalDistributionValidOnly(ax, vB_all, vA_all, thr, xlab)
% Replicates your "valid region only" overlay + catastrophic annotation
grid(ax,'on'); hold(ax,'on');

thr2 = 4 * thr;

vB_all = vB_all(isfinite(vB_all));
vA_all = vA_all(isfinite(vA_all));

vB_valid = vB_all(vB_all <= thr2);
vA_valid = vA_all(vA_all <= thr2);

nB_cat = sum(vB_all > thr2);
nA_cat = sum(vA_all > thr2);

fracB_cat = nB_cat / max(1,numel(vB_all));
fracA_cat = nA_cat / max(1,numel(vA_all));

edges = linspace(0, thr2, 40);

% Before = grey, After = blue (nice “manuscript-like” contrast)
histogram(ax, vB_valid, edges, 'Normalization','probability', ...
    'FaceAlpha',0.35,'EdgeAlpha',0.15);
histogram(ax, vA_valid, edges, 'Normalization','probability', ...
    'FaceAlpha',0.35,'EdgeAlpha',0.15);

xline(ax, thr,  '--k','threshold');
xline(ax, thr2, ':k','4×threshold');

xlabel(ax, xlab);
ylabel(ax, 'Probability');
title(ax, sprintf('%s (≤ 4×thr = %.2f)', xlab, thr2));

legend(ax, {'before','after','threshold','4×threshold'}, 'Location','best');

txt = sprintf(['Catastrophic failures (>4×thr)\n' ...
    'before: %.1f%% (%d/%d)\n' ...
    'after:  %.1f%% (%d/%d)'], ...
    100*fracB_cat, nB_cat, numel(vB_all), ...
    100*fracA_cat, nA_cat, numel(vA_all));

text(ax, 0.98, 0.98, txt, 'Units','normalized', ...
    'HorizontalAlignment','right', 'VerticalAlignment','top', ...
    'FontSize',9, 'BackgroundColor',[1 1 1 0.85], 'EdgeColor',[0.7 0.7 0.7]);

hold(ax,'off');
end

function plotAccuracyVsRadius(ax, gridPts, maskBase, maskAft, vBase, vAft, P, thr, ylab)
hold(ax,'on'); grid(ax,'on');

rBase = vecnorm(gridPts(maskBase,:),2,2);
rAft  = vecnorm(gridPts(maskAft,:),2,2);

eBase = vBase(maskBase);
eAft  = vAft(maskAft);

% cap for visibility only (do not hide catastrophes; just clip display)
eCap = 2 * thr;
eBaseP = min(eBase, eCap);
eAftP  = min(eAft,  eCap);

scatter(ax, rBase, eBaseP, 14, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha',0.25);
scatter(ax, rAft,  eAftP,  14, [0.1 0.5 0.9], 'filled', 'MarkerFaceAlpha',0.35);

% distance bins for summary curves
if isfield(P,'Rin') && isfield(P,'R') && isfinite(P.Rin) && isfinite(P.R) && P.R > P.Rin
    rbins = linspace(P.Rin, P.R, 10);
else
    rr = [rBase(:); rAft(:)];
    rbins = linspace(min(rr), max(rr), 10);
end

[rmedB, rloB, rhiB] = radialQuantiles(rBase, eBase, rbins);
[rmedA, rloA, rhiA] = radialQuantiles(rAft,  eAft,  rbins);

plot(ax, rbins(1:end-1), rmedB, 'k--', 'LineWidth',1.5);
plot(ax, rbins(1:end-1), rmedA, 'b-',  'LineWidth',2);

fillBand(ax, rbins, rloB, rhiB, [0.4 0.4 0.4], 0.15);
fillBand(ax, rbins, rloA, rhiA, [0.1 0.5 0.9], 0.20);

yline(ax, thr, '--', 'threshold');

xlabel(ax, 'Distance from array centre (m)');
ylabel(ax, ylab);
title(ax, 'Accuracy vs distance from array centre');

legend(ax, {'before (points)','after (points)', ...
    'before median','after median','before 10–90%','after 10–90%'}, ...
    'Location','northwest');

hold(ax,'off');
end

function plotMicDriftAxes(axX, axY, axZ, itVals, dx, dy, dz, accepted, isShake)

plotOne(axX, itVals, dx, accepted, isShake, '\Delta x (m)', 'Mic drift: X');
plotOne(axY, itVals, dy, accepted, isShake, '\Delta y (m)', 'Mic drift: Y');
plotOne(axZ, itVals, dz, accepted, isShake, '\Delta z (m)', 'Mic drift: Z');

    function plotOne(ax, itVals, D, accepted, isShake, ylab, ttl)
        hold(ax,'on'); grid(ax,'on');

        nM = size(D,2);

        % --- mic lines ---
        hLines = gobjects(nM,1);
        for m = 1:nM
            hLines(m) = plot(ax, itVals, D(:,m), 'LineWidth', 1.2);
        end

        % --- event markers (only if present) ---
        hAcc = gobjects(1); hShake = gobjects(1);
        hasAcc = any(accepted);
        hasShake = any(isShake);

        if hasAcc
            hAcc = scatter(ax, itVals(accepted), zeros(sum(accepted),1), 18, 'filled');
        end
        if hasShake
            hShake = scatter(ax, itVals(isShake), zeros(sum(isShake),1), 30, 'd');
        end

        xlabel(ax,'Iteration (it)');
        ylabel(ax, ylab);
        title(ax, ttl);

        % --- legend: only include what exists ---
        h = hLines;
        labs = arrayfun(@(m) sprintf('M%d', m), 1:nM, 'UniformOutput', false);

        if hasAcc
            h(end+1) = hAcc; %#ok<AGROW>
            labs{end+1} = 'accepted';
        end
        if hasShake
            h(end+1) = hShake; %#ok<AGROW>
            labs{end+1} = 'shake/reset';
        end

        lgd = legend(ax, h, labs, 'Location','best', 'Interpreter','none');
        lgd.Color = 'none';        % transparent background
        lgd.Box   = 'off';         % remove border
        lgd.TextColor = [0 0 0];   % ensure readability
        hold(ax,'off');
    end
end

function thr_ang = deriveAngularThreshold(angVals)
% Robust heuristic: define a "threshold-like" scale for angular plots.
angVals = angVals(isfinite(angVals));
if isempty(angVals)
    thr_ang = 5; % fallback
    return;
end
m  = median(angVals);
p95 = prctile(angVals,95);
thr_ang = max(1e-3, m + 0.5*(p95 - m));
end

function makeMicStartEndTableFigure(M0, M1, AFTER)
% Table with start/end xyz + drifts for fast adoption
if isempty(M0) || isempty(M1) || size(M0,1) ~= size(M1,1)
    return;
end
N = size(M0,1);

d = M1 - M0;
dR = vecnorm(d,2,2);

T = table((1:N)', ...
    M0(:,1), M0(:,2), M0(:,3), ...
    M1(:,1), M1(:,2), M1(:,3), ...
    d(:,1),  d(:,2),  d(:,3),  dR, ...
    'VariableNames', {'Mic','x_start','y_start','z_start','x_end','y_end','z_end','dx','dy','dz','dR'});

fT = figure('Name','Mic start/end table','Color',[1 1 1], 'Position',[120 120 980 320]);
uit = uitable(fT, 'Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, ...
    'Units','normalized', 'Position',[0.01 0.01 0.98 0.98]);

titleStr = sprintf('Microphone coordinates: start vs end (%s)', string(AFTER));
annotation(fT,'textbox',[0.01 0.92 0.98 0.08], ...
    'String', titleStr, 'EdgeColor','none', 'FontWeight','bold', 'FontSize',12);
end

%% ===================== Core helper functions  =====================

function [vals, labelStr] = pickMetric(stats, METRIC)
if string(METRIC) == "ang"
    vals = stats.ang_deg;
    labelStr = 'Angular error (deg)';
else
    vals = stats.err_cm;
    labelStr = 'Position error (cm)';
end
vals = vals(:);
end

function idx = maxPassTieBreak(pass, p95, J)
maxP = max(pass);
cand = find(pass == maxP);
if numel(cand) == 1, idx = cand; return; end
[~,ord] = sort(p95(cand),'ascend');
cand = cand(ord);
bestP95 = p95(cand(1));
cand2 = cand(p95(cand) == bestP95);
if numel(cand2) == 1, idx = cand2; return; end
[~,ord2] = sort(J(cand2),'descend');
idx = cand2(ord2(1));
end

function plot3D(ax, P, M, gridPts, vals, usedMask, labelStr, thr, cmax)
cla(ax); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax,3);
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');

[sx,sy,sz] = sphere(30);
surf(ax, P.Lmax*sx, P.Lmax*sy, P.Lmax*sz, 'FaceAlpha',0.05,'EdgeAlpha',0.12);
surf(ax, P.R*sx,    P.R*sy,    P.R*sz,    'FaceAlpha',0.03,'EdgeAlpha',0.08);
plot3(ax,0,0,0,'pentagram','MarkerSize',6,'LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','none');

cmap = lines(size(M,1));
for ii = 1:size(M,1)
    scatter3(ax, M(ii,1),M(ii,2),M(ii,3), 40, cmap(ii,:), 'filled','MarkerEdgeColor','k','LineWidth',1.5);
    text(ax, M(ii,1),M(ii,2),M(ii,3), sprintf('  %d',ii), 'FontWeight','bold','FontSize',10,'Color',cmap(ii,:));
end

pts = gridPts(usedMask,:);
vv  = vals(usedMask);
vvCap = min(vv, cmax);
scatter3(ax, pts(:,1),pts(:,2),pts(:,3), 28, vvCap, 'filled', 'MarkerFaceAlpha',0.75);

colormap(ax, turbo);
cb = colorbar(ax);
cb.Label.String = labelStr;
cb.Label.FontWeight = 'bold';
clim(ax,[0, cmax]);

text(ax, 0.02, 0.98, sprintf('thr=%.2f | 2×thr=%.2f', thr, cmax), ...
    'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontWeight','bold','BackgroundColor',[1 1 1 0.6],'Margin',4);

hold(ax,'off');
end

function [med, qlo, qhi] = radialQuantiles(r, e, rbins)
n = numel(rbins) - 1;
med = nan(n,1); qlo = nan(n,1); qhi = nan(n,1);
for i = 1:n
    m = r>=rbins(i) & r<rbins(i+1) & isfinite(e);
    if any(m)
        med(i) = median(e(m));
        qlo(i) = prctile(e(m),10);
        qhi(i) = prctile(e(m),90);
    end
end
end

function fillBand(ax, x, lo, hi, col, alpha)
x  = x(:); lo = lo(:); hi = hi(:);
if numel(x) == numel(lo)+1
    x = (x(1:end-1) + x(2:end))/2;
end
ok = isfinite(x) & isfinite(lo) & isfinite(hi);
x = x(ok); lo = lo(ok); hi = hi(ok);
if numel(x) < 2, return; end
fill(ax, [x; flipud(x)], [lo; flipud(hi)], col, 'FaceAlpha', alpha, 'EdgeColor','none');
end

function [M0, dM, dx, dy, dz] = micDriftXYZ(H)
% Returns baseline mic positions, per-iteration mic deltas, and dx/dy/dz matrices.
nI = numel(H.iter);
M0 = H.iter(1).M;
nM = size(M0,1);
dM = zeros(nI, nM, 3);

for i = 1:nI
    Mi = H.iter(i).M;
    if isempty(Mi) || size(Mi,1) ~= nM
        continue;
    end
    Di = Mi - M0;
    dM(i,:,:) = reshape(Di, [1 nM 3]);
end

dx = squeeze(dM(:,:,1));
dy = squeeze(dM(:,:,2));
dz = squeeze(dM(:,:,3));
end

function cb = mkCB(parent,r,c,txt,val)
cb = uicheckbox(parent,'Text',txt,'Value',val);
cb.Layout.Row = r; cb.Layout.Column = c;
end

function dirs = ruleSet2Dirs()
dirs = [
    1 0 0; -1 0 0;
    0 1 0;  0 -1 0;
    0 0 1;  0 0 -1;
    1 1 0;  1 -1 0; -1 1 0; -1 -1 0;
    1 0 1;  1 0 -1;
    0 1 1;  0 1 -1
    ];
dirs = dirs ./ vecnorm(dirs,2,2);
end


function M = randomArrayNonCoplanar(N, Lmax, minPairDist)
for t = 1:5000
    M = randn(N,3);
    M = M ./ vecnorm(M,2,2);
    radii = (0.25 + 0.75*rand(N,1)) * Lmax;
    M = M .* radii;
    M = M - mean(M,1);
    if isValidGeometry(M, Lmax, minPairDist)
        return;
    end
end
error('Could not generate a valid non-coplanar geometry (try loosening constraints).');
end

function ok = isValidGeometry(M, Lmax, minPairDist)
ok = true;
Mc = M - mean(M,1);

% within arm limit
r = vecnorm(Mc,2,2);
if any(r > Lmax*(1+1e-9))
    ok = false;
    return;
end

% min pairwise spacing
if minPairDist > 0
    D = squareform(pdist(Mc));
    D(D==0) = inf;
    if any(D(:) < minPairDist)
        ok = false;
        return;
    end
end

% non-coplanar
if rank(Mc,1e-6) < 3
    ok = false;
    return;
end
end

function gridPts = buildShellGrid(R, Rin, step, quadrants)
% Generate 3D grid within spherical shell, then apply quadrant filtering.
ax = -R:step:R;
[X,Y,Z] = ndgrid(ax,ax,ax);
pts = [X(:) Y(:) Z(:)];
r = vecnorm(pts,2,2);

% Shell constraint
shellMask = (r <= R) & (r >= Rin);
pts = pts(shellMask,:);

if nargin < 4 || isempty(quadrants)
    gridPts = pts;
    return;
end

% Quadrants:
% Front: X > 0, Back: X < 0
% Right: Y > 0, Left: Y < 0
% Top: Z > 0, Bottom: Z < 0
quadMask = false(size(pts,1),1);

if quadrants.frontTopLeft
    quadMask = quadMask | (pts(:,1) > 0 & pts(:,2) < 0 & pts(:,3) > 0);
end
if quadrants.frontTopRight
    quadMask = quadMask | (pts(:,1) > 0 & pts(:,2) > 0 & pts(:,3) > 0);
end
if quadrants.frontBottomLeft
    quadMask = quadMask | (pts(:,1) > 0 & pts(:,2) < 0 & pts(:,3) < 0);
end
if quadrants.frontBottomRight
    quadMask = quadMask | (pts(:,1) > 0 & pts(:,2) > 0 & pts(:,3) < 0);
end
if quadrants.backTopLeft
    quadMask = quadMask | (pts(:,1) < 0 & pts(:,2) < 0 & pts(:,3) > 0);
end
if quadrants.backTopRight
    quadMask = quadMask | (pts(:,1) < 0 & pts(:,2) > 0 & pts(:,3) > 0);
end
if quadrants.backBottomLeft
    quadMask = quadMask | (pts(:,1) < 0 & pts(:,2) < 0 & pts(:,3) < 0);
end
if quadrants.backBottomRight
    quadMask = quadMask | (pts(:,1) < 0 & pts(:,2) > 0 & pts(:,3) < 0);
end

gridPts = pts(quadMask,:);
end

function Mout = shakeGeometry(M, sigma, Lmax, minPairDist)
Mout = M + sigma*randn(size(M));
Mout = Mout - mean(Mout,1);

r = vecnorm(Mout,2,2);
mask = r > Lmax;
Mout(mask,:) = Mout(mask,:) .* (Lmax ./ r(mask)) * 0.999;

if ~isValidGeometry(Mout, Lmax, minPairDist)
    Mout = M;
end
end


function J = scoreJ(stats, lambdaMean, failPenalty_cm)
meanNorm = min(stats.meanErr / failPenalty_cm, 1);
J = stats.passRate - lambdaMean * meanNorm;
end

function tdoa = estimateTDOA_xcorr(signals, fs)
% signals: [Nsamp x Nmics]
num_mics = size(signals,2);
tdoa = zeros(num_mics-1,1);
ref = signals(:,1);
for i = 2:num_mics
    [c,lags] = xcorr(signals(:,i), ref, 'coeff');
    [~,idx] = max(abs(c));
    tdoa(i-1) = lags(idx)/fs;
end
end

