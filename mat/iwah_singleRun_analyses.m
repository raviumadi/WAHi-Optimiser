%% iwah_analyse_single_run.m
% Adds: distance-from-centre analysis as an extra subplot:
%   - overlays BEFORE vs AFTER histograms *within distance bins*
%   - uses transparency

clear; clc;

%% ===================== USER OPTIONS =====================
RUN_FILE = "../data/single_runs/iwah_gui_run__4mics__20260201_151117.mat";  % pyramid random rng

AFTER    = "final";   % "final" | "best"
METRIC   = "pos";     % "pos" (cm) | "ang" (deg)

CAT_MODE = "thresholdPlus4Tol";  % "thresholdPlus4Tol" | "4xThreshold"
USE_GRID_FROM = "baseline";      % "baseline" | "after"
IGNORE_UNUSED = true;

% Distance-binned distribution settings
DIST = struct();
DIST.nBins = 4;                 % number of radial shells
DIST.binEdges = [];             % leave empty to auto-compute between Rin..R
DIST.useRinR = true;            % if true, uses [P.Rin, P.R] for auto edges; else uses min/max of data
DIST.histBins = 35;             % number of histogram bins for metric
DIST.alpha = 0.35;              % transparency for overlays
DIST.showLegends = true;

%% ===================== LOAD =====================
S = load(RUN_FILE);
if ~isfield(S,'H'); error('File does not contain variable H'); end
H = S.H;

P = H.meta.P;
if ~isfield(P,'threshold_cm')
    P.threshold_cm = P.best_cm + P.tol_cm;
end

iBase  = 1;
iFinal = numel(H.iter);

passAll = arrayfun(@(x) x.stats.passRate, H.iter(:))';
p95All  = arrayfun(@(x) x.stats.p95Err,    H.iter(:))';
JAll    = arrayfun(@(x) x.J,              H.iter(:))';
iBest   = maxPassTieBreak(passAll, p95All, JAll);

switch AFTER
    case "best",  iAfter = iBest;
    otherwise,    iAfter = iFinal;
end

base = H.iter(iBase);
aft  = H.iter(iAfter);

switch USE_GRID_FROM
    case "after",    gridPts = aft.gridPts;
    otherwise,       gridPts = base.gridPts;
end

%% ===================== THRESHOLDS =====================
thr  = P.threshold_cm;
cmax = 2*thr;

switch CAT_MODE
    case "4xThreshold"
        catThr = 4*thr;
    otherwise
        catThr = thr + 4*P.tol_cm;
end

[valsBase, labelStr] = pickMetric(base.stats, METRIC);
[valsAft,  ~]        = pickMetric(aft.stats,  METRIC);

maskBase = true(size(valsBase));
maskAft  = true(size(valsAft));
if IGNORE_UNUSED
    if isfield(base.stats,'usedMask'), maskBase = base.stats.usedMask(:); end
    if isfield(aft.stats,'usedMask'),  maskAft  = aft.stats.usedMask(:);  end
end

%% ===================== PASS TRAJECTORY =====================
itVals = arrayfun(@(x) x.it, H.iter(:))';
accepted = false(numel(H.iter),1);
note = strings(numel(H.iter),1);
for i = 1:numel(H.iter)
    if isfield(H.iter(i),'accepted'), accepted(i) = logical(H.iter(i).accepted); end
    if isfield(H.iter(i),'note'),     note(i) = string(H.iter(i).note); end
end
isShake = contains(note,"shake/reset");

%% ===================== DISTANCE BINNING PREP =====================
% distance from array centre (origin, after recentering)
rAll = vecnorm(gridPts,2,2);

if isempty(DIST.binEdges)
    if DIST.useRinR && isfield(P,'Rin') && isfield(P,'R')
        r0 = P.Rin; r1 = P.R;
    else
        % use actual range of used points only
        rUsed = rAll(maskBase | maskAft);
        r0 = min(rUsed); r1 = max(rUsed);
    end
    DIST.binEdges = linspace(r0, r1, DIST.nBins+1);
end

% metric bins for histograms (shared across panels for comparability)
vAll = [valsBase(maskBase); valsAft(maskAft)];
vAll = vAll(isfinite(vAll));
if isempty(vAll)
    error('No finite metric values found for histograms (check masks).');
end
vMaxForH = max(cmax, prctile(vAll, 99));              % avoid one insane outlier ruining the view
metricEdges = linspace(0, vMaxForH, DIST.histBins+1); % fixed edges => comparable overlays

%% ===================== FIGURE =====================
fig = figure('Name','iWAH Single-run analysis','Color',[1 1 1], ...
    'Position',[60 60 1650 980]);

t = tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');

% (1) 3D BEFORE
ax1 = nexttile(t,1);
plot3D(ax1, P, base.M, gridPts, valsBase, maskBase, labelStr, thr, cmax, catThr);
title(ax1, sprintf('BEFORE (it=%.1f): pass=%.1f%%, mean=%.1f, p95=%.1f', ...
    base.it, 100*base.stats.passRate, base.stats.meanErr, base.stats.p95Err));
xlim([-0.5 2]); ylim([-2 2]); zlim([-0.5 2]); 

% (2) 3D AFTER
ax2 = nexttile(t,2);
plot3D(ax2, P, aft.M, gridPts, valsAft, maskAft, labelStr, thr, cmax, catThr);
title(ax2, sprintf('AFTER (%s, it=%.1f): pass=%.1f%%, mean=%.1f, p95=%.1f', ...
    AFTER, aft.it, 100*aft.stats.passRate, aft.stats.meanErr, aft.stats.p95Err));
xlim([-0.5 2]); ylim([-2 2]); zlim([-0.5 2]); 

% (3) Pass trajectory
ax3 = nexttile(t,3);
grid(ax3,'on'); hold(ax3,'on');
plot(ax3, itVals, 100*passAll, 'LineWidth', 1.5);
scatter(ax3, itVals(accepted), 100*passAll(accepted), 45, 'filled');
scatter(ax3, itVals(isShake),  100*passAll(isShake),  55, 'd');
% yline(ax3, 100*P.targetPass, '--', 'Target');
xlabel(ax3,'Iteration (it)'); ylabel(ax3,'Pass rate (%)');
title(ax3,'Pass-rate progression');
legend(ax3, {'pass','accepted','shake/reset','target'}, 'Location','best');
hold(ax3,'off');

% (4) Global distribution before vs after (VALID region only)
ax4 = nexttile(t,4);
grid(ax4,'on'); hold(ax4,'on');

thr  = P.threshold_cm;
thr2 = 4*thr;

% --- separate valid vs catastrophic ---
vB_all = valsBase(maskBase);
vA_all = valsAft(maskAft);

vB_valid = vB_all(vB_all <= thr2);
vA_valid = vA_all(vA_all <= thr2);

nB_cat = sum(vB_all > thr2);
nA_cat = sum(vA_all > thr2);

fracB_cat = nB_cat / max(1,numel(vB_all));
fracA_cat = nA_cat / max(1,numel(vA_all));

% --- histogram (valid region only) ---
edges = linspace(0, thr2, 40);

histogram(ax4, vB_valid, edges, ...
    'Normalization','probability', ...
    'FaceAlpha',0.35,'EdgeAlpha',0.15);

histogram(ax4, vA_valid, edges, ...
    'Normalization','probability', ...
    'FaceAlpha',0.35,'EdgeAlpha',0.15);

% --- reference lines ---
xline(ax4, thr,  '--k','threshold');
xline(ax4, thr2, ':k','4×threshold');

xlabel(ax4, labelStr);
ylabel(ax4,'Probability');

title(ax4, sprintf('Global accuracy (≤ 4×thr = %.1f cm)', thr2));

legend(ax4, {'before','after','threshold','4×threshold'}, ...
    'Location','best');

% --- catastrophic failure annotation ---
txt = sprintf(['Catastrophic failures (>4×thr)\n' ...
               'before: %.1f%% (%d/%d)\n' ...
               'after:  %.1f%% (%d/%d)'], ...
               100*fracB_cat, nB_cat, numel(vB_all), ...
               100*fracA_cat, nA_cat, numel(vA_all));

text(ax4, 0.98, 0.98, txt, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',9, ...
    'BackgroundColor',[1 1 1 0.85], ...
    'EdgeColor',[0.7 0.7 0.7]);

hold(ax4,'off');
% (5) Mic geometry scatter before vs after
ax5 = nexttile(t,5);
grid(ax5,'on'); axis(ax5,'equal'); hold(ax5,'on'); view(ax5,3);
xlabel(ax5,'X (m)'); ylabel(ax5,'Y (m)'); zlabel(ax5,'Z (m)');
title(ax5,'Mic geometry: before (dull) vs after (bright)');
beforeCol = repmat([0.55 0.55 0.55], size(base.M,1), 1);
afterCol  = lines(size(aft.M,1));
scatter3(ax5, base.M(:,1), base.M(:,2), base.M(:,3), 120, beforeCol, 'filled', ...
    'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth',1.2);
scatter3(ax5, aft.M(:,1),  aft.M(:,2),  aft.M(:,3),  120, afterCol, 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.4);
legend(ax5, {'before','after'}, 'Location','best');
hold(ax5,'off');

% (6) NEW: Accuracy distributions as a function of distance (r)
% === Subplot: Error vs distance from array centre ===
ax6 = nexttile(t,6);
hold(ax6,'on'); grid(ax6,'on');

rBase = vecnorm(gridPts(maskBase,:),2,2);
rAft  = vecnorm(gridPts(maskAft,:),2,2);

eBase = valsBase(maskBase);
eAft  = valsAft(maskAft);

% Cap for visibility (do NOT hide catastrophes, just clip display)
eCap = 2 * P.threshold_cm;
eBaseP = min(eBase, eCap);
eAftP  = min(eAft,  eCap);

scatter(ax6, rBase, eBaseP, 14, [0.6 0.6 0.6], 'filled', ...
    'MarkerFaceAlpha',0.25);
scatter(ax6, rAft,  eAftP,  14, [0.1 0.5 0.9], 'filled', ...
    'MarkerFaceAlpha',0.35);

% Distance bins for summary curves
rbins = linspace(P.Rin, P.R, 10);
[rmedB, rloB, rhiB] = radialQuantiles(rBase, eBase, rbins);
[rmedA, rloA, rhiA] = radialQuantiles(rAft,  eAft,  rbins);

plot(ax6, rbins(1:end-1), rmedB, 'k--', 'LineWidth',1.5);
plot(ax6, rbins(1:end-1), rmedA, 'b-',  'LineWidth',2);

fillBand(ax6, rbins, rloB, rhiB, [0.4 0.4 0.4], 0.15);
fillBand(ax6, rbins, rloA, rhiA, [0.1 0.5 0.9], 0.20);

yline(ax6, P.threshold_cm, '--', 'threshold');

xlabel(ax6,'Distance from array centre (m)');
ylabel(ax6,'Position error (cm)');
title(ax6,'Accuracy vs distance from array centre');

legend(ax6, ...
    {'before (points)','after (points)', ...
     'before median','after median','before 10–90%','after 10–90%'}, ...
    'Location','northwest');

%% ===================== LOCAL HELPERS =====================

function plot3D(ax, P, M, gridPts, vals, usedMask, labelStr, thr, cmax, catThr)
cla(ax); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax,3);
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');

[sx,sy,sz] = sphere(30);
surf(ax, P.Lmax*sx, P.Lmax*sy, P.Lmax*sz, 'FaceAlpha',0.05,'EdgeAlpha',0.12);
surf(ax, P.R*sx,    P.R*sy,    P.R*sz,    'FaceAlpha',0.03,'EdgeAlpha',0.08);

plot3(ax,0,0,0,'pentagram','MarkerSize',6,'LineWidth',2,...
    'MarkerFaceColor','r','MarkerEdgeColor','none');

cmap = lines(size(M,1));
for ii = 1:size(M,1)
    scatter3(ax, M(ii,1),M(ii,2),M(ii,3), 40, cmap(ii,:), 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.5);
    text(ax, M(ii,1),M(ii,2),M(ii,3), sprintf('  %d',ii), ...
        'FontWeight','bold','FontSize',10,'Color',cmap(ii,:));
end

pts = gridPts(usedMask,:);
vv  = vals(usedMask);

exceedMask = vv > cmax;
catMask    = vv > catThr;

vvCap = min(vv, cmax);
scatter3(ax, pts(:,1),pts(:,2),pts(:,3), 28, vvCap, 'filled', 'MarkerFaceAlpha',0.75);

if any(exceedMask)
    scatter3(ax, pts(exceedMask,1),pts(exceedMask,2),pts(exceedMask,3), 20, 'm', 'filled', 'MarkerEdgeColor', 'none');
end
if any(catMask)
    scatter3(ax, pts(catMask,1),pts(catMask,2),pts(catMask,3), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.25);
end

colormap(ax, turbo);
cb = colorbar(ax);
cb.Label.String = labelStr;
cb.Label.FontWeight = 'bold';
caxis(ax,[0, cmax]);

text(ax, 0.02, 0.98, sprintf('thr=%.1f | 2×thr=%.1f | cat=%.1f', thr, cmax, catThr), ...
    'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'FontWeight','bold', 'BackgroundColor',[1 1 1 0.6], 'Margin',4);

hold(ax,'off');
end

function [vals, labelStr] = pickMetric(stats, METRIC)
if METRIC == "ang"
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
x  = x(:);  lo = lo(:);  hi = hi(:);

% If x is bin edges and lo/hi are per-bin, convert to centres
if numel(x) == numel(lo)+1
    x = (x(1:end-1) + x(2:end))/2;
end

if ~(numel(x)==numel(lo) && numel(lo)==numel(hi))
    error('fillBand:SizeMismatch', ...
        'x (%d), lo (%d), hi (%d) must match (or x can be edges = lo+1).', ...
        numel(x), numel(lo), numel(hi));
end

ok = isfinite(x) & isfinite(lo) & isfinite(hi);
x = x(ok); lo = lo(ok); hi = hi(ok);

if numel(x) < 2, return; end

fill(ax, [x; flipud(x)], [lo; flipud(hi)], col, ...
    'FaceAlpha', alpha, 'EdgeColor','none');
end