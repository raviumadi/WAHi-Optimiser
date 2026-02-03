function [H, runInfo] = iwah_run_optimiser_core(P, OPT, Q, CSV, seed)
% iwah_run_optimiser_core
% Core optimiser runner for batch experiments (NO plotting).
%
% Inputs:
%   P, OPT, Q : your parameter structs
%   CSV       : struct with .useCSVStart, .path, repair options, etc.
%   seed      : rng seed used for this run (stored for reproducibility)
%
% Outputs:
%   H         : full history struct (same as your standalone)
%   runInfo   : lightweight metadata

%% ===================== INITIALISE RNG =====================
if nargin < 5 || isempty(seed)
    seed = randi(2^31-1);
end
rng(seed,'twister');

% Fix f0/f1 ordering
if P.f0 > P.f1
    tmp = P.f0; P.f0 = P.f1; P.f1 = tmp;
end
if P.Rin >= P.R, error('P.Rin must be < P.R'); end
if P.step <= 0, error('P.step must be > 0'); end
if ~any(structfun(@(x)x,Q)), error('No quadrants enabled.'); end

% Build grid
gridAll = buildShellGrid(P.R, P.Rin, P.step, Q);
if isempty(gridAll)
    error('Empty grid. Check R/Rin/step/quadrants.');
end

% Base directions (Rule-Set-2)
baseDirs = [
    1 0 0; -1 0 0;
    0 1 0;  0 -1 0;
    0 0 1;  0 0 -1;
    1 1 0;  1 -1 0; -1 1 0; -1 -1 0;
    1 0 1;  1 0 -1;
    0 1 1;  0 1 -1
];
baseDirs = baseDirs ./ vecnorm(baseDirs,2,2);

%% ===================== INITIAL GEOMETRY (CSV OR RANDOM) =====================
if isfield(CSV,'useCSVStart') && CSV.useCSVStart
    [M, startNote] = loadMicGeometryCSV(CSV.path, CSV.hasHeader, CSV.recenter);

    if CSV.enforceConstraints && ~isValidGeometry(M, P.Lmax, P.minPairDist)
        if isfield(CSV,'autoRepair') && CSV.autoRepair
            [Mfix, ok, rep] = repairMicGeometry(M, P.Lmax, P.minPairDist, CSV);
            if ok
                M = Mfix;
                startNote = sprintf('%s (repaired: %s)', startNote, rep.note);
            else
                if isfield(CSV,'repairFallbackToRandom') && CSV.repairFallbackToRandom
                    M = randomArrayNonCoplanar(size(M,1), P.Lmax, P.minPairDist);
                    M = M - mean(M,1);
                    startNote = sprintf('%s (repair failed -> random fallback)', startNote);
                else
                    error('CSV geometry violates constraints and repair failed (%s).', rep.note);
                end
            end
        else
            error('CSV geometry violates constraints (Lmax/minPairDist/coplandarity).');
        end
    end

    P.N = size(M,1); % adopt N from CSV
else
    M = randomArrayNonCoplanar(P.N, P.Lmax, P.minPairDist);
    M = M - mean(M,1);
    startNote = "random start";
end

%% ===================== HISTORY STORAGE =====================
H = struct();
H.meta = struct();
H.meta.rngSeed   = seed;
H.meta.timeStart = datestr(now);
H.meta.P         = P;
H.meta.OPT       = OPT;
H.meta.Q         = Q;

H.meta.start = struct();
H.meta.start.source = char(startNote);
H.meta.start.csvPath = '';
if isfield(CSV,'useCSVStart') && CSV.useCSVStart, H.meta.start.csvPath = CSV.path; end
H.meta.start.M = M;

H.iter = struct( ...
    'it',{}, 'M',{}, 'micIdx',{}, 'stepSize',{}, 'dirsTried',{}, ...
    'P',{}, 'OPT',{}, 'Q',{}, 'gridPts',{}, ...
    'stats',{}, 'diag',{}, 'J',{}, ...
    'bestLocal',{}, 'accepted',{}, 'note',{} );

%% ===================== LOCALISER PARAM TEMPLATE =====================
paramTemplate = struct();
paramTemplate.fs = P.fs;
paramTemplate.d  = P.d;
paramTemplate.f0 = P.f0;
paramTemplate.f1 = P.f1;
paramTemplate.tail = P.tail;
paramTemplate.snr_db = P.snr_db;
paramTemplate.callType = P.callType;
paramTemplate.velocity = P.velocity;

%% ===================== BASELINE EVAL =====================
paramBase = paramTemplate;
paramBase.mic_positions = M;
locBase = BatCallLocaliser(paramBase);

[bestStats, bestDiag] = evaluateFoA_full(locBase, gridAll, P.threshold_cm, ...
    P.minDistToAnyMic, P.epsDist, P.failPenalty_cm, P.c);

bestJ = scoreJ(bestStats, P.lambdaMean, P.failPenalty_cm);

H.iter(end+1) = makeIterRecord(0, M, NaN, NaN, [], P, OPT, Q, gridAll, bestStats, bestDiag, bestJ, [], true, char(startNote));

%% ===================== OPTIMISER LOOP =====================
stall = 0;

for it = 1:OPT.maxIters
    if OPT.maxIters > 1
        stepFrac = max(OPT.minStepFrac, OPT.initStepFrac * (1 - (it-1)/(OPT.maxIters-1)));
    else
        stepFrac = OPT.initStepFrac;
    end
    stepSize = stepFrac * P.Lmax;

    micIdx = mod(it-1, size(M,1)) + 1;

    dirs = baseDirs(1:min(OPT.probesPerMic, size(baseDirs,1)),:);

    v = M(micIdx,:);
    nv = norm(v);
    radial = [1 0 0];
    if nv > 1e-9, radial = v / nv; end
    dirs = [dirs; radial; -radial];

    bestLocalJ = -inf;
    bestLocalM = [];
    bestLocalStats = [];
    bestLocalDiag = [];

    for k = 1:size(dirs,1)
        Mtry = M;
        Mtry(micIdx,:) = Mtry(micIdx,:) + stepSize * dirs(k,:);
        Mtry = Mtry - mean(Mtry,1);

        r = vecnorm(Mtry,2,2);
        if r(micIdx) > P.Lmax
            Mtry(micIdx,:) = Mtry(micIdx,:) * (P.Lmax/r(micIdx)) * 0.999;
            Mtry = Mtry - mean(Mtry,1);
        end

        if ~isValidGeometry(Mtry, P.Lmax, P.minPairDist)
            continue;
        end

        paramTry = paramTemplate;
        paramTry.mic_positions = Mtry;
        locTry = BatCallLocaliser(paramTry);

        [stTry, dgTry] = evaluateFoA_full(locTry, gridAll, P.threshold_cm, ...
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
    if ~isempty(bestLocalM) && (bestLocalJ > bestJ + OPT.improveEpsJ)
        M = bestLocalM;
        bestJ = bestLocalJ;
        bestStats = bestLocalStats;
        bestDiag  = bestLocalDiag;
        stall = 0;
        accepted = true;
        note = "accepted (improved J)";
    else
        stall = stall + 1;
        note = sprintf('rejected (stall=%d)', stall);
    end

    bestLocal = struct();
    bestLocal.J = bestLocalJ;
    bestLocal.M = bestLocalM;
    bestLocal.stats = bestLocalStats;
    bestLocal.diag  = bestLocalDiag;

    H.iter(end+1) = makeIterRecord(it, M, micIdx, stepSize, dirs, P, OPT, Q, gridAll, bestStats, bestDiag, bestJ, bestLocal, accepted, note);

    if bestStats.passRate >= P.targetPass
        break;
    end

    if stall >= OPT.stallLimit
        Mshake = shakeGeometry(M, OPT.shakeFrac*P.Lmax, P.Lmax, P.minPairDist);

        if ~isValidGeometry(Mshake, P.Lmax, P.minPairDist)
            Mshake = randomArrayNonCoplanar(size(M,1), P.Lmax, P.minPairDist);
            Mshake = Mshake - mean(Mshake,1);
        end

        paramShake = paramTemplate;
        paramShake.mic_positions = Mshake;
        locShake = BatCallLocaliser(paramShake);

        [stShake, dgShake] = evaluateFoA_full(locShake, gridAll, P.threshold_cm, ...
            P.minDistToAnyMic, P.epsDist, P.failPenalty_cm, P.c);
        Jshake = scoreJ(stShake, P.lambdaMean, P.failPenalty_cm);

        M = Mshake;
        bestStats = stShake;
        bestDiag  = dgShake;
        bestJ     = Jshake;
        stall = 0;

        H.iter(end+1) = makeIterRecord(it+0.5, M, NaN, NaN, [], P, OPT, Q, gridAll, bestStats, bestDiag, bestJ, [], true, 'shake/reset');
    end
end

H.meta.timeEnd = datestr(now);
H.meta.final = struct('M',M,'stats',bestStats,'diag',bestDiag,'J',bestJ);

runInfo = struct();
runInfo.seed = seed;
runInfo.N = size(M,1);
runInfo.startNote = char(startNote);
if isfield(CSV,'useCSVStart') && CSV.useCSVStart
    runInfo.csvPath = CSV.path;
else
    runInfo.csvPath = '';
end
end

%% ===================== LOCAL FUNCTIONS =====================

function rec = makeIterRecord(it, M, micIdx, stepSize, dirs, P, OPT, Q, gridPts, stats, diag, J, bestLocal, accepted, note)
rec = struct();
rec.it = it;
rec.M = M;                     % [N x 3] mic positions (m)
rec.micIdx = micIdx;           % mic index moved (or NaN)
rec.stepSize = stepSize;       % step size used (m)
rec.dirsTried = dirs;          % directions attempted (rows)
rec.P = P;                     % full inputs (copy)
rec.OPT = OPT;                 % optimiser settings (copy)
rec.Q = Q;                     % quadrant selection (copy)
rec.gridPts = gridPts;         % evaluated grid points (m)
rec.stats = stats;             % stats struct (includes full err vectors + passMask)
rec.diag = diag;               % diag struct (fail counts etc.)
rec.J = J;                     % objective
rec.bestLocal = bestLocal;     % best local candidate info (may be empty)
rec.accepted = accepted;       % accepted into current state?
rec.note = note;               % string note
end

function updatePlot(ax, P, PLOT, M, gridPts, stats)

% Ensure correct figure is current (important if other figures open)
fig = ancestor(ax,'figure');
if ~isempty(fig) && isvalid(fig)
    figure(fig);
end

cla(ax,'reset');  % reset is sometimes more reliable than cla(ax)
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax,3);

% shells
[sx,sy,sz] = sphere(30);
surf(ax, P.Lmax*sx, P.Lmax*sy, P.Lmax*sz, 'FaceAlpha',0.05,'EdgeAlpha',0.12);
surf(ax, P.R*sx,    P.R*sy,    P.R*sz,    'FaceAlpha',0.03,'EdgeAlpha',0.08);

% mic positions
cmap = lines(size(M,1));
for ii = 1:size(M,1)
    scatter3(ax, M(ii,1),M(ii,2),M(ii,3), 120, cmap(ii,:), 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.5);
    text(ax, M(ii,1),M(ii,2),M(ii,3), sprintf('  M%d',ii), ...
        'FontWeight','bold','FontSize',10,'Color',cmap(ii,:));
end

% origin
plot3(ax,0,0,0,'pentagram','MarkerSize',12,'LineWidth',2,...
    'MarkerFaceColor','r','MarkerEdgeColor','k');

% metric values
if PLOT.metric == "pos"
    vals = stats.err_cm;
    ttl = 'Position Error (cm)';
else
    vals = stats.ang_deg;
    ttl = 'Angular Error (deg)';
end

scatter3(ax, gridPts(:,1),gridPts(:,2),gridPts(:,3), 28, vals, ...
    'filled', 'MarkerFaceAlpha', PLOT.markerAlpha);

colormap(ax, turbo);

% Create one colourbar per axes only, reuse it
cb = findall(fig,'Type','ColorBar','-and','Axes',ax);
if isempty(cb)
    cb = colorbar(ax);
end
cb.Label.String = ttl;
cb.FontWeight = 'bold';

if isfield(PLOT,'lockCaxisToThreshold') && PLOT.lockCaxisToThreshold && PLOT.metric == "pos"
    caxis(ax,[0, max(P.threshold_cm, eps)]);
end

xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');
hold(ax,'off');
end

function [M, note] = loadMicGeometryCSV(path, hasHeader, doRecenter)
if ~isfile(path)
    error('CSV file not found: %s', path);
end

if hasHeader
    T = readtable(path);
    if width(T) < 3
        error('CSV table must contain at least 3 columns for x,y,z.');
    end
    M = table2array(T(:,1:3));
else
    M = readmatrix(path);
end

if size(M,2) ~= 3
    error('CSV must have exactly 3 columns (x,y,z).');
end
if size(M,1) < 4
    error('CSV must have at least 4 microphones (rows).');
end
if any(~isfinite(M(:)))
    error('CSV contains NaN/Inf.');
end

if doRecenter
    M = M - mean(M,1);
end

note = sprintf('CSV start: %s', path);
drawnow; pause(0.001);
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
error('Could not generate valid geometry');
end

function ok = isValidGeometry(M, Lmax, minPairDist)
ok = true;
Mc = M - mean(M,1);

r = vecnorm(Mc,2,2);
if any(r > Lmax*(1+1e-9))
    ok = false; return;
end

if minPairDist > 0
    D = squareform(pdist(Mc));
    D(D==0) = inf;
    if any(D(:) < minPairDist)
        ok = false; return;
    end
end

if rank(Mc,1e-6) < 3
    ok = false; return;
end
end

function gridPts = buildShellGrid(R, Rin, step, quadrants)
ax = -R:step:R;
[X,Y,Z] = ndgrid(ax,ax,ax);
pts = [X(:) Y(:) Z(:)];
r = vecnorm(pts,2,2);

shellMask = (r <= R) & (r >= Rin);
pts = pts(shellMask,:);

if nargin < 4 || isempty(quadrants)
    gridPts = pts; return;
end

quadMask = false(size(pts,1),1);

% Front: X > 0, Back: X < 0
% Left:  Y < 0, Right: Y > 0
% Top:   Z > 0, Bottom: Z < 0
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

function [stats, diag] = evaluateFoA_full(loc, gridPts, threshold_cm, minDistToAnyMic, epsDist, failPenalty, c)
nP = size(gridPts,1);

err_cm  = nan(nP,1);
ang_deg = nan(nP,1);

simFail = 0;
tdoaFail = 0;
solveFail = 0;
tooClose = 0;

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
        err_cm(i)  = failPenalty;
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

        e = norm(est - src) * 100; % cm
        if ~isfinite(e)
            solveFail = solveFail + 1;
            continue;
        end
        err_cm(i) = e;

        % angular error (az/el) wrt mic 1
        vT = src - mic_pos(1,:);
        vE = est - mic_pos(1,:);
        if norm(vT) < eps || norm(vE) < eps
            ang_deg(i) = 180;
        else
            azT = atan2d(vT(2),vT(1));
            elT = asind(vT(3)/norm(vT));
            azE = atan2d(vE(2),vE(1));
            elE = asind(vE(3)/norm(vE));

            dAz = mod((azE-azT)+180,360)-180;
            dEl = elE-elT;
            ang_deg(i) = hypot(dAz,dEl);
        end

    catch
        solveFail = solveFail + 1;
        continue;
    end

    % if mod(i,25)==0
    %     drawnow nocallbacks;
    % end
end

% Penalise failures for USED points only
err_f = err_cm;
err_f(isnan(err_f) & usedMask) = failPenalty;

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
stats.usedMask = usedMask;          % IMPORTANT for later analysis

diag = struct();
diag.nP = nP;
diag.gridUsed = sum(usedMask);
diag.tooClose = tooClose;
diag.simFail = simFail;
diag.tdoaFail = tdoaFail;
diag.solveFail = solveFail;
end

function J = scoreJ(stats, lambdaMean, failPenalty)
meanNorm = min(stats.meanErr / failPenalty, 1);
J = stats.passRate - lambdaMean * meanNorm;
end

function tdoa = estimateTDOA_xcorr(signals, fs)
num_mics = size(signals,2);
tdoa = zeros(num_mics-1,1);
ref = signals(:,1);
for i = 2:num_mics
    [c,lags] = xcorr(signals(:,i), ref, 'coeff');
    [~,idx] = max(abs(c));
    tdoa(i-1) = lags(idx)/fs;
end
end

function pos = localiseTDOA_global(tdoa, mic_pos, c)
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

function [Mout, ok, rep] = repairMicGeometry(Min, Lmax, minPairDist, CSV)
% Try to minimally modify a CSV geometry so it satisfies:
%  - all points within Lmax radius (after centring)
%  - min pair distance >= minPairDist
%  - rank >= 3 (non-coplanar)
%
% Strategy:
%  1) recenter
%  2) optional global scaling to satisfy Lmax
%  3) iterative tiny jitter + clamp + recenter until valid

rep = struct();
rep.note = '';
rep.nTries = 0;
rep.finalJitter = NaN;

M = Min;
M = M - mean(M,1);

% (A) Prefer a global scale if any mic exceeds Lmax (this is minimal distortion)
if isfield(CSV,'repairPreferScale') && CSV.repairPreferScale
    r = vecnorm(M,2,2);
    rmax = max(r);
    if rmax > Lmax
        s = (Lmax / rmax) * 0.999;   % tiny margin
        M = M * s;
        M = M - mean(M,1);
    end
end

% Optional hard clamp (also low distortion if only a few points exceed)
if isfield(CSV,'repairClampLmax') && CSV.repairClampLmax
    M = clampToLmax(M, Lmax);
    M = M - mean(M,1);
end

% If valid already, done
if isValidGeometry(M, Lmax, minPairDist)
    Mout = M; ok = true;
    rep.note = 'scale/clamp';
    rep.nTries = 0;
    rep.finalJitter = 0;
    return;
end

% (B) Iterative jitter (minimal random changes)
maxTries = getField(CSV,'repairMaxTries',200);
j0 = getField(CSV,'repairJitter0',1e-3);
grow = getField(CSV,'repairJitterGrowth',1.15);
doClamp = getField(CSV,'repairClampLmax',true);

j = j0;
for t = 1:maxTries
    Mt = M + j*randn(size(M));       % tiny perturbation
    Mt = Mt - mean(Mt,1);

    if doClamp
        Mt = clampToLmax(Mt, Lmax);
        Mt = Mt - mean(Mt,1);
    end

    if isValidGeometry(Mt, Lmax, minPairDist)
        Mout = Mt; ok = true;
        rep.note = sprintf('jitter repair (tries=%d)', t);
        rep.nTries = t;
        rep.finalJitter = j;
        return;
    end

    j = j * grow;
end

% Fail
Mout = M;
ok = false;
rep.note = sprintf('failed after %d tries (final jitter=%.4g m)', maxTries, j);
rep.nTries = maxTries;
rep.finalJitter = j;

end

function M = clampToLmax(M, Lmax)
r = vecnorm(M,2,2);
mask = r > Lmax;
if any(mask)
    M(mask,:) = M(mask,:) .* ((Lmax ./ r(mask)) * 0.999);
end
end

function v = getField(S, name, defaultVal)
if isstruct(S) && isfield(S,name)
    v = S.(name);
else
    v = defaultVal;
end
end