%% run_iwah_experiment_batch.m
% Batch experiment runner for iWAH optimiser
% - runs multiple CSV geometries + random-start geometries at N={4,6,8,10,12}
% - repeats each condition K times
% - optional parallelism using parfor
%
% Output: many .mat files, each containing H + metaRun info

clear; clc;

%% ===================== EXPERIMENT SETTINGS =====================
EXP = struct();
EXP.repeats = 5;
EXP.useParfor = true;        % set false to debug serially
EXP.outDir = fullfile(pwd, "../data/runs/runsXX"); % <-- name your run destination folder
if ~exist(EXP.outDir,'dir'); mkdir(EXP.outDir); end

% CSV configs to test (edit this list)
EXP.csvList = {
    fullfile("../configs","4mics_Pyramid.csv")
    fullfile("../configs","4mics_Tetrahedron.csv")
    fullfile("../configs","6mics_Octahedron.csv")
    fullfile("../configs","8mics_double_tetrahedorn.csv")
    fullfile("../configs","10mics_pentagonal_antiprism.csv")
    fullfile("../configs","12mics_icosahedron.csv")
};

% Random-start N sweep
EXP.randomNs = [4 6 8 10 12];

% Optional: also run "CSV-derived N" (default true)
EXP.useCsvNativeN = true;    % if CSV has N rows, use that N

%% ===================== BASE PARAMETERS (YOUR CURRENT VALUES) =====================
P = struct();
P.N = 8;                 % overwritten per condition
P.Lmax = 0.50;
P.minPairDist = 0.125/2;

P.R = 2.0;
P.Rin = 0.30;
P.step = 0.250;

P.best_cm = 10;
P.tol_cm  = 5;
P.threshold_cm = P.best_cm + P.tol_cm;

P.fs = 192000;
P.d  = 5e-3;
P.f0 = 30000;
P.f1 = 70000;
P.snr_db = 40;
P.tail = 100;
P.callType = 'FM';
P.velocity = [0 0 0];

P.minDistToAnyMic = 0.08;
P.epsDist         = 0.03;
P.failPenalty_cm  = 1000;
P.c               = 343;

P.lambdaMean = 0.15;
P.targetPass = 0.95;

OPT = struct();
OPT.maxIters      = 25;
OPT.probesPerMic  = 14;
OPT.initStepFrac  = 0.25;
OPT.minStepFrac   = 0.03;
OPT.stallLimit    = 12;
OPT.shakeFrac     = 0.08;
OPT.improveEpsJ   = 1e-6;

Q = struct();
Q.frontTopLeft     = true;
Q.frontTopRight    = true;
Q.frontBottomLeft  = false;
Q.frontBottomRight = false;
Q.backTopLeft      = false;
Q.backTopRight     = false;
Q.backBottomLeft   = false;
Q.backBottomRight  = false;

CSV = struct();
CSV.hasHeader = false;
CSV.recenter = true;
CSV.enforceConstraints = true;

CSV.autoRepair = true;
CSV.repairMaxTries = 200;
CSV.repairJitter0 = 1e-3;
CSV.repairJitterGrowth = 1.15;
CSV.repairClampLmax = true;
CSV.repairPreferScale = true;
CSV.repairFallbackToRandom = false;

%% ===================== BUILD CONDITION LIST =====================
conds = {};

% (A) CSV conditions
for i = 1:numel(EXP.csvList)
    pth = EXP.csvList{i};
    [~,name,~] = fileparts(pth);
    cond = struct();
    cond.mode = "csv";
    cond.csvPath = pth;
    cond.configTag = string(name); % filename without .csv
    cond.useCsvNativeN = EXP.useCsvNativeN;
    conds{end+1} = cond; %#ok<SAGROW>
end

% (B) Random-start conditions for N sweep
for n = EXP.randomNs
    cond = struct();
    cond.mode = "random";
    cond.csvPath = "";
    cond.configTag = "randomStart";
    cond.N = n;
    conds{end+1} = cond; %#ok<SAGROW>
end

% Expand repeats
jobs = {};
for c = 1:numel(conds)
    for r = 1:EXP.repeats
        job = conds{c};
        job.rep = r;
        jobs{end+1} = job; %#ok<SAGROW>
    end
end

fprintf("Total runs: %d\n", numel(jobs));

%% ===================== OPTIONAL PARALLEL POOL =====================
if EXP.useParfor
    if isempty(gcp('nocreate'))
        parpool; % uses default settings
    end
end

%% ===================== RUN =====================
if EXP.useParfor
    parfor j = 1:numel(jobs)
        runOneJob(jobs{j}, P, OPT, Q, CSV, EXP.outDir);
    end
else
    for j = 1:numel(jobs)
        runOneJob(jobs{j}, P, OPT, Q, CSV, EXP.outDir);
    end
end

fprintf("\nBatch complete. Output folder:\n%s\n", EXP.outDir);

%% ===================== LOCAL: one job runner =====================
function runOneJob(job, P, OPT, Q, CSV, outDir)
    % Create a per-run seed (no fixed seed) but still stored for reproducibility
    seed = makeNonRepeatingSeed();
    rng(seed, 'twister');

    % Build per-run config structs
    CSVrun = CSV;
    if job.mode == "csv"
        CSVrun.useCSVStart = true;
        CSVrun.path = job.csvPath;
        % For CSV we allow N to be taken from file (inside optimiser)
        Prun = P;
    else
        CSVrun.useCSVStart = false;
        CSVrun.path = "";
        Prun = P;
        Prun.N = job.N;
    end

    % Run optimiser (no plotting)
    [H, runInfo] = iwah_run_optimiser_core(Prun, OPT, Q, CSVrun, seed);

    % Attach run identifiers for analysis later
    runInfo.configTag = job.configTag;
    runInfo.rep = job.rep;
    runInfo.mode = char(job.mode);

    if job.mode == "csv"
        [~,nm,~] = fileparts(job.csvPath);
        runInfo.csvFile = nm;
    else
        runInfo.csvFile = "";
    end

    % Filename: include configTag + N + rep + datetime stamp
    Nfinal = size(H.meta.final.M,1);
    ts = datestr(now,'yyyymmdd_HHMMSS');
    fn = sprintf("iwah_optimiser_run__%s__%dmics__rep%02d__%s.mat", ...
        runInfo.configTag, Nfinal, job.rep, ts);

    save(fullfile(outDir, fn), 'H', 'runInfo', '-v7.3');
end

function seed = makeNonRepeatingSeed()
    % robust-ish seed for parallel contexts
    t = now;                         % days
    frac = t - floor(t);
    seed = mod(uint32(frac * 1e9) + uint32(feature('getpid')) + uint32(randi(2^31-1)), 2^31-1);
    seed = double(seed);
end