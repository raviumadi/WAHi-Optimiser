%% economy_responsiveness_table.m
clear; clc;

%% -------- INPUT --------
posCSV = '../results/summary_stats_pos.csv';
angCSV = '../results/summary_stats_ang.csv';
outCSV = '../results/economy_responsiveness_final.csv';

%% -------- LOAD --------
Tpos = readtable(posCSV);
Tang = readtable(angCSV);

%% -------- KEY FOR JOIN --------
keyPos = strcat(string(Tpos.startType), "|", string(Tpos.condition), "|N=", string(Tpos.N));
keyAng = strcat(string(Tang.startType), "|", string(Tang.condition), "|N=", string(Tang.N));

Tpos.key = keyPos;
Tang.key = keyAng;

%% -------- ECONOMY METRICS --------
Tpos.delta_pos = Tpos.pass_after - Tpos.pass_before;
Tpos.econ_pos  = Tpos.delta_pos ./ Tpos.N;

Tang.delta_ang = Tang.pass_after - Tang.pass_before;
Tang.econ_ang  = Tang.delta_ang ./ Tang.N;

%% -------- REDUCE & JOIN --------
A = Tpos(:,{'key','startType','condition','N','delta_pos','econ_pos'});
B = Tang(:,{'key','delta_ang','econ_ang'});

T = innerjoin(A,B,'Keys','key');

%% -------- RANKINGS --------
[~,ordPos] = sort(T.econ_pos,'descend');
rank_pos = zeros(height(T),1);
rank_pos(ordPos) = 1:height(T);

[~,ordAng] = sort(T.econ_ang,'descend');
rank_ang = zeros(height(T),1);
rank_ang(ordAng) = 1:height(T);

T.rank_pos = rank_pos;
T.rank_ang = rank_ang;

%% -------- FINAL SORT (BY POSITION RANK) --------
T = sortrows(T,'rank_pos','ascend');

%% -------- FINAL TABLE --------
FinalTable = T(:,{ ...
    'startType','condition','N', ...
    'delta_pos','econ_pos','rank_pos', ...
    'delta_ang','econ_ang','rank_ang' ...
});

%% -------- WRITE --------
writetable(FinalTable,outCSV);

disp('Final economy–responsiveness table (sorted by position ranking):');
disp(FinalTable);
fprintf('\nWrote: %s\n', outCSV);