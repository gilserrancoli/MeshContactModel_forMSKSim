% % List of movements
% clear all;
% close all;
% clc
% % movements = { 'ngait_og1', 'bouncy4', 'mtpgait3', 'ngait_tm_fast1', 'ngait_og5', 'bouncy7','mtpgait9','ngait_tm_set1'};
% movements = { 'ngait_og1'};
% 
% % Threshold values
% % thresholds = [2, 1, 0.5];
% thresholds=[2];
% femfaces_values = [171];
% 
% % for i = 1:length(movements)
% %     nametrial_id = movements{i};
% %     for j = 1:length(thresholds) 
% %         for useReducedPolynomials = [0]
% %             err_poly = thresholds(j);
% %                 for femfaces = femfaces_values  % Iterate over femfaces values (171 and 258)
% % 
% %                 % Create Options structure and assign femfaces value
% %                 Options.nfacesFem = femfaces;
% % 
% %                 TrackSim_3D_GC_v2(nametrial_id,useReducedPolynomials,err_poly,Options);
% %                 end
% % 
% %         end
% %     end
% % end
% for i = 1:length(movements)
%     nametrial_id = movements{i};
%     for j = 1:length(thresholds) 
%         for useReducedPolynomials = [1]
% 
%             err_poly = thresholds(j); % Keep err_poly assignment before the femfaces loop
% 
%             for femfaces = femfaces_values  % Iterate over femfaces values (171 and 258)
% 
%                 % Create Options structure and assign femfaces value
%                 Options.nfacesFem = femfaces;  %  Add femfaces to Options struct
%                 Options.useReducedPolynomials = useReducedPolynomials;
%                 Options.err_poly = err_poly;
% 
%                 % Call the function with Options struct
%                 TrackSim_3D_GC_v2(nametrial_id, useReducedPolynomials, err_poly, Options);
% 
%             end
%         end
%     end
% end
% clear all;
% close all;
% clc
% 
% movements = {'mtpgait9'};
% thresholds=[1];
% femfaces_values = [188];
% % 
% % Define ranges of weights to try
% % Qs_weights_to_try  = [50,40,30];
% % KCF_weights_to_try = [20,25,30];
% % GRF_weights_to_try = [20,25,30];
% % MsAc_weights_to_try =10;
% % % Initialize result storage
% % results = [];
% % trial_number = 0;
% % 
% % 
% % for i = 1:length(movements)
% %     nametrial_id = movements{i};
% %     for j = 1:length(thresholds) 
% %         for useReducedPolynomials = [1]
% %             err_poly = thresholds(j);
% %             for femfaces = femfaces_values
% %                 % Set options
% %                 Options.nfacesFem = femfaces;
% %                 Options.useReducedPolynomials = useReducedPolynomials;
% %                 Options.err_poly = err_poly;
% % 
% %                 % Create output folder
% %                 outputFolder = fullfile('Results', nametrial_id);
% %                 if ~exist(outputFolder, 'dir')
% %                     mkdir(outputFolder);
% %                 end
% % 
% % 
% %                 % Tuning loop
% %                 for w1 = 1:length(Qs_weights_to_try)
% %                     for w2 = 1:length(KCF_weights_to_try)
% %                         for w3 = 1:length(GRF_weights_to_try)
% %                             for w4 = 1:length(MsAc_weights_to_try)
% %                                 trial_number = trial_number + 1;
% %                                 % Check skip condition
% %                                 if Qs_weights_to_try(w1) == 20 && ...
% %                                    KCF_weights_to_try(w2) == 30 && ...
% %                                    GRF_weights_to_try(w3) == 30
% %                                     continue; % Skip this combination
% %                                 end
% % 
% %                                 % Set Weights
% %                                 W.Qs = Qs_weights_to_try(w1);
% %                                 W.KCF = KCF_weights_to_try(w2);
% %                                 W.GRF = GRF_weights_to_try(w3);
% %                                 W.a = MsAc_weights_to_try(w4);
% %                                 W.Qdots = 10;
% %                                 W.GRM = 10;
% %                                 W.ID_act = 0;
% %                                 W.minPelvisRes = 2;
% %                                 W.u = 0.03;
% %                                 W.u_qd2dot = 0.003;
% %                                 W.u_vA = 0.52;
% % 
% %                                 % Unique save name
% %                                 savename_suffix = sprintf('_Qs%d_KCF%d_a%d_GRF%d_T%d', ...
% %                                     W.Qs, W.KCF, W.a, W.GRF, trial_number);
% % 
% %                                 % Run simulation
% %                                 Results_3D = TrackSim_3D_GC_v2(nametrial_id, useReducedPolynomials, err_poly, Options, W, savename_suffix);
% % 
% %                                 % Save simulation output
% %                                 save_filename = fullfile(outputFolder, ['Result' savename_suffix '.mat']);
% %                                 save(save_filename, 'Results_3D');
% % 
% %                                 % Extract tracking signals
% %                                 exp_hip_flex = Results_3D.NMesh_50.Qs_toTrack(:,10);
% %                                 sim_hip_flex = Results_3D.Simulated.Qs_opt(:,10);
% %                                 exp_hip_add = Results_3D.NMesh_50.Qs_toTrack(:,11);
% %                                 sim_hip_add = Results_3D.Simulated.Qs_opt(:,11);
% %                                 exp_hip_rot = Results_3D.NMesh_50.Qs_toTrack(:,12);
% %                                 sim_hip_rot = Results_3D.Simulated.Qs_opt(:,12);
% %                                 exp_knee_flex = Results_3D.NMesh_50.Qs_toTrack(:,14);
% %                                 sim_knee_flex = Results_3D.Simulated.Qs_opt(:,14);
% % 
% %                                 % Calculate RMSE and R2
% %                                 rmse = @(y, yhat) sqrt(mean((y - yhat).^2));
% %                                 r2 = @(y, yhat) 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
% % 
% %                                 % Store metrics
% %                                 results(trial_number).TrialID = trial_number;
% %                                 results(trial_number).W_Qs = W.Qs;
% %                                 results(trial_number).W_KCF = W.KCF;
% %                                 results(trial_number).W_GRF = W.GRF;
% %                                 results(trial_number).W_a = W.a;
% %                                 results(trial_number).RMSE_HipFlex = rmse(exp_hip_flex, sim_hip_flex);
% %                                 results(trial_number).R2_HipFlex = r2(exp_hip_flex, sim_hip_flex);
% %                                 results(trial_number).RMSE_HipAdd = rmse(exp_hip_add, sim_hip_add);
% %                                 results(trial_number).R2_HipAdd = r2(exp_hip_add, sim_hip_add);
% %                                 results(trial_number).RMSE_HipRot = rmse(exp_hip_rot, sim_hip_rot);
% %                                 results(trial_number).R2_HipRot = r2(exp_hip_rot, sim_hip_rot);
% %                                 results(trial_number).RMSE_Knee = rmse(exp_knee_flex, sim_knee_flex);
% %                                 results(trial_number).R2_Knee = r2(exp_knee_flex, sim_knee_flex);
% %                                 results(trial_number).SaveName = save_filename;
% %                             end
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % end
% % 
% % % After all simulations
% % ResultsTable = struct2table(results);
% % disp(ResultsTable);
% % sortedTable = sortrows(ResultsTable, {'RMSE_HipFlex','RMSE_HipAdd','RMSE_HipRot','RMSE_Knee'}, 'ascend');
% % % disp('Top 50 Best Trials:');
% % % disp(sortedTable(1:50,:));
% % Define specific weight combinations
% % weight_combinations = [
% %     20, 30, 25;
% %     30, 30, 40;
% %     30, 35, 40;
% %     40, 25, 30
% % ];
% 
% weight_combinations = [
%     20, 40, 25
%     ]; %20, 50, 25
% MsAc_weights_to_try = 10;  % Assuming this is constant
% trial_number = 0;
% 
% for i = 1:length(movements)
%     nametrial_id = movements{i};
%     for j = 1:length(thresholds)
%         for useReducedPolynomials = [1]
%             err_poly = thresholds(j);
%             for femfaces = femfaces_values
%                 % Set options
%                 Options.nfacesFem = femfaces;
%                 Options.useReducedPolynomials = useReducedPolynomials;
%                 Options.err_poly = err_poly;
% 
%                 % Create output folder
%                 outputFolder = fullfile('Results', nametrial_id);
%                 if ~exist(outputFolder, 'dir')
%                     mkdir(outputFolder);
%                 end
% 
%                 % Loop through defined weight combinations
%                 for k = 1:size(weight_combinations, 1)
%                     trial_number = trial_number + 1;
% 
%                     W.Qs  = weight_combinations(k, 1);
%                     W.KCF = weight_combinations(k, 2);
%                     W.GRF = weight_combinations(k, 3);
%                     W.a   = MsAc_weights_to_try;
% 
%                     W.Qdots = 10;
%                     W.GRM = 10;
%                     W.ID_act = 0;
%                     W.minPelvisRes = 0.2;
%                     W.u = 0.03;
%                     W.u_qd2dot = 0.003;
%                     W.u_qd2dot_kneesecdof=50;
%                     W.u_vA = 0.52;
% 
%                     % Unique save name
%                     savename_suffix = sprintf('_Qs%d_KCF%d_a%d_GRF%d_T%d', ...
%                         W.Qs, W.KCF, W.a, W.GRF, trial_number);
% 
%                     % Run simulation
%                     Results_3D = TrackSim_3D_GC_v2(nametrial_id, useReducedPolynomials, err_poly, Options, W, savename_suffix);
% 
%                     % Save simulation output
%                     save_filename = fullfile(outputFolder, ['Result' savename_suffix '.mat']);
%                     save(save_filename, 'Results_3D');
% 
%                     % Extract and compute RMSE & R2
%                     exp_hip_flex = Results_3D.NMesh_40.Qs_toTrack(:,10);
%                     sim_hip_flex = Results_3D.Simulated.Qs_opt(:,10);
%                     exp_hip_add  = Results_3D.NMesh_40.Qs_toTrack(:,11);
%                     sim_hip_add  = Results_3D.Simulated.Qs_opt(:,11);
%                     exp_hip_rot  = Results_3D.NMesh_40.Qs_toTrack(:,12);
%                     sim_hip_rot  = Results_3D.Simulated.Qs_opt(:,12);
%                     exp_knee_flex = Results_3D.NMesh_40.Qs_toTrack(:,14);
%                     sim_knee_flex = Results_3D.Simulated.Qs_opt(:,14);
% 
%                     rmse = @(y, yhat) sqrt(mean((y - yhat).^2));
%                     r2   = @(y, yhat) 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
% 
%                     results(trial_number).TrialID = trial_number;
%                     results(trial_number).W_Qs = W.Qs;
%                     results(trial_number).W_KCF = W.KCF;
%                     results(trial_number).W_GRF = W.GRF;
%                     results(trial_number).W_a = W.a;
%                     results(trial_number).RMSE_HipFlex = rmse(exp_hip_flex, sim_hip_flex);
%                     results(trial_number).R2_HipFlex   = r2(exp_hip_flex, sim_hip_flex);
%                     results(trial_number).RMSE_HipAdd  = rmse(exp_hip_add, sim_hip_add);
%                     results(trial_number).R2_HipAdd    = r2(exp_hip_add, sim_hip_add);
%                     results(trial_number).RMSE_HipRot  = rmse(exp_hip_rot, sim_hip_rot);
%                     results(trial_number).R2_HipRot    = r2(exp_hip_rot, sim_hip_rot);
%                     results(trial_number).RMSE_Knee    = rmse(exp_knee_flex, sim_knee_flex);
%                     results(trial_number).R2_Knee      = r2(exp_knee_flex, sim_knee_flex);
%                     results(trial_number).SaveName     = save_filename;
%                 end
%             end
%         end
%     end
% end
% 
% % Convert results to table and sort
% ResultsTable = struct2table(results);
% disp(ResultsTable);
% sortedTable = sortrows(ResultsTable, {'RMSE_HipFlex','RMSE_HipAdd','RMSE_HipRot','RMSE_Knee'}, 'ascend');


clear all; close all; clc

% movements = { 'ngait_og1', 'bouncy4', 'mtpgait3', 'ngait_tm_fast1', 'ngait_og5', 'bouncy7','mtpgait9','ngait_tm_set1'};
movements = { 'ngait_og1'};

thresholds = [1];

% Faces (paired by index)
femfaces_values = [342];
tifaces_values  = [100];
% Damping factor
damping_Pairs = [0.001 1000];  % moment damping for knee add/rot
% Parameter grids
kmax_list   = {1e4};   % Options.kInmaxpen 1e4 1e3 or max "Max"
kpress_list = [1e4];   % Options.kInpress 1e3 1e4 1e5  best 1e4
kcheck_list = [1e2];   % Options.kInCheckContacts 1e2 1e4 best 1e3
rad_list    = [1];             % Options.rad4Pairs (use [1 0.5] to sweep both)


% Use contact model from ANN (1) or from collision detection algorithm (0)
Options.useANNforKneeCont=0;
    Options.prunedANN=0;

if Options.useANNforKneeCont==1
    if Options.prunedANN==1
        suffix_typeKneeCont='_ANN_pruned';
    else
        suffix_typeKneeCont='_ANN';
    end
else
    suffix_typeKneeCont='';
end

% Optimization weights (you can keep testing weights, but we won't put them in filenames)
weight_combinations = [20, 40, 25];
MsAc_weights_to_try = 10;

trial_number = 0;
results = [];

for i = 1:length(movements)
    nametrial_id = movements{i};
    outputFolder = fullfile('Results', nametrial_id);
    if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

    for j = 1:length(thresholds)
        err_poly = thresholds(j);

        for useReducedPolynomials = [1]
            for p = 1:numel(femfaces_values)
                % Faces
                Options.nfacesFem = femfaces_values(p);
                Options.nfacesTib = tifaces_values(p);
                Options.useReducedPolynomials = useReducedPolynomials;
                Options.err_poly = err_poly;

                fprintf('\n==== Faces: Femur %d, Tibia %d | Move: %s ====\n', ...
                    Options.nfacesFem, Options.nfacesTib, nametrial_id);

                % Sweep model 
                for rad = rad_list
                    Options.rad4Pairs = rad;

                    for ik = 1:numel(kmax_list)
                        kmax = kmax_list{ik};
                        % Keep MellowMax;  kmax = "Max"
                        % Always use the MellowMax DLL family
                        Options.maxsmoothness = 'MellowMax';

                        if (isstring(kmax) || ischar(kmax))
                            if strcmpi(kmax, "Max")
                                Options.kInmaxpen = "Max";   % string token  kmaxMax 
                            else
                                Options.kInmaxpen = str2double(kmax); %  "1e3" -> 1000
                            end
                        else
                            Options.kInmaxpen = kmax; % already numeric
                        end

                        for kpress = kpress_list
                            Options.kInpress = kpress;

                            for kcheck = kcheck_list
                                Options.kInCheckContacts = kcheck;
                                Options.maxsmoothness = Options.maxsmoothness; % matches your DLLs
                                for dp = 1:size(damping_Pairs,1)
                                %% 
                                    Options.dampM = damping_Pairs(dp,1);
                                    Options.dampF = damping_Pairs(dp,2);

                                        % Build a case-based suffix 
                                        % map rad to text
                                        if Options.rad4Pairs == 0.5
                                           radStr = 'rad05';
                                        else
                                            radStr = 'rad1';
                                        end
        
                                        % string for kmax
                                        if (isstring(Options.kInmaxpen) || ischar(Options.kInmaxpen)) && strcmpi(Options.kInmaxpen,"Max")
                                            kmaxStr = 'Max';
                                        else
                                            kmaxStr = sprintf('%g', Options.kInmaxpen);
                                            kmaxStr = strrep(kmaxStr,'+','');
                                        end
                                        kpressStr = sprintf('%g', Options.kInpress);        kpressStr = strrep(kpressStr,'+','');
                                        kcheckStr = sprintf('%g', Options.kInCheckContacts);kcheckStr = strrep(kcheckStr,'+','');

                                         % damping strings for filenames 
                                        dMStr = sprintf('%g', Options.dampM);  dMStr = strrep(dMStr,'+',''); dMStr = strrep(dMStr,'.','p');
                                        dFStr = sprintf('%g', Options.dampF);  dFStr = strrep(dFStr,'+',''); dFStr = strrep(dFStr,'.','p');
                            
                                        savename_suffix = sprintf('_Fem%d_Tib%d_kmax%s_kpress%s_kcheck%s_%s_dM%s_dF%s%s', ...
                                            Options.nfacesFem, Options.nfacesTib, kmaxStr, kpressStr, kcheckStr, radStr, dMStr, dFStr, suffix_typeKneeCont);
        
                                        % Loop weight  
                                        for k = 1:size(weight_combinations, 1)
                                            trial_number = trial_number + 1;
        
                                            W.Qs  = weight_combinations(k, 1);
                                            W.KCF = weight_combinations(k, 2);
                                            W.GRF = weight_combinations(k, 3);
                                            W.a   = MsAc_weights_to_try;
        
                                            W.Qdots = 10; W.GRM = 10; W.ID_act = 0;
                                            W.minPelvisRes = 0.2; W.u = 0.03; W.u_qd2dot = 0.003;
                                            W.u_qd2dot_kneesecdof = 50; W.u_vA = 0.52;
        
                                            % Run simulation 
                                            Results_3D = TrackSim_3D_GC_v2( ...
                                                nametrial_id, useReducedPolynomials, err_poly, Options, W, savename_suffix);
        
                                            % Save using the case suffix (no weights)
                                            save_filename = fullfile(outputFolder, ['Result' savename_suffix '.mat']);
                                            save(save_filename, 'Results_3D');
        
        
                                            if isfield(Results_3D, 'NMesh_50')
                                                NMesh = Results_3D.NMesh_50;
                                            elseif isfield(Results_3D, 'NMesh_40')
                                                NMesh = Results_3D.NMesh_40;
                                            else
                                                error('Neither NMesh_50 nor NMesh_40 found in Results_3D.');
                                            end
        
                                            exp_hip_flex  = NMesh.Qs_toTrack(:,10);
                                            sim_hip_flex  = Results_3D.Simulated.Qs_opt(:,10);
                                            exp_hip_add   = NMesh.Qs_toTrack(:,11);
                                            sim_hip_add   = Results_3D.Simulated.Qs_opt(:,11);
                                            exp_hip_rot   = NMesh.Qs_toTrack(:,12);
                                            sim_hip_rot   = Results_3D.Simulated.Qs_opt(:,12);
                                            exp_knee_flex = NMesh.Qs_toTrack(:,14);
                                            sim_knee_flex = Results_3D.Simulated.Qs_opt(:,14);
        
                                            rmse = @(y, yhat) sqrt(mean((y - yhat).^2));
                                            r2   = @(y, yhat) 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
        
                                            r.TrialID = trial_number;
                                            r.Movement = nametrial_id;
                                            r.FemurFaces = Options.nfacesFem;
                                            r.TibiaFaces = Options.nfacesTib;
                                            r.kInmaxpen = Options.kInmaxpen;
                                            r.kInpress  = Options.kInpress;
                                            r.kInCheckContacts = Options.kInCheckContacts;
                                            r.rad4Pairs = Options.rad4Pairs;
                                            r.W_Qs = W.Qs; r.W_KCF = W.KCF; r.W_GRF = W.GRF; r.W_a = W.a;
        
                                            r.RMSE_HipFlex = rmse(exp_hip_flex,  sim_hip_flex);
                                            r.R2_HipFlex   = r2(  exp_hip_flex,  sim_hip_flex);
                                            r.RMSE_HipAdd  = rmse(exp_hip_add,   sim_hip_add);
                                            r.R2_HipAdd    = r2(  exp_hip_add,   sim_hip_add);
                                            r.RMSE_HipRot  = rmse(exp_hip_rot,   sim_hip_rot);
                                            r.R2_HipRot    = r2(  exp_hip_rot,   sim_hip_rot);
                                            r.RMSE_Knee    = rmse(exp_knee_flex, sim_knee_flex);
                                            r.R2_Knee      = r2(  exp_knee_flex, sim_knee_flex);
        
                                            r.SaveName     = save_filename;
                                            % results(trial_number) = r; 
                                        
                                         end % weights
                                    end % damping
                                end % kcheck
                            end % kpress
                        end % kmax
                    end % rad
                end % faces
            end % useReducedPolynomials
        end % thresholds
end % movements

% ---- Results table ----
% ResultsTable = struct2table(results);
% disp(ResultsTable);
% sortedTable = sortrows(ResultsTable, {'RMSE_HipFlex','RMSE_HipAdd','RMSE_HipRot','RMSE_Knee'}, 'ascend');




% clear; close all; clc
% 
% %% Inputs 
% movements   = {'ngait_og1', 'bouncy4', 'mtpgait3', 'ngait_tm_fast1', 'bouncy7','mtpgait9','ngait_tm_set1'};      % add more as needed
% thresholds  = [1];
% 
% % ---- Neutral (baseline) values ----
% neutral.femfaces = 258;
% neutral.tifaces  = 75;
% neutral.kmax     = 1e4;           % Options.kInmaxpen
% neutral.kpress   = 1e4;           % Options.kInpress
% neutral.kcheck   = 1e3;           % Options.kInCheckContacts
% neutral.rad      = 1;             % Options.rad4Pairs
% 
% % ---- Sweep lists ----
% faces_pairs   = [258 75; 342 100; 171 49];               % paired faces sweep
% kpress_list   = [1e3 1e4 1e5];                           % sweep this with others neutral
% kmax_list     = {1e4, 2e4, 'Max'};                       % includes 'Max'
% kcheck_list   = [1e2 1e3 1e4];
% rad_list      = [1];                                     % add 0.5 if you want a rad sweep
% 
% % ---- Weights (unchanged) ----
% weight_combinations = [20, 40, 25];   % [W.Qs, W.KCF, W.GRF]
% MsAc_weights_to_try = 10;
% 
% %% Build experiment scenarios (OFAT)
% % Each scenario defines which single field(s) to sweep; others stay neutral
% experiments = {};
% 
% % 0) Baseline only
% experiments{end+1} = struct('name','baseline', 'type','baseline');
% 
% % 1) Faces sweep (paired fem/tib)
% experiments{end+1} = struct('name','faces', 'type','faces', ...
%     'values', arrayfun(@(i) struct('fem',faces_pairs(i,1),'tib',faces_pairs(i,2)), ...
%                        1:size(faces_pairs,1)));
% 
% % 2) kpress sweep
% experiments{end+1} = struct('name','kpress', 'type','kpress', 'values', num2cell(kpress_list));
% 
% % 3) kmax sweep
% experiments{end+1} = struct('name','kmax', 'type','kmax', 'values', kmax_list);
% 
% % 4) kcheck sweep
% experiments{end+1} = struct('name','kcheck', 'type','kcheck', 'values', num2cell(kcheck_list));
% 
% % 5) rad sweep (only if more than one)
% if numel(rad_list) > 1 || rad_list(1) ~= neutral.rad
%     experiments{end+1} = struct('name','rad', 'type','rad', 'values', num2cell(rad_list));
% end
% 
% %% Run
% trial_number = 0;
% 
% for i = 1:numel(movements)
%     nametrial_id = movements{i};
%     outputFolder = fullfile('Results', nametrial_id);
%     if ~exist(outputFolder,'dir'); mkdir(outputFolder); end
% 
%     for j = 1:numel(thresholds)
%         err_poly = thresholds(j);
% 
%         for useReducedPolynomials = 1  % set [0 1] if you also want to sweep this
%             for ex = 1:numel(experiments)
% 
%                 expt = experiments{ex};
% 
%                 % Determine list of variants for this experiment
%                 switch expt.type
%                     case 'baseline'
%                         variant_list = {[]};   % single run with neutral values
%                     otherwise
%                         variant_list = expt.values;
%                 end
% 
%                 for v = 1:numel(variant_list)
% 
%                     % ---- Start from neutral Options each time ----
%                     Options = struct();
%                     Options.nfacesFem = neutral.femfaces;
%                     Options.nfacesTib = neutral.tifaces;
%                     Options.kInmaxpen = neutral.kmax;
%                     Options.kInpress  = neutral.kpress;
%                     Options.kInCheckContacts = neutral.kcheck;
%                     Options.rad4Pairs = neutral.rad;
% 
%                     Options.useReducedPolynomials = useReducedPolynomials;
%                     Options.err_poly = err_poly;
%                     Options.maxsmoothness = 'MellowMax';  % keep your DLL family
% 
%                     % ---- Override only the parameter under test ----
%                     suffix_bits = {};
% 
%                     switch expt.type
%                         case 'baseline'
%                             % no overrides
%                         case 'faces'
%                             pair = variant_list{v};
%                             Options.nfacesFem = pair.fem;
%                             Options.nfacesTib = pair.tib;
%                             suffix_bits{end+1} = sprintf('Fem%d_Tib%d', pair.fem, pair.tib);
% 
%                         case 'kpress'
%                             Options.kInpress = variant_list{v};
%                             suffix_bits{end+1} = sprintf('kpress%g', Options.kInpress);
% 
%                         case 'kmax'
%                             kv = variant_list{v};
%                             if ischar(kv) || (isstring(kv) && strcmpi(kv,"Max"))
%                                 Options.kInmaxpen = "Max";
%                                 suffix_bits{end+1} = 'kmaxMax';
%                             else
%                                 Options.kInmaxpen = kv;
%                                 suffix_bits{end+1} = sprintf('kmax%g', Options.kInmaxpen);
%                             end
% 
%                         case 'kcheck'
%                             Options.kInCheckContacts = variant_list{v};
%                             suffix_bits{end+1} = sprintf('kcheck%g', Options.kInCheckContacts);
% 
%                         case 'rad'
%                             Options.rad4Pairs = variant_list{v};
%                             if Options.rad4Pairs == 0.5, suffix_bits{end+1} = 'rad05';
%                             else,                      suffix_bits{end+1} = 'rad1';
%                             end
%                     end
% 
%                     % Always include faces and rad in suffix for traceability
%                     faceStr = sprintf('Fem%d_Tib%d', Options.nfacesFem, Options.nfacesTib);
%                     if Options.rad4Pairs == 0.5, radStr = 'rad05'; else, radStr = 'rad1'; end
% 
%                     % kmax string
%                     if (ischar(Options.kInmaxpen) || isstring(Options.kInmaxpen)) && strcmpi(Options.kInmaxpen,"Max")
%                         kmaxStr = 'Max';
%                     else
%                         kmaxStr = sprintf('%g', Options.kInmaxpen);
%                     end
% 
%                     base_suffix = sprintf('%s_kmax%s_kpress%g_kcheck%g_%s', ...
%                         faceStr, kmaxStr, Options.kInpress, Options.kInCheckContacts, radStr);
% 
%                     % If this is a sweep, append only the changing piece(s) too (nice for filtering)
%                     if ~isempty(suffix_bits)
%                         savename_suffix = ['_' base_suffix '_' strjoin(suffix_bits,'_')];
%                     else
%                         savename_suffix = ['_' base_suffix];
%                     end
% 
%                     % ---- Run weights for this case ----
%                     for k = 1:size(weight_combinations,1)
%                         trial_number = trial_number + 1;
% 
%                         W.Qs  = weight_combinations(k,1);
%                         W.KCF = weight_combinations(k,2);
%                         W.GRF = weight_combinations(k,3);
%                         W.a   = MsAc_weights_to_try;
% 
%                         % keep your extra weights
%                         W.Qdots = 10; W.GRM = 10; W.ID_act = 0;
%                         W.minPelvisRes = 0.2; W.u = 0.03; W.u_qd2dot = 0.003;
%                         W.u_qd2dot_kneesecdof = 50; W.u_vA = 0.52;
% 
%                         fprintf('\n==== [%s] Move:%s | %s ====\n', ...
%                             expt.name, nametrial_id, savename_suffix(2:end));
% 
%                         Results_3D = TrackSim_3D_GC_v2( ...
%                             nametrial_id, useReducedPolynomials, err_poly, Options, W, savename_suffix);
% 
%                         save_filename = fullfile(outputFolder, ['Result' savename_suffix '.mat']);
%                         save(save_filename, 'Results_3D');
% 
%                         % --- Metrics (unchanged) ---
%                         if isfield(Results_3D,'NMesh_50')
%                             NMesh = Results_3D.NMesh_50;
%                         elseif isfield(Results_3D,'NMesh_40')
%                             NMesh = Results_3D.NMesh_40;
%                         else
%                             error('Neither NMesh_50 nor NMesh_40 found in Results_3D.');
%                         end
% 
%                         exp_hip_flex  = NMesh.Qs_toTrack(:,10);
%                         sim_hip_flex  = Results_3D.Simulated.Qs_opt(:,10);
%                         exp_hip_add   = NMesh.Qs_toTrack(:,11);
%                         sim_hip_add   = Results_3D.Simulated.Qs_opt(:,11);
%                         exp_hip_rot   = NMesh.Qs_toTrack(:,12);
%                         sim_hip_rot   = Results_3D.Simulated.Qs_opt(:,12);
%                         exp_knee_flex = NMesh.Qs_toTrack(:,14);
%                         sim_knee_flex = Results_3D.Simulated.Qs_opt(:,14);
% 
%                         rmse = @(y,yhat) sqrt(mean((y - yhat).^2));
%                         r2   = @(y,yhat) 1 - sum((y - yhat).^2)/sum((y - mean(y)).^2);
% 
%                         r.TrialID = trial_number;
%                         r.Movement = nametrial_id;
%                         r.FemurFaces = Options.nfacesFem;
%                         r.TibiaFaces = Options.nfacesTib;
%                         r.kInmaxpen = Options.kInmaxpen;
%                         r.kInpress  = Options.kInpress;
%                         r.kInCheckContacts = Options.kInCheckContacts;
%                         r.rad4Pairs = Options.rad4Pairs;
%                         r.W_Qs = W.Qs; r.W_KCF = W.KCF; r.W_GRF = W.GRF; r.W_a = W.a;
% 
%                         r.RMSE_HipFlex = rmse(exp_hip_flex,  sim_hip_flex);
%                         r.R2_HipFlex   = r2(  exp_hip_flex,  sim_hip_flex);
%                         r.RMSE_HipAdd  = rmse(exp_hip_add,   sim_hip_add);
%                         r.R2_HipAdd    = r2(  exp_hip_add,   sim_hip_add);
%                         r.RMSE_HipRot  = rmse(exp_hip_rot,   sim_hip_rot);
%                         r.R2_HipRot    = r2(  exp_hip_rot,   sim_hip_rot);
%                         r.RMSE_Knee    = rmse(exp_knee_flex, sim_knee_flex);
%                         r.R2_Knee      = r2(  exp_knee_flex, sim_knee_flex);
% 
%                         r.SaveName = save_filename;
%                         % TODO: append r to a struct array or table if you want a summary file.
%                     end
%                 end
%             end
%         end
%     end
% end
