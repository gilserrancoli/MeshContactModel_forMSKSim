%% Code to Run all tracking simulations for the sensitivity analysis of the 
% smooth mesh-based knee contact pressure model
% Authors: Mohanad Harba and Gil Serrancolí
% Date: 01-09-2026

clear all; close all; clc

%List of movements
movements = { 'ngait_og1', 'bouncy4', 'mtpgait3', 'ngait_tm_fast1', 'ngait_og5', 'bouncy7','mtpgait9','ngait_tm_set1'};


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

