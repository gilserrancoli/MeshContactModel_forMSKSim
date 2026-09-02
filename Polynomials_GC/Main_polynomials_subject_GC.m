% This function generates polynomials to approximate muscle-tendon lengths
% and moment arms. The code is from Wouter Aerts and is adapted to be 
% used with CasADi.
%
% Author: Antoine Falisse
%
% Datum: 03/04/2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% User inputs
runPolynomialfit = 1;
saveQdot = 0;
savePolynomials = 0;

%% Extract time and angles from dummy motion

pathmain = pwd;
name_dummymotion = 'dummy_motion.mot';
path_dummymotion = [pathmain,'\MuscleAnalysis\dummy_motion\'];
path_resultsMA = [pathmain,'\MuscleAnalysis\ResultsMA_subject_GC\'];

dummy_motion = importdata([path_dummymotion,name_dummymotion]);
% 15 dofs (mtp locked)
% Order of dofs: hip flex r, hip add r, hip rot r, knee add r, knee rot r, 
% knee flex r, knee tx r, knee ty r, knee tz r, ankle flex r, subtalar r, 
% hip flex l, hip add l, hip rot l, knee add l, knee rot l, knee flex l, 
% knee tx l, knee ty l, knee tz l, ankle flex l, lumbar
% ext, lumbar bend, lumbar rot, subtalar r, subtalar l
q = dummy_motion.data(:,[8:10, 13, 17:18 ,30:32]).*(pi/180);

% Generate random numbers between -1000 and 1000 (°/s) 
if saveQdot
    a = -1000;
    b = 1000;
    r1 = (b-a).*rand(size(q,1),1) + a;
    r2 = (b-a).*rand(size(q,1),1) + a;
    r3 = (b-a).*rand(size(q,1),1) + a;
    r4 = (b-a).*rand(size(q,1),1) + a;
    r5 = (b-a).*rand(size(q,1),1) + a;
    r6 = (b-a).*rand(size(q,1),1) + a;
    r7 = (b-a).*rand(size(q,1),1) + a;
    r8 = (b-a).*rand(size(q,1),1) + a;
    r9 = (b-a).*rand(size(q,1),1) + a;
    r = [r1,r2,r3,r4,r5,r6,r7,r8,r9];
    qdot = zeros(size(q));
    qdot = r.*(pi/180);
    dummy_qdot = qdot;
    save([path_dummymotion,'dummy_qdot.mat'],'dummy_qdot');
end
load([path_dummymotion,'dummy_qdot.mat']);
qdot = dummy_qdot(:,:);

%% Import data
% lMT
lMT = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_Length.sto']);
% hip flexion r
MA.hip.flex = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_hip_flexion.sto']);
% hip adduction r
MA.hip.add = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_hip_adduction.sto']);
% hip rotation r
MA.hip.rot = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_hip_rotation.sto']);
% knee adduction r 
MA.knee.add = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_adduction.sto']);
% knee rotation r 
MA.knee.rot = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_rotation.sto']);
% knee flexion r 
MA.knee.flex = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_flexion.sto']);
% knee tx r 
MA.knee.tx = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_tx.sto']);
% knee ty r 
MA.knee.ty = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_ty.sto']);
% knee tz r 
MA.knee.tz = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_knee_tz.sto']);
% ankle flexion r
MA.ankle.flex = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_ankle_angle.sto']);
% subtalar r
MA.sub = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_subtalar_angle.sto']);
% lumbar extension
MA.trunk.ext = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_lumbar_extension.sto']);
% lumbar bending
MA.trunk.ben = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_lumbar_bending.sto']);
% lumbar rotation
MA.trunk.rot = importdata([path_resultsMA,'JW All Body-scaled_MuscleAnalysis_MomentArm_lumbar_rotation.sto']);


%% Organize MuscleData
if runPolynomialfit
   %     MuscleData.dof_names = dummy_motion.colheaders([5:9,18,15:17]); 
    MuscleData.dof_names = dummy_motion.colheaders([8:18,30:32]); %right leg deofs and lumbar dofs 
    muscleNames_forLengths = {'addbrev','addlong','addmagProx','addmagMid',...
        'addmagDist','addmagIsch','bflh','bfsh','edl','ehl','fdl','fhl',...
        'gaslat','gasmed','gem','glmax1','glmax2','glmax3','glmed1',...
        'glmed2','glmed3','glmin1','glmin2','glmin3','grac','iliacus',...
        'pect','perbrev','perlong','pertert','piri','psoas','quadfem',...
        'recfem_withreallength','sart','semimem','semiten','soleus','tfl',...
        'tibant','tibpost','vasint_withreallength','vaslat_withreallength',...
        'vasmed_withreallength','ercspn_r','intobl_r','extobl_r',...
        'ercspn_l','intobl_l','extobl_l'};
    muscleNames_forMA = {'addbrev','addlong','addmagProx','addmagMid',...
        'addmagDist','addmagIsch','bflh','bfsh','edl','ehl','fdl','fhl',...
        'gaslat','gasmed','gem','glmax1','glmax2','glmax3','glmed1',...
        'glmed2','glmed3','glmin1','glmin2','glmin3','grac','iliacus',...
        'pect','perbrev','perlong','pertert','piri','psoas','quadfem',...
        'recfem','sart','semimem','semiten','soleus','tfl',...
        'tibant','tibpost','vasint','vaslat',...
        'vasmed','ercspn_r','intobl_r','extobl_r',...
        'ercspn_l','intobl_l','extobl_l'};
    MuscleData.muscle_names_forLengths = muscleNames_forLengths;
    MuscleData.muscle_names_forMA = muscleNames_forMA;
    for m = 1:length(muscleNames_forLengths)
        MuscleData.lMT(:,m)     = lMT.data(:,strcmp(lMT.colheaders,muscleNames_forLengths{m}));       % lMT    
        MuscleData.dM(:,m,1)    = MA.hip.flex.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));    % hip_flex
        MuscleData.dM(:,m,2)    = MA.hip.add.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % hip_add
        MuscleData.dM(:,m,3)    = MA.hip.rot.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % hip_rot
        MuscleData.dM(:,m,4)    = MA.knee.add.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));    % knee add
        MuscleData.dM(:,m,5)    = MA.knee.rot.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));    % knee rot
        MuscleData.dM(:,m,6)    = MA.knee.flex.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % knee flex
        MuscleData.dM(:,m,7)    = MA.knee.tx.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % knee tx
        MuscleData.dM(:,m,8)    = MA.knee.ty.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % knee ty
        MuscleData.dM(:,m,9)    = MA.knee.tz.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));     % knee tz
        MuscleData.dM(:,m,10)    = MA.ankle.flex.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));  % ankle
        MuscleData.dM(:,m,11)    = MA.sub.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));         % sub 
        MuscleData.dM(:,m,12)    = MA.trunk.ext.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % trunk ext
        MuscleData.dM(:,m,13)    = MA.trunk.ben.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % trunk ben
        MuscleData.dM(:,m,14)    = MA.trunk.rot.data(:,strcmp(lMT.colheaders,muscleNames_forMA{m}));   % trunk rot
    end
    MuscleData.q = q;
    MuscleData.qdot = qdot;
end

%% Call PolynomialFit
if runPolynomialfit
    [muscle_spanning_joint_INFO,MuscleInfo] = ...
        PolynomialFit_GC(MuscleData);
    %Call PolynomialFit for second knee dofs
    [MuscleInfo]=PolyFitKnee(MuscleData,MuscleInfo);

    if savePolynomials
        save MuscleData_subject_GC MuscleData
        save muscle_spanning_joint_INFO_subject_GC muscle_spanning_joint_INFO
        save MuscleInfo_subject_GC MuscleInfo
    end
end


%% Create CasADi functions
import casadi.*
% Order mobilities: hip_flex, hip_add, hip_rot, knee_angle, ankle-angle, 
load muscle_spanning_joint_INFO_subject_GC.mat
load MuscleInfo_subject_GC.mat
NMuscle = length(MuscleInfo.muscle);
q_leg_trunk = 9;
qin     = SX.sym('qin',1,q_leg_trunk);
qdotin  = SX.sym('qdotin',1,q_leg_trunk);
lMT     = SX(NMuscle,1);
vMT     = SX(NMuscle,1);
dM      = SX(NMuscle,q_leg_trunk);
for i=1:NMuscle     
    index_dof_crossing_withknee  = find(muscle_spanning_joint_INFO(i,:)==1);
    index_dof_crossing_withknee(index_dof_crossing_withknee==4)=[];
    index_dof_crossing_withknee(index_dof_crossing_withknee==5)=[];
    index_dof_crossing_withknee(index_dof_crossing_withknee==7)=[];
    index_dof_crossing_withknee(index_dof_crossing_withknee==8)=[];
    index_dof_crossing_withknee(index_dof_crossing_withknee==9)=[];
    index_dof_crossing_withoutknee= find(muscle_spanning_joint_INFO(i,[1:3 6 10:end])==1);
    order                        = MuscleInfo.muscle(i).order;
    [mat,diff_mat_q]             = n_art_mat_3_cas_SX(qin(1,index_dof_crossing_withoutknee),order);
    lMT(i,1)                     = mat*MuscleInfo.muscle(i).coeff;
    vMT(i,1)                     = 0;
    dM(i,1:q_leg_trunk)          = 0;
    nr_dof_crossing              = length(index_dof_crossing_withoutknee); 
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing_withoutknee(dof_nr)) = (-(diff_mat_q(:,dof_nr)))'*MuscleInfo.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing_withoutknee(dof_nr))*qdotin(1,index_dof_crossing_withoutknee(dof_nr)));
    end 
end
f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});

%% Check results
load MuscleData_subject_GC.mat
lMT_out_r = zeros(size(q,1),NMuscle);
vMT_out_r = zeros(size(q,1),NMuscle);
dM_out_r = zeros(size(q,1),NMuscle,q_leg_trunk);
for i = 1:size(q,1)
    [out1_r,out2_r,out3_r] = f_lMT_vMT_dM(MuscleData.q(i,:),MuscleData.qdot(i,:));
    lMT_out_r(i,:) = full(out1_r);
    vMT_out_r(i,:) = full(out2_r);
    dM_out_r(i,:,1) = full(out3_r(:,1));
    dM_out_r(i,:,2) = full(out3_r(:,2));
    dM_out_r(i,:,3) = full(out3_r(:,3));
    dM_out_r(i,:,4) = full(out3_r(:,4));
    dM_out_r(i,:,5) = full(out3_r(:,5));   
    dM_out_r(i,:,6) = full(out3_r(:,6));
    dM_out_r(i,:,7) = full(out3_r(:,7));
    dM_out_r(i,:,8) = full(out3_r(:,8));
    dM_out_r(i,:,9) = full(out3_r(:,9)); 
end

%% lMT
% right
figure()
subplot(4,4,1)
scatter(MuscleData.q(:,4),lMT_out_r(:,10)); hold on;
scatter(MuscleData.q(:,4),MuscleData.lMT(:,10));
xlabel('q knee');
title('BFSH');
subplot(4,4,2)
scatter(MuscleData.q(:,4),lMT_out_r(:,29)); hold on;
scatter(MuscleData.q(:,4),MuscleData.lMT(:,29));
xlabel('q knee');
title('VM');
subplot(4,4,3)
scatter(MuscleData.q(:,4),lMT_out_r(:,30)); hold on;
scatter(MuscleData.q(:,4),MuscleData.lMT(:,30));
xlabel('q knee');
title('VI');
subplot(4,4,4)
scatter(MuscleData.q(:,4),lMT_out_r(:,31)); hold on;
scatter(MuscleData.q(:,4),MuscleData.lMT(:,31));
xlabel('q knee');
title('VL');
subplot(4,4,5)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,34)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,34));
xlabel('q knee');
ylabel('q ankle');
title('GM');
subplot(4,4,6)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,35)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,35));
xlabel('q knee');
ylabel('q ankle');
title('GL');
subplot(4,4,7)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,36)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,36));
xlabel('q knee');
ylabel('q ankle');
title('GM');
subplot(4,4,8)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,37)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,37));
xlabel('q knee');
ylabel('q ankle');
title('GL');
subplot(4,4,9)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,38)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,38));
xlabel('q knee');
ylabel('q ankle');
title('GM');
subplot(4,4,10)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,39)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,39));
xlabel('q knee');
ylabel('q ankle');
title('GL');
subplot(4,4,11)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,40)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,40));
xlabel('q knee');
ylabel('q ankle');
title('GM');
subplot(4,4,12)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,41)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,41));
xlabel('q knee');
ylabel('q ankle');
title('GL');
subplot(4,4,13)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,42)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,42));
xlabel('q knee');
ylabel('q ankle');
title('GM');
subplot(4,4,14)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out_r(:,43)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,43));
xlabel('q knee');
ylabel('q ankle');
title('GL');
legend('Polynomial','Model');
title('lMT right');

%% Assert results
for i = 1:NMuscle  
    assertLMT(:,i) = abs(lMT_out_r(:,i) - MuscleData.lMT(:,i));
    assertdM.hip.flex(:,i) = abs(dM_out_r(:,i,1) - MuscleData.dM(:,i,1));
    assertdM.hip.add(:,i) = abs(dM_out_r(:,i,2) - MuscleData.dM(:,i,2));
    assertdM.hip.rot(:,i) = abs(dM_out_r(:,i,3) - MuscleData.dM(:,i,3));
    assertdM.knee(:,i) = abs(dM_out_r(:,i,4) - MuscleData.dM(:,i,4));
    assertdM.ankle(:,i) = abs(dM_out_r(:,i,5) - MuscleData.dM(:,i,5));
    assertdM.sub(:,i) = abs(dM_out_r(:,i,6) - MuscleData.dM(:,i,6));
    assertdM.lumb.ext(:,i) = abs(dM_out_r(:,i,7) - MuscleData.dM(:,i,7));
    assertdM.lumb.bend(:,i) = abs(dM_out_r(:,i,8) - MuscleData.dM(:,i,8));
    assertdM.lumb.rot(:,i) = abs(dM_out_r(:,i,9) - MuscleData.dM(:,i,9));
end

assertLMTmax_r = max(max(assertLMT));
assertdM.hip.flexmax = max(max(assertdM.hip.flex));
assertdM.hip.addmax = max(max(assertdM.hip.add));
assertdM.hip.rotmax = max(max(assertdM.hip.rot));
assertdM.kneemax = max(max(assertdM.knee));
assertdM.anklemax = max(max(assertdM.ankle));
assertdM.lumb.extmax = max(max(assertdM.lumb.ext));
assertdM.lumb.bendmax = max(max(assertdM.lumb.bend));
assertdM.lumb.rotmax = max(max(assertdM.lumb.rot));

 