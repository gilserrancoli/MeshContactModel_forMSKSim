%% Muscle-driven simulations including mesh-based contact
clear all;
close all;

import casadi.*


solveProblem = true; %Set true to solve the optimal control problem.
analyseResults = true; %Set true to analyze the results.
loadResults = false; %Set true to load the results of the optimization.
saveResults = true; %Set true to save the results of the optimization.
checkBoundsIG = false; % visualize guess-bounds
writeMotionFiles = true; %Set true to write motion files for use in OpenSim GUI
saveOptimalTrajectories = true; %Set true to save optimal trajectories
writeIKmotion=true; %Set true to write .mot file

subject = 'subject_GC';
parallelMode = 'thread';
NThreads = 8; % Number of threads used in parallel.

pathmain = pwd;
[pathRepo,~,~] = fileparts(pathmain);
pathSettings = [pathRepo,'/Settings'];
addpath(genpath(pathSettings));

 % Load settings
 %Use Tangential force with static and dynamic coefficients of 0.8 and 0.8,
    %or 0.1 and 0.1
    Options.useContactModelwithTangentialForce0808=0;
    Options.useContactModelwithTangentialForce0101=0;
 %Related to joint force and moment residuals
    Options.useKCFresiduals=0;
        Options.useKCFconstraintsasBounds=0;
    Options.useJointResMom=0; % originally 1
        Options.useJRMconstraitsasBounds=1;
    Options.usePelvisResMom=0;
    Options.minBothContactForces=0;
    Options.mindiffKCFGRF=0;
    Options.dampingInKneeSec=1; %use damping at the knee secondary dofs
    Options.defineGRFphases=0;
%Related to polynomials approximating muscle-tendon lengths
    Options.useReducedPolynomials=1;
        Options.useReoptimizedPoly=1;
    Options.usedCompiled_lMT_vMT_External=0;
    Options.err_poly = 1; % Define option before running the code
%Related to mesh-based contact model
    %Choose .dll
    Options.maxsmoothness='MellowMax'; %options: logSum, MellowMax, nosmooth
        Options.kInmaxpen=1e4; %with logsum valid numbers 1e3 for TrackSim_2_kneeCont_kmax1e3.dll, 1e5 TrackSim_2_kneeCont_logSum_kmax1e5_kpress5e4.dll, 5e5 TrackSim_2_kneeCont_kmax1e5.dll, 1e4 for TrackSim_2_kneeCont_logSum_kmax1e4_kpress5e4.dll,5e4 for TrackSim_2_kneeCont_logSum_kmax5e4_kpress5e4.dll 
        % 'nosmooth' for non-smooth version TrackSim_2_kneeCont_maxnosmooth.dll
        Options.kInpress=1e4; %with mellowmax valid for 5e4 and 5e5
    Options.kInCheckContacts=1e3; %valid values 1e3, 5e3 and 1e4
    Options.nfacesTib=100;
    Options.nfacesFem=342;
    Options.rad4Pairs=1;
    Options.numpairs=0; %1378, 499... old way to search for selected .dll
    Options.minPressuredev=0;
    Options.minPressures=1;
%Related to IG
    IGi=2; % Data-informed initial guess
    Options.multiple_IG=0   ;
    Options.IGn=1;
    Options.force2createNewIGs=0;
    Options.noiseApplied2NewIGs=1; %percentage
    Options.preoptimizekneedofs=3; %0 is considering all dofs as 
    % experimental data or 0 values, 1 is preoptimizing knee ty, 2 is 
    % preoptimizing knee adduction and knee ty to match knee contact
    % forces, 3 is preoptimizing knee adduction and ty to start with
    % 200 N of knee contact force at each compartment
%Related to how external function is handled
        Options.KCFasinputstoExternalFunction=1;   
hessi=1; %approximate Hession hessi = 2 is the exact hessian
deri=1; %use of AD Recorder
Options.nametrial='ngait_og1';
nametrial.id=Options.nametrial;

% Fixed parameters
W.a         = 4;    % weight muscle activations    was 10
W.u         = 0.01; 
W.Ak        = 50; %weight joint accelerations 
W.AkArms    = 1; %weight joint accelerations for arms was 1
W.ArmE      = 100; %weight arm excitations %initially 1e6
W.passMom   = 0.01; %weight passive torques %original 1000?
W.passFor   = 1e-3;
W.u_vA      = 0.0052;
W.E         = 0.08; % weight metabolic energy rate %original 500 was 0.1
W.Qskneesec = 10; %weight minimize knee secondary dofs (knee add and rot, tx and tz)
W.minBothKCF = 0; %weight to minimize the squared sum of both KCFs
W.diffKCFGRF= 0;
W.minpressuredev=2.5; % was 0.005 tried with 2.5
W.minpressures=3; %initial 3e-4,

exp_E       = 2; % power metabolic energy
IGsel       = 2; % initial guess identifier
cm          = 1; % contact model identifier
IGm         = 1; % initial guess mode identifier
IGcase      = 0; % initial guess case identifier
mE          = 0; % metabolic energy model identifier
coCont      = 0; % co-contraction identifier
v_tgt       = 1.33; %average speed

if Options.useKCFresiduals
    W.KCF_res=1e-1;
else
    W.KCF_res=0;
end
if Options.useJointResMom
    W.JRM_res=1;
else
    W.JRM_res=0;
end
tol_ipopt   = 4;    % tolerance (means 1e-(tol_ipopt))
N           = 50;   % number of mesh intervals

if Options.nfacesTib==49
    Options.nfacestib1=26;
    Options.nfacestib2=23;
elseif Options.nfacesTib==100
    Options.nfacestib1=51;
    Options.nfacestib2=49;
end

setup.derivatives = 'AD_Recorder'; % Algorithmic differentiation / Recorder   
 % Identifiers for experimental data (only used for initial guess or bound)
nametrial.ID    = ['ID_',nametrial.id];
nametrial.GRF   = ['GRF_',nametrial.id];
nametrial.IK    = ['IK_',nametrial.id];
nametrial.KCF   = ['KCF_',nametrial.id];
switch nametrial.id
    case 'ngait_og1'
        time_opt = [2.2292 3.3740]; % right leg cycle [1.680833, 2.775] --> need left GRF at the beginning
    case 'ngait_og1_old'
        time_opt = [2.2292 3.3740]; % right leg cycle [1.680833, 2.775] --> need left GRF at the beginning
    case 'ngait_og5'
        time_opt = [2.7158 3.919]; 
    case 'bouncy4'
        time_opt = [3.5433 4.825];
    case 'bouncy7'
        time_opt = [0.7175 1.965];
    case 'mtpgait3'
        time_opt = [2.1508 3.324];
    case 'mtpgait9'
        time_opt = [2.2583 3.482];
    case 'ngait_tm_fast1'
        time_opt = [4.7067 5.7983]; %one cycle apparently neat
    case 'ngait_tm_set1'
        time_opt = [2.2000 3.5033];
end  

%set save name 
if Options.useReducedPolynomials
    if Options.useReoptimizedPoly
        poly=['red' num2str(Options.err_poly) '_reopt'];
    else
        poly=['red' num2str(Options.err_poly)];
    end
else
    poly=['full' num2str(Options.err_poly)];
end
% The filename used to save the results depends on the settings 
savename = ['N',num2str(N),'_',num2str(Options.IGn),'_poly',poly];
savename2 = [savename,'_with_',num2str(Options.noiseApplied2NewIGs),'_percent_noise']; %used to save the gait trial name plus its noise when multiple IGs are applied

%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. 
% We use different external functions. A first external function extracts 
% several parameters of the bodies to which the contact spheres are attached.
% This is just to access some paramaters of the model in a post-processing
% phase.
pathExternalFunctions = [pathRepo,'/ExternalFunctions_GC'];
% Loading external functions. 
cd(pathExternalFunctions);
switch setup.derivatives
    case 'AD_Recorder'   
        if hessi == 1
            F1 = external('F','TrackSim_1_kneeCont.dll');
            if strcmp(Options.maxsmoothness,'MellowMax')
                if Options.kInmaxpen==1e4
                    if (Options.nfacesTib==49)&&(Options.nfacesFem==171)&&(Options.rad4Pairs==0.5)
                        if (Options.kInpress==5e5)&&(Options.kInCheckContacts==1e3)
                            F2 = external('F','PredSim_2_kneeCont_kmaxpen1e4_kpress5e5_checkContact1e3_49x171_rad4pairs05.dll');
                            F2_debug=external('F','PredSim_2_kneeCont_kmaxpen1e4_kpress5e5_checkContact1e3_49x171_rad4pairs05_debug.dll');
                        
                        end
                    elseif (Options.nfacesTib==100)&&(Options.nfacesFem==342)&&(Options.rad4Pairs==1)
                            if (Options.kInpress==1e4)&&(Options.kInCheckContacts==1e3)
                                if Options.KCFasinputstoExternalFunction
                                    F2_skeletal=external('F','PredSim_2_kneeCont_KCFasinput.dll');
                                    F2_skeletal_debug=external('F','PredSim_2_kneeCont_KCFasinput_forDebug.dll');
                                    F2 = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress1e4_checkContact1e3_100x342_rad1.dll');
                                    F2_debug = external('F','TrackSim_2_kneeCont_MellowMax_kmax1e4_kpress1e4_checkContact1e3_100x342_rad1_debug.dll');
                                end
                            end
                    end
                end
            end
        elseif hessi == 2
            disp('Memory issue with exact Hessian; case not available')
        end
    case 'AD_ADOLC' 
        disp('ADOL-C cases not available');
    case 'FD'
        %deprecated
        F1 = external('F','TrackSim_1_kneeCont.dll',struct('enable_fd',true,...
            'enable_forward',false,'enable_reverse',false,...
            'enable_jacobian',false,'fd_method','forward'));
        F2 = external('F','TrackSim_2_kneeCont.dll',struct('enable_fd',true,...
            'enable_forward',false,'enable_reverse',false,...
            'enable_jacobian',false,'fd_method','forward'));        
end
cd(pathmain);

%% Indices external function
% External function: F2
% First, joint torques. 
jointi.pelvis.tilt  = 1; 
jointi.pelvis.list  = 2; 
jointi.pelvis.rot   = 3; 
jointi.pelvis.tx    = 4;
jointi.pelvis.ty    = 5;
jointi.pelvis.tz    = 6;
jointi.hip_flex.l   = 7;
jointi.hip_add.l    = 8;
jointi.hip_rot.l    = 9;
jointi.hip_flex.r   = 10;
jointi.hip_add.r    = 11;
jointi.hip_rot.r    = 12;
jointi.knee_flex.l  = 13;
jointi.knee_flex.r  = 14;
jointi.knee_add.r   = 15;
jointi.knee_rot.r   = 16;
jointi.knee_tx.r    = 17;
jointi.knee_ty.r    = 18;
jointi.knee_tz.r    = 19;
jointi.ankle.l      = 20;
jointi.ankle.r      = 21;
jointi.subt.l       = 22;
jointi.subt.r       = 23;
jointi.trunk.ext    = 24;
jointi.trunk.ben    = 25;
jointi.trunk.rot    = 26;
jointi.sh_flex.l    = 27;
jointi.sh_add.l     = 28;
jointi.sh_rot.l     = 29;
jointi.sh_flex.r    = 30;
jointi.sh_add.r     = 31;
jointi.sh_rot.r     = 32;
jointi.elb.l        = 33;
jointi.elb.r        = 34;
% jointi.prosup.l     = 35; %fixed values
% jointi.prosup.r     = 36;
% Calcaneus
calcOr.r    = 35:36;
calcOr.l    = 37:38;
calcOr.all  = [calcOr.r,calcOr.l];
NcalcOr     = length(calcOr.all);
% Femurs
femurOr.r   = 39:40;
femurOr.l   = 41:42;
femurOr.all = [femurOr.r,femurOr.l];
NfemurOr    = length(femurOr.all);
% Hands
handOr.r    = 43:44;
handOr.l    = 45:46;
handOr.all  = [handOr.r,handOr.l];
NhandOr     = length(handOr.all);
% Tibias
tibiaOr.r   = 47:48;
tibiaOr.l   = 49:50;
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
NtibiaOr    = length(tibiaOr.all);


% Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all 
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
residuals_noarmsi  = jointi.pelvis.tilt:jointi.trunk.rot; % all but arms
residuals_acti      = [jointi.hip_flex.l:jointi.knee_flex.r jointi.ankle.r:jointi.elb.r]; % all but gr-pelvis and knee
residual_bptyi      = [jointi.pelvis.tilt:jointi.pelvis.tx,...
    jointi.pelvis.tz:jointi.elb.r]; % all but pelvis_ty
knee_secDOF         = jointi.knee_add.r:jointi.knee_tz.r;
% Number of degrees of freedom for later use
nq.all              = length(residualsi); % all 
nq.abs              = length(ground_pelvisi); % ground-pelvis
nq.act              = nq.all-nq.abs;% all but ground-pelvis
nq.trunk            = length(trunki); % trunk
nq.arms             = length(armsi); % arms
nq.leg              = 9; % #joints needed for polynomials
nq.legright         = 14;
nq.legleft          = 9;
Qsi                 = 1:2:2*nq.all; % indices Qs only
% % Second, GRFs
% GRFi.r              = 35:37;
% GRFi.l              = 38:40;
% GRFi.all            = [GRFi.r,GRFi.l];
% nGRF                = length(GRFi.all);
% % Third, GRMs
% GRMi.r              = 41:43;
% GRMi.l              = 44:46;
% GRMi.all            = [GRMi.r,GRMi.l];
% nGRM                = length(GRMi.all);
KCFi.L              = 35;
KCFi.M              = 36;

%% Model info
body_mass = 62.6;
body_weight = body_mass*9.81;

%% Collocation scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points. 
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Muscle-tendon parameters 
% Muscles from one leg and from the back
muscleNames = {'addbrev','addlong','addmagProx','addmagMid',...
        'addmagDist','addmagIsch','bflh','bfsh','edl','ehl','fdl','fhl',...
        'gaslat','gasmed','gem','glmax1','glmax2','glmax3','glmed1',...
        'glmed2','glmed3','glmin1','glmin2','glmin3','grac','iliacus',...
        'pect','perbrev','perlong','pertert','piri','psoas','quadfem',...
        'recfem','sart','semimem','semiten','soleus','tfl',...
        'tibant','tibpost','vasint','vaslat',...
        'vasmed','ercspn_r','intobl_r','extobl_r',...
        'ercspn_l','intobl_l','extobl_l'};
% Muscle indices for later use
pathmusclemodel = [pathRepo,'/MuscleModel'];
addpath(genpath(pathmusclemodel));    
% (1:end-3), since we do not want to count twice the back muscles
musi = MuscleIndices_3D_GC(muscleNames(1:end-3));
% Total number of muscles
NMuscle = length(musi)*2;
% Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation 
% angles; Row 5: maximal contraction velocities
load([pathmusclemodel,'/MTparameters_',subject,'.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];
%Slightly modify lM0 to avoid huge passive forces at the initial guess
MTparameters_m(2,56)=0.09; %edl l lM0;
MTparameters_m(2,59)=0.08; %fhl l lM0;
MTparameters_m(2,60)=0.1; %gaslat l lM0
MTparameters_m(2,61)=0.1; %gasmed l lM0
MTparameters_m(2,62)=0.04; %gem l lM0
MTparameters_m(2,75)=0.07; %perbrev l lM0
MTparameters_m(2,76)=0.1; %perlong l lM0
MTparameters_m(2,77)=0.09; %pertert l lM0
MTparameters_m(2,78)=0.04; %piri l lM0
MTparameters_m(2,81)=0.1; %recfem l lM0
MTparameters_m(2,85)=0.1; %soleus l lM0
MTparameters_m(2,87)=0.1; %tibant l lM0
MTparameters_m(2,88)=0.06; %tibpost l lM0
MTparameters_m(2,90)=0.105; %vaslat l lM0

MTparameters_m(2,9)=0.09; %edl r lM0;
MTparameters_m(2,12)=0.08; %fhl r lM0;
MTparameters_m(2,13)=0.1; %gaslat r lM0
MTparameters_m(2,14)=0.1; %gasmed r lM0
MTparameters_m(2,15)=0.04; %gem r lM0
MTparameters_m(2,28)=0.07; %perbrev r lM0
MTparameters_m(2,29)=0.1; %perlong r lM0
MTparameters_m(2,30)=0.09; %pertert r lM0
MTparameters_m(2,31)=0.04; %piri r lM0
MTparameters_m(2,34)=0.1; %recfem r lM0
MTparameters_m(2,38)=0.1; %soleus r lM0
MTparameters_m(2,40)=0.1; %tibant r lM0
MTparameters_m(2,41)=0.06; %tibpost r lM0
MTparameters_m(2,43)=0.105; %vaslat r lM0

%These are just to check if semimem can operate closer to lMtilde=1
MTparameters_m(2,36)=0.06; %semimem_l
MTparameters_m(2,36+47)=0.06; %semimem_r
MTparameters_m(3,36)=0.345; %semimem_l
MTparameters_m(3,36+47)=0.345; %semimem_r
MTparameters_m(3,38)=0.29; %soleus_l
MTparameters_m(3,38+47)=0.29; %soleus_r

MTparameters_m(1,:)=MTparameters_m(1,:)*1.5;

% Indices of the muscles actuating the different joints for later use
pathpolynomial = [pathRepo,'/Polynomials_GC'];
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
[~,mai] = MomentArmIndices_3D(muscleNames(1:end-3),...
    tl.muscle_spanning_joint_INFO(1:end-3,:));

% By default, the tendon stiffness is 35 and the shift is 0.
aTendon = 35*ones(NMuscle,1);
shift = zeros(NMuscle,1);

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Contact model parameters
% The locations of the contact spheres (x and z coordinates), the radii and 
% the other parameters of the foot-ground contact model are 
% needed to compute GRFs. 
% Indices variable parameters
loci.s1.r.x = 1;
loci.s1.r.z = 2;
loci.s2.r.x = 3;
loci.s2.r.z = 4;
loci.s3.r.x = 5;
loci.s3.r.z = 6;
loci.s4.r.x = 7;
loci.s4.r.z = 8;
loci.s5.r.x = 9;
loci.s5.r.z = 10;
loci.s6.r.x = 11;
loci.s6.r.z = 12;
radi.s1 = 13;
radi.s2 = 14;
radi.s3 = 15;
radi.s4 = 16;
radi.s5 = 17;
radi.s6 = 18;
% Fixed parameters
dissipation = 2; %check this value
normal = [0,1,0];
transitionVelo = 0.2;
staticFriction = 0.8;
dynamicFriction = 0.8;
viscousFriction = 0.5; 
stif = 1000000;  
loc.s1.y    = -0.021859;
loc.s2.y    = -0.021859;
loc.s3.y    = -0.021859;
loc.s4.y    = -0.0214476;
loc.s5.y    = -0.021859;
loc.s6.y    = -0.0214476;

%% Metabolic energy model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
% (1:end-3), since we do not want to count twice the back muscles
tension = getSpecificTensions(muscleNames(1:end-3)); 
tensions = [tension;tension];
% (1:end-3), since we do not want to count twice the back muscles
pctst = getSlowTwitchRatios(muscleNames(1:end-3)); 
pctsts = [pctst;pctst];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% acoording to the order of 94 muscles
 
%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));
pathContactModel = [pathRepo,'/Contact'];
addpath(genpath(pathContactModel));
% We load some variables for the polynomial approximations
load([pathpolynomial,'/muscle_spanning_joint_INFO_',subject,'.mat']);
if Options.useReducedPolynomials
    if Options.useReoptimizedPoly
        switch Options.err_poly 
            case 0.5
                load([pathpolynomial,'/MuscleInfo_red_05_ropt','.mat']);
                %MuscleInfo = MuscleInfo_red_05;
            case 1
                load([pathpolynomial,'/MuscleInfo_red_1_ropt','.mat']);
                %MuscleInfo = MuscleInfo_red_1;
            case 2
                load([pathpolynomial,'/MuscleInfo_red_2_ropt','.mat']);
                %MuscleInfo = MuscleInfo_red_2;
            case 3
                load([pathpolynomial,'/MuscleInfo_red_3_ropt','.mat']);
                %MuscleInfo = MuscleInfo_red_3;
        end
    else
        switch Options.err_poly 
            case 0.5
                load([pathpolynomial,'/MuscleInfo_red_05','.mat']);
                %MuscleInfo = MuscleInfo_red_05;
            case 1
                load([pathpolynomial,'/MuscleInfo_red_1','.mat']);
                %MuscleInfo = MuscleInfo_red_1;
            case 2
                load([pathpolynomial,'/MuscleInfo_red_2','.mat']);
                %MuscleInfo = MuscleInfo_red_2;
            case 3
                load([pathpolynomial,'/MuscleInfo_red_3','.mat']);
                %MuscleInfo = MuscleInfo_red_3;
        end
    end
else
    switch Options.err_poly 
        case 0.5
            load([pathpolynomial,'/MuscleInfo_full_05','.mat']);
            %MuscleInfo = MuscleInfo_full_05;
        case 1
            load([pathpolynomial,'/MuscleInfo_full_1','.mat']);
            %MuscleInfo = MuscleInfo_full_1;
        case 2
            load([pathpolynomial,'/MuscleInfo_full_2','.mat']);
            %MuscleInfo = MuscleInfo_full_2;
        case 3
            load([pathpolynomial,'/MuscleInfo_full_3','.mat']);
            %MuscleInfo = MuscleInfo_full_3;

    end
end

% For the polynomials, we want all independent muscles. So we do not need
% the muscles from both legs, since we assume bilateral symmetry, but want
% all muscles from the back (indices 48:50).
musi_pol = [musi,48,49,50];
NMuscle_pol = NMuscle/2+3;
CasADiFunctions_3D_GC_pred
cd(pathExternalFunctions);
if Options.usedCompiled_lMT_vMT_External
    if Options.useReducedPolynomials
        keyboard;
        %deprecated
        f_lMT_vMT_dM = external('f_lMT_vMT_dM','f_lMT_vMT_reduced.dll');
    else
        f_lMT_vMT_dM = external('f_lMT_vMT_dM','f_lMT_vMT_full.dll');
    end
end
cd(pathmain);

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk
pathPassiveMoments = [pathRepo,'/PassiveMoments'];
addpath(genpath(pathPassiveMoments));
PassiveMomentsData

%% Experimental data
% We extract experimental data to set bounds and initial guesses if needed
% (not within the cost function or constraints)
pathData = [pathRepo,'/OpenSimModel/',subject];
joints = {'pelvis_rz','pelvis_rx','pelvis_ry','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion','hip_adduction','hip_rotation',...
    'knee_flexion_l','knee_flexion','knee_adduction','knee_rotation',...
    'knee_tx','knee_ty','knee_tz','ankle_angle_l','ankle_angle',...
    'subtalar_angle_l','subtalar_angle',...
    'lumbar_extension','lumbar_bending','lumbar_rotation',...
    'arm_flex_l','arm_add_l','arm_rot_l',...
    'arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));
% Extract joint kinematics
pathIK = [pathData,'/IK/',nametrial.IK,'.mat'];
Qs = getIK_GC(pathIK,joints);
% Extract ground reaction forces and moments
pathGRF = [pathData,'/GRF/',nametrial.GRF,'.mat'];
GRF = getGRF_GC(pathGRF);
% Extract joint kinetics
pathID = [pathData,'/ID/',nametrial.ID,'.mat'];
ID = getID_GC(pathID,joints);
%Extract knee contact forces
pathKCF=[pathData, '/KCF/',nametrial.KCF,'.mat'];
KCF=get_KCF(pathKCF);
% Interpolation experimental data
sampfreq_kin=1/120;
time_expi.ID(1) = find((ID.time<(time_opt(1)+sampfreq_kin/2))&(ID.time>=(time_opt(1)-sampfreq_kin/2)));
time_expi.ID(2) = find((ID.time<(time_opt(2)+sampfreq_kin/2))&(ID.time>=(time_opt(2)-sampfreq_kin/2)));
time_expi.KCF(1) = find((KCF.time<(time_opt(1)+sampfreq_kin/2))&(KCF.time>=(time_opt(1)-sampfreq_kin/2)));
time_expi.KCF(2) = find((KCF.time<(time_opt(2)+sampfreq_kin/2))&(KCF.time>=(time_opt(2)-sampfreq_kin/2)));
sampfreq_GRF=1/1200;
time_expi.GRF(1) = find((GRF.time<(time_opt(1)+sampfreq_GRF/2))&(GRF.time>=(time_opt(1)-sampfreq_GRF/2)));
time_expi.GRF(2) = find((GRF.time<(time_opt(2)+sampfreq_GRF/2))&(GRF.time>=(time_opt(2)-sampfreq_GRF/2)));
% step = (ID.time(time_expi.ID(2))-ID.time(time_expi.ID(1)))/N;
step=(time_opt(2)-time_opt(1))/N;
interval = time_opt(1):step:time_opt(2);
tgrid_col=zeros(N*d,1);
tgrid_col(1:d:end)=interval(1:end-1)+tau_root(2)*step;
tgrid_col(2:d:end)=interval(1:end-1)+tau_root(3)*step;
tgrid_col(3:d:end)=interval(1:end-1)+tau_root(4)*step;
ID.allinterp = interp1(ID.all(:,1),ID.all,interval);
ID.allinterp_col = interp1(ID.all(:,1),ID.all,tgrid_col);
Qs.allinterpfilt = interp1(Qs.allfilt(:,1),Qs.allfilt,interval);
Qs.allinterpfilt_col = interp1(Qs.allfilt(:,1),Qs.allfilt,tgrid_col);
KCF.allinterpfilt = interp1(KCF.data(:,1),KCF.data(:,2:3),interval);
KCF.allinterpfilt_col=interp1(interval,KCF.allinterpfilt,tgrid_col);
GRF.val.allinterp = interp1(round(GRF.val.all(:,1),4),...
    GRF.val.all,round(interval,4));
GRF.val.allinterp_col=interp1(interval,GRF.val.allinterp,tgrid_col);
GRF.MorGF.allinterp = interp1(round(GRF.MorGF.all(:,1),4),...
    GRF.MorGF.all,round(interval,4));
GRF.MorGF.allinterp_col=interp1(interval,GRF.MorGF.allinterp,tgrid_col);
GRF.pos.allinterp = interp1(round(GRF.pos.all(:,1),4),...
    GRF.pos.all,round(interval,4));
GRF.pos.allinterp_col=interp1(interval,GRF.pos.allinterp,tgrid_col);
GRF.Mcop.allinterp = interp1(round(GRF.Mcop.all(:,1),4),...
    GRF.Mcop.all,round(interval,4));
GRF.Mcop.allinterp_col=interp1(interval,GRF.Mcop.allinterp,tgrid_col);
    
%% Bounds
pathBounds = [pathRepo,'/Bounds'];
addpath(genpath(pathBounds));
dev_cm.loc=0; %not used here, just to use the same function to get bounds
dev_cm.rad=0; %not used here, just to use the same function to get bounds
[bounds,scaling] = getBounds_3D_GC(Qs,NMuscle,nq,jointi,dev_cm,GRF,'pred');
%% Bounds for Final time
bounds.tf.lower = 0.5;
bounds.tf.upper = 1.5;
%% Hard bounds
% We impose the initial position of pelvis_tx to be 0
bounds.QsQdots_0.lower = bounds.QsQdots.lower;
bounds.QsQdots_0.upper = bounds.QsQdots.upper;
bounds.QsQdots_0.lower(2*jointi.pelvis.tx-1) = 0;
bounds.QsQdots_0.upper(2*jointi.pelvis.tx-1) = 0;

%% Scale Qs experimental 
Qs_res_interpfilt_scaled(:,jointi.knee_ty.r)= (Qs.allinterpfilt(:,jointi.knee_ty.r+1)-scaling.knee_ty.a)/scaling.knee_ty.b;
Qs_res_interpfilt_scaled(:,[1:jointi.knee_tx.r jointi.knee_tz.r:nq.all])=    Qs.allinterpfilt(:, [1:jointi.knee_tx.r jointi.knee_tz.r:nq.all]+1)./scaling.Qs([1:jointi.knee_tx.r jointi.knee_tz.r:end]);          
Qs_res_interpfilt_scaled_aux=Qs_res_interpfilt_scaled(:,[1:4 6:14 20:34]); %exclude knee internal dofs
Qs_res_interpfilt_scaled_aux_col=interp1(interval(1:N+1),Qs_res_interpfilt_scaled_aux,tgrid_col);

for i=2:size(Qs.allinterpfilt,2)
    Q_spline(i-1)=spline(Qs.allinterpfilt(:,1),Qs.allinterpfilt(:,i));
    Qdot_spline(i-1)=fnder(Q_spline(i-1),1);
    Qdots.allinterpfilt(:,1)=Qs.allinterpfilt(:,1);
    Qdots.allinterpfilt(:,i)=ppval(Qdot_spline(i-1),Qdots.allinterpfilt(:,1));
end
Qdots_res_interpfilt_scaled=Qdots.allinterpfilt(:,2:end)./scaling.Qdots;
Qdots_res_interpfilt_scaled_aux=Qdots_res_interpfilt_scaled(:,[1:4 6:14 20:34]); %exclude knee internal dofs
Qdots_res_interpfilt_scaled_aux_col=interp1(interval(1:N+1),Qdots_res_interpfilt_scaled_aux,tgrid_col);

%% Initial guess
pathIG = [pathRepo,'/Guess'];
addpath(genpath(pathIG));
if (Options.multiple_IG==0)
else
    if ~exist([pathIG '\multipleIGs_v2', savename2, '.mat'], 'file')&(Options.multiple_IG==1)
        Options.force2createNewIGs=1; 
    else
        load([pathIG '\multipleIGs_v2', savename2, '.mat']);
        if ~isfield(IGs,nametrial.id)
            Options.force2createNewIGs=1; 
        end
    end
end
if (Options.multiple_IG==0)|(Options.force2createNewIGs==1)
    % The initial guess depends on the settings
    if IGi == 1 % Quasi-random initial guess  
        guess = getGuess_3D_QR(Qs,nq,N,NMuscle,jointi,time_opt,scaling);
    elseif IGi == 2 % Data-informed initial guess  
        if Options.KCFasinputstoExternalFunction
            guess = getGuess_3D_DI_GC_v2(Qs,nq,N,NMuscle,jointi,scaling,d,KCF,GRF,F1,{F2_skeletal_debug F2},f_contactForce,Options,'pred');
        else
           guess = getGuess_3D_DI_GC_v2(Qs,nq,N,NMuscle,jointi,scaling,d,KCF,GRF,F1,F2,f_contactForce,Options,'pred');
        end
    elseif IGi == 3 % Data-informed initial guess with muscle information obtained 
        % using the muscle redundancy solver 
        % (https://simtk.org/projects/optcntrlmuscle)
        load([pathIG,'/mvar.mat'],'mvar');
        mvarinterp.a = interp1(mvar.t(:,1),mvar.a,interval);
        mvarinterp.FTtilde = interp1(mvar.t(:,1),mvar.FTtilde,interval);
        mvarinterp.vA = interp1(mvar.t(1:end-1,1),mvar.vA,interval);
        mvarinterp.dFTtilde = interp1(mvar.t(1:end-1,1),mvar.dFTtilde,interval);
        guess = getGuess_3D_DIm(Qs,mvarinterp,nq,N,jointi,scaling);
    end
end
if size(guess.a_a,1)==N
    guess.a_a(N+1,:)=guess.a_a(N,:);
end

if Options.force2createNewIGs==1
    guess_i(1)=guess;
    IGnoise=0.01*Options.noiseApplied2NewIGs;
    for i=1:9
        guess_i(i+1).Qs=guess.Qs+IGnoise*randn([N+1 nq.all]);
        guess_i(i+1).Qdots=guess.Qdots+IGnoise*randn([N+1 nq.all]);
        guess_i(i+1).QsQdots=guess.QsQdots+IGnoise*randn([N+1 nq.all*2]);
        guess_i(i+1).Qdotdots=guess.Qdotdots+IGnoise*randn([N+1 nq.all]);
        guess_i(i+1).a=guess.a+IGnoise*randn([N+1 NMuscle]);
        guess_i(i+1).a(guess_i(i+1).a<1e-4)=1e-4;
        guess_i(i+1).vA=guess.vA+IGnoise*randn([N NMuscle]);
        guess_i(i+1).FTtilde=guess.FTtilde+IGnoise*randn([N+1 NMuscle]);
        guess_i(i+1).FTtilde(guess_i(i+1).FTtilde<1e-4)=1e-4;
        guess_i(i+1).dFTtilde=guess.dFTtilde+IGnoise*randn([N NMuscle]);
        guess_i(i+1).a_a=guess.a_a+IGnoise*randn(size(guess.a_a));
        guess_i(i+1).e_a=guess.e_a+IGnoise*randn(size(guess.e_a));
        guess_i(i+1).params=guess.params+IGnoise*randn(size(guess.params));
        guess_i(i+1).KCFresiduals=guess.KCFresiduals+IGnoise*randn(size(guess.KCFresiduals));
        guess_i(i+1).JRM_res=guess.JRM_res+IGnoise*randn(size(guess.JRM_res));
        guess_i(i+1).a_col=guess.a_col+IGnoise*randn(size(guess.a_col));
        guess_i(i+1).a_col(guess_i(i+1).a_col<1e-4)=1e-4;
        guess_i(i+1).FTtilde_col=guess.FTtilde_col+IGnoise*randn(size(guess.FTtilde_col));
        guess_i(i+1).QsQdots_col=guess.QsQdots_col+IGnoise*randn(size(guess.QsQdots_col));
        guess_i(i+1).a_a_col=guess.a_a_col+IGnoise*randn(size(guess.a_a_col));
    end
    IGs.(nametrial.id).guess_i=guess_i;
    save([pathIG '\multipleIGs_v2',savename2, '.mat'],'IGs');
end
if Options.multiple_IG==1
    guess=IGs.(nametrial.id).guess_i(Options.IGn);
end

%% Initial guess for Final time
% The final time is function of the imposed speed
all_speeds = 0.73:0.1:2.73;
all_tf = 2*(0.70:-((0.70-0.35)/(length(all_speeds)-1)):0.35);
idx_speed = find(all_speeds==v_tgt);
if isempty(idx_speed)
    idx_speed = find(all_speeds > v_tgt,1,'first');
end
guess.tf = all_tf(idx_speed);

%% Parameters contact model
% Original values
B_locSphere_s1_r    = [0.00190115788407966, -0.00382630379623308];
B_locSphere_s2_r    = [0.148386399942063, -0.028713422052654];
B_locSphere_s3_r    = [0.133001170607051, 0.0516362473449566];
B_locSphere_s4_r    = [0.06, -0.0187603084619177];    
B_locSphere_s5_r    = [0.0662346661991635, 0.0263641606741698];
B_locSphere_s6_r    = [0.045, 0.0618569567549652];
IG_rad              = 0.032*ones(1,6); 
paramsCM_nsc = [B_locSphere_s1_r,B_locSphere_s2_r,...
    B_locSphere_s3_r,B_locSphere_s4_r,B_locSphere_s5_r,B_locSphere_s6_r,...
    IG_rad];

%% This allows visualizing the initial guess and the bounds
if checkBoundsIG
    pathPlots = [pathRepo,'/Plots'];
    addpath(genpath(pathPlots));
    plot_BoundsVSInitialGuess_3D_GC
end



%% Formulate the NLP
if solveProblem
    % Create opti instance
    opti = casadi.Opti();
    % Define static parameters
    % Final time
    tf = opti.variable();
    opti.subject_to(bounds.tf.lower < tf < bounds.tf.upper);
    opti.set_initial(tf, guess.tf);
    % Define states 
    % Muscle activations at mesh points
    a = opti.variable(NMuscle,N+1);
    opti.subject_to(bounds.a.lower'*ones(1,N+1) < a < ...
        bounds.a.upper'*ones(1,N+1));
    opti.set_initial(a, guess.a');
    % Muscle activations at collocation points
    a_col = opti.variable(NMuscle,d*N);
    opti.subject_to(bounds.a.lower'*ones(1,d*N) < a_col < ...
        bounds.a.upper'*ones(1,d*N));
    opti.set_initial(a_col, guess.a_col');    
    % Muscle-tendon forces at mesh points
    FTtilde = opti.variable(NMuscle,N+1);
    opti.subject_to(bounds.FTtilde.lower'*ones(1,N+1) < FTtilde < ...
        bounds.FTtilde.upper'*ones(1,N+1));
    opti.set_initial(FTtilde, guess.FTtilde');
    % Muscle-tendon forces at collocation points
    FTtilde_col = opti.variable(NMuscle,d*N);
    opti.subject_to(bounds.FTtilde.lower'*ones(1,d*N) < FTtilde_col < ...
        bounds.FTtilde.upper'*ones(1,d*N));
    opti.set_initial(FTtilde_col, guess.FTtilde_col');    
    % Qs at mesh points
    Qs = opti.variable(nq.all,N+1);
    % We want to constraint the pelvis_tx position at the first mesh point,
    % and avoid redundant bounds
    lboundsQsk = bounds.QsQdots.lower(1:2:end)'*ones(1,N+1);
    lboundsQsk(jointi.pelvis.tx,1) = ...
        bounds.QsQdots_0.lower(2*jointi.pelvis.tx-1);
    uboundsQsk = bounds.QsQdots.upper(1:2:end)'*ones(1,N+1);
    uboundsQsk(jointi.pelvis.tx,1) = ...
        bounds.QsQdots_0.upper(2*jointi.pelvis.tx-1);
    opti.subject_to(lboundsQsk < Qs < uboundsQsk);
    opti.set_initial(Qs, guess.QsQdots(:,1:2:end)');
    % Qs at collocation points
    Qs_col = opti.variable(nq.all,d*N);
    opti.subject_to(bounds.QsQdots.lower(1:2:end)'*ones(1,d*N) < Qs_col < ...
        bounds.QsQdots.upper(1:2:end)'*ones(1,d*N));
    opti.set_initial(Qs_col, guess.QsQdots_col(:,1:2:end)');   
    % Qdots at mesh points
    Qdots = opti.variable(nq.all,N+1);
    opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,N+1) < Qdots < ...
        bounds.QsQdots.upper(2:2:end)'*ones(1,N+1));
    opti.set_initial(Qdots, guess.QsQdots(:,2:2:end)');    
    % Qdots at collocation points
    Qdots_col = opti.variable(nq.all,d*N);
    opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,d*N) < Qdots_col < ...
        bounds.QsQdots.upper(2:2:end)'*ones(1,d*N));
    opti.set_initial(Qdots_col, guess.QsQdots_col(:,2:2:end)');
    % Arm activations at mesh points
    a_a = opti.variable(nq.arms,N+1);
    opti.subject_to(bounds.a_a.lower'*ones(1,N+1) < a_a < ...
        bounds.a_a.upper'*ones(1,N+1));
    opti.set_initial(a_a, guess.a_a');  
    % Arm activations at collocation points
    a_a_col = opti.variable(nq.arms,d*N);
    opti.subject_to(bounds.a_a.lower'*ones(1,d*N) < a_a_col < ...
        bounds.a_a.upper'*ones(1,d*N));
    opti.set_initial(a_a_col, guess.a_a_col');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define controls
    % Time derivative of muscle activations (states) at mesh points
    vA = opti.variable(NMuscle, N);
    opti.subject_to(bounds.vA.lower'*ones(1,N) < vA < ...
        bounds.vA.upper'*ones(1,N));
    opti.set_initial(vA, guess.vA');   
    % Arm excitations
    e_a = opti.variable(nq.arms, N);
    opti.subject_to(bounds.e_a.lower'*ones(1,N) < e_a < ...
        bounds.e_a.upper'*ones(1,N));
    opti.set_initial(e_a, guess.e_a');
    if Options.useKCFresiduals||Options.useJointResMom||Options.usePelvisResMom
        disp('not writting for predictive simulations');
        keyboard;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define "slack" controls
    % Time derivative of muscle-tendon forces (states) at collocation points
    dFTtilde_col = opti.variable(NMuscle, d*N);
    opti.subject_to(bounds.dFTtilde.lower'*ones(1,d*N) < dFTtilde_col < ...
        bounds.dFTtilde.upper'*ones(1,d*N));
    opti.set_initial(dFTtilde_col, guess.dFTtilde_col');
    % Time derivative of Qdots (states) at collocation points
    A_col = opti.variable(nq.all, d*N);
    opti.subject_to(bounds.Qdotdots.lower'*ones(1,d*N) < A_col < ...
        bounds.Qdotdots.upper'*ones(1,d*N));
    opti.set_initial(A_col, guess.Qdotdots_col');          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parallel formulation
    % Define CasADi variables for static parameters
    tfk = MX.sym('tfk');
    % Define CasADi variables for states
    ak = MX.sym('ak',NMuscle);
    aj = MX.sym('akmesh',NMuscle,d);
    akj = [ak aj];
    FTtildek = MX.sym('FTtildek',NMuscle); 
    FTtildej = MX.sym('FTtildej',NMuscle,d);
    FTtildekj = [FTtildek FTtildej];
    Qsk = MX.sym('Qsk',nq.all);
    Qsj = MX.sym('Qsj',nq.all,d);
    Qskj = [Qsk Qsj];
    Qdotsk = MX.sym('Qdotsk',nq.all);
    Qdotsj = MX.sym('Qdotsj',nq.all,d);
    Qdotskj = [Qdotsk Qdotsj];
    a_ak = MX.sym('a_ak',nq.arms);
    a_aj = MX.sym('a_akmesh',nq.arms,d);
    a_akj = [a_ak a_aj];
    % Define CasADi variables for controls
    vAk = MX.sym('vAk',NMuscle);
    e_ak = MX.sym('e_ak',nq.arms);
    % Define CasADi variables for "slack" controls
    dFTtildej = MX.sym('dFTtildej',NMuscle,d);
    Aj = MX.sym('Aj',nq.all,d);  
    J = 0; % Initialize cost function
    JE=0;
    Ja=0;
    Jak=0;
    JpassMom=0;
    JpassFor=0;
    JvAk=0;
    JdFTtilde=0;
    JArmE=0;
    JAjarms=0;
    JQskneesec=0;
    JminbothKCF=0;
    JdiffKCFGRF=0;
    eq_constr = {}; % Initialize equality constraint vector
    ineq_constr1 = {}; % Initialize inequality constraint vector 1
    ineq_constr2 = {}; % Initialize inequality constraint vector 2
    ineq_constr3 = {}; % Initialize inequality constraint vector 3
    ineq_constr4 = {}; % Initialize inequality constraint vector 4
    ineq_constr5 = {}; % Initialize inequality constraint vector 5    
    ineq_constr6 = {}; % Initialize inequality constraint vector 5    
    ineq_constr7 = {}; % Initialize inequality constraint vector 5    
    g_names_coll = {}; % Initialize names of constraints at collocation points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Time step
    h = tfk/N; 
    % Loop over collocation points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unscale variables        
    Qskj_nsc = Qskj.*(scaling.QsQdots(1:2:end)'*ones(1,size(Qskj,2)/2));
    Qskj_nsc(jointi.knee_ty.r,:) = Qskj(jointi.knee_ty.r,:).*scaling.knee_ty.b+scaling.knee_ty.a;
    Qdotskj_nsc = Qdotskj.*(scaling.QsQdots(2:2:end)'*...
        ones(1,size(Qdotskj,2)/2));
    FTtildekj_nsc = FTtildekj.*(scaling.FTtilde'*ones(1,size(FTtildekj,2)));
    dFTtildej_nsc = dFTtildej.*scaling.dFTtilde;
    Aj_nsc = Aj.*(scaling.Qdotdots'*ones(1,size(Aj,2)));  
    vAk_nsc = vAk.*scaling.vA;              
    
    QsQdotskj_nsc = MX(nq.all*2, d+1);
    QsQdotskj_nsc(1:2:end,:) = Qskj_nsc;
    QsQdotskj_nsc(2:2:end,:) = Qdotskj_nsc;
    %%%%%%%%%%%%%%%%% to debug
    lMtildej_all=MX.zeros(94,3);
    lMTj_lr_all=MX.zeros(94,3);
    
    Tj_all=MX.zeros(F2_skeletal.nnz_out,3);
    Tau_passj_all_d=MX.zeros(20,3);
    FTj_all=MX.zeros(94,3);
    MAj_knee_ty_r_all=MX.zeros(13,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:d            
        % Call external functions
        % The second external function (F2) returns joint torques, GRFs,
        % and GRMs based on joint states, controls, contact forces.
        if Options.KCFasinputstoExternalFunction
            q_knee=QsQdotskj_nsc([jointi.knee_flex.r jointi.knee_add.r jointi.knee_rot.r...
                    jointi.knee_tx.r jointi.knee_ty.r jointi.knee_tz.r]*2-1,j+1);
            out_skeletal=F2(q_knee);
            KCF_Mj=out_skeletal(1);
            KCF_Lj=out_skeletal(2);
            SumForces=out_skeletal(3:5);
            SumMoments=out_skeletal(6:8);
            if Options.defineGRFphases
                [Tj] = F2_skeletal([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j);SumForces;SumMoments]);
                [Tj_debug] = F2_skeletal_debug([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j);SumForces;SumMoments]);
                GRFrj=Tj_debug(35:37);
                GRFlj=Tj_debug(38:40);
            else
                [Tj] = F2_skeletal([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j);SumForces;SumMoments]);
            end
        else
            [Tj] = F2([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);
            [Tj_debug] = F2_debug([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);                        
        end
        % Get muscle-tendon lengths, velocities, and moment arms
        % Left leg
        qinj_l = [Qskj_nsc(jointi.hip_flex.l,j+1),...
            Qskj_nsc(jointi.hip_add.l,j+1), ...
            Qskj_nsc(jointi.hip_rot.l,j+1), ...
            0,0,...
            Qskj_nsc(jointi.knee_flex.l,j+1), ...
            0,0.042,0,...
            Qskj_nsc(jointi.ankle.l,j+1),...
            Qskj_nsc(jointi.subt.l,j+1),...
            Qskj_nsc(jointi.trunk.ext,j+1),...
            Qskj_nsc(jointi.trunk.ben,j+1),...
            Qskj_nsc(jointi.trunk.rot,j+1)];  
        qdotinj_l = [Qdotskj_nsc(jointi.hip_flex.l,j+1),...
            Qdotskj_nsc(jointi.hip_add.l,j+1),...
            Qdotskj_nsc(jointi.hip_rot.l,j+1),...
            0,0,...
            Qdotskj_nsc(jointi.knee_flex.l,j+1),...
            0,0,0,...
            Qdotskj_nsc(jointi.ankle.l,j+1),...
            Qdotskj_nsc(jointi.subt.l,j+1),...
            Qdotskj_nsc(jointi.trunk.ext,j+1),...
            Qdotskj_nsc(jointi.trunk.ben,j+1),...
            Qdotskj_nsc(jointi.trunk.rot,j+1)];  
        [lMTj_l,vMTj_l,MAj_l] =  f_lMT_vMT_dM(qinj_l,qdotinj_l);    
        MAj.hip_flex.l   =  MAj_l(mai(1).mus.l',1);
        MAj.hip_add.l    =  MAj_l(mai(2).mus.l',2);
        MAj.hip_rot.l    =  MAj_l(mai(3).mus.l',3);
        MAj.knee_flex.l  =  MAj_l(mai(6).mus.l',6);
        MAj.ankle.l      =  MAj_l(mai(10).mus.l',10);  
        MAj.subt.l       =  MAj_l(mai(11).mus.l',11); 
        % For the back muscles, we want left and right together: left
        % first, right second. In MuscleInfo, we first have the right
        % muscles (44:46) and then the left muscles (47:49). Since the back
        % muscles only depend on back dofs, we do not care if we extract
        % them "from the left or right leg" so here we just picked left.
        MAj.trunk_ext    =  MAj_l([48:50,mai(12).mus.l]',12);
        MAj.trunk_ben    =  MAj_l([48:50,mai(13).mus.l]',13);
        MAj.trunk_rot    =  MAj_l([48:50,mai(14).mus.l]',14);
        % Right leg
        qinj_r = [Qskj_nsc(jointi.hip_flex.r,j+1),...
            Qskj_nsc(jointi.hip_add.r,j+1),...
            Qskj_nsc(jointi.hip_rot.r,j+1),...
            Qskj_nsc(jointi.knee_add.r,j+1),...
            Qskj_nsc(jointi.knee_rot.r,j+1),...
            Qskj_nsc(jointi.knee_flex.r,j+1),...
            Qskj_nsc(jointi.knee_tx.r,j+1),...
            Qskj_nsc(jointi.knee_ty.r,j+1),...
            Qskj_nsc(jointi.knee_tz.r,j+1),...
            Qskj_nsc(jointi.ankle.r,j+1),...
            Qskj_nsc(jointi.subt.r,j+1),...
            Qskj_nsc(jointi.trunk.ext,j+1),...
            Qskj_nsc(jointi.trunk.ben,j+1),...
            Qskj_nsc(jointi.trunk.rot,j+1)];  
        qdotinj_r = [Qdotskj_nsc(jointi.hip_flex.r,j+1),...
            Qdotskj_nsc(jointi.hip_add.r,j+1),...
            Qdotskj_nsc(jointi.hip_rot.r,j+1),...
            Qdotskj_nsc(jointi.knee_add.r,j+1),...
            Qdotskj_nsc(jointi.knee_rot.r,j+1),...
            Qdotskj_nsc(jointi.knee_flex.r,j+1),...
            Qdotskj_nsc(jointi.knee_tx.r,j+1),...
            Qdotskj_nsc(jointi.knee_ty.r,j+1),...
            Qdotskj_nsc(jointi.knee_tz.r,j+1),...
            Qdotskj_nsc(jointi.ankle.r,j+1),...
            Qdotskj_nsc(jointi.subt.r,j+1),...
            Qdotskj_nsc(jointi.trunk.ext,j+1),...
            Qdotskj_nsc(jointi.trunk.ben,j+1),...
            Qdotskj_nsc(jointi.trunk.rot,j+1)];      
        [lMTj_r,vMTj_r,MAj_r] = f_lMT_vMT_dM(qinj_r,qdotinj_r);
        % Here we take the indices from left since the vector is 1:49
        MAj.hip_flex.r   =  MAj_r(mai(1).mus.l',1);
        MAj.hip_add.r    =  MAj_r(mai(2).mus.l',2);
        MAj.hip_rot.r    =  MAj_r(mai(3).mus.l',3);
        MAj.knee_add.r   =  MAj_r(mai(4).mus.l',4);
        MAj.knee_rot.r   =  MAj_r(mai(5).mus.l',5);
        MAj.knee_flex.r  =  MAj_r(mai(6).mus.l',6);
        MAj.knee_tx.r    =  MAj_r(mai(7).mus.l',7);
        MAj.knee_ty.r    =  MAj_r(mai(8).mus.l',8);
        MAj.knee_tz.r    =  MAj_r(mai(9).mus.l',9);
        MAj.ankle.r      =  MAj_r(mai(10).mus.l',10);
        MAj.subt.r       =  MAj_r(mai(11).mus.l',11);
        % Both legs
        % In MuscleInfo, we first have the right back muscles (44:46) and 
        % then the left back muscles (47:49). Here we re-organize so that
        % we have first the left muscles and then the right muscles.
        lMTj_lr     = [lMTj_l([1:44,48:50],1);lMTj_r(1:47,1)];
        vMTj_lr     = [vMTj_l([1:44,48:50],1);vMTj_r(1:47,1)];   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get muscle-tendon forces and derive Hill-equilibrium
        [Hilldiffj,FTj,Fcej,Fpassj,Fisoj,vMmaxj,massMj] = ...
            f_forceEquilibrium_FtildeState_all_tendon(akj(:,j+1),...
                FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
                lMTj_lr,vMTj_lr,tensions,aTendon,shift); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get metabolic energy rate if in the cost function   
        if W.E ~= 0    
            % Get muscle fiber lengths
            [~,lMtildej] = f_FiberLength_TendonForce_tendon(...
                FTtildekj_nsc(:,j+1),lMTj_lr,aTendon,shift); 
            % Get muscle fiber velocities
            [vMj,~] = f_FiberVelocity_TendonForce_tendon(...
                FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
                lMTj_lr,vMTj_lr,aTendon,shift);
            % Get metabolic energy rate
            if mE == 0 % Bhargava et al. (2004)
                [e_totj,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                    akj(:,j+1),akj(:,j+1),lMtildej,...
                    vMj,Fcej,Fpassj,massMj,pctsts,Fisoj,...
                    MTparameters_m(1,:)',body_mass,10);
            elseif mE == 1 % Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2003all(...
                    akj(:,j+1),akj(:,j+1),lMtildej,...
                    vMtildeUmbj,vMj,Fcej,massMj,...
                    pctsts,vMmaxj,Fisoj,body_mass,10);
            elseif mE == 2 % Umberger (2010)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2010all(...
                    akj(:,j+1),akj(:,j+1),lMtildej,...
                    vMtildeUmbj,vMj,Fcej,massMj,...
                    pctsts,vMmaxj,Fisoj,body_mass,10);  
            elseif mE == 3 % Uchida et al. (2016)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2016all(...
                    akj(:,j+1),akj(:,j+1),lMtildej,...
                    vMtildeUmbj,vMj,Fcej,massMj,...
                    pctsts,vMmaxj,Fisoj,body_mass,10);  
            elseif mE == 4 % Umberger (2010) treating muscle lengthening 
                % heat rate as Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2010all_hl(...
                    akj(:,j+1),akj(:,j+1),lMtildej,...
                    vMtildeUmbj,vMj,Fcej,massMj,...
                    pctsts,vMmaxj,Fisoj,body_mass,10); 
            elseif mE == 5 % Umberger (2010) treating negative mechanical 
                % work as Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbj = vMj./(MTparameters_m(2,:)');
                [e_totj,~,~,~,~] = fgetMetabolicEnergySmooth2010all_neg(...
                    akj(:,j+1),akj(:,j+1),lMtildej,...
                    vMtildeUmbj,vMj,Fcej,massMj,...
                    pctsts,vMmaxj,Fisoj,body_mass,10); 
            end                
        end        
        %%%% just to debug
        Tj_all(:,j)=Tj;
        lMtildej_all(:,j)=lMtildej;
        lMTj_lr_all(:,j)=lMTj_lr;
        FTj_all(:,j)=FTj;
        MAj_knee_ty_r_all(:,j)=MAj.knee_ty.r;
        
        % Get passive torques
        Tau_passj.hip.flex.l    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Qskj_nsc(jointi.hip_flex.l,j+1),...
            Qdotskj_nsc(jointi.hip_flex.l,j+1));
        Tau_passj.hip.flex.r    = f_PassiveMoments(k_pass.hip.flex,...
            theta.pass.hip.flex,Qskj_nsc(jointi.hip_flex.r,j+1),...
            Qdotskj_nsc(jointi.hip_flex.r,j+1));
        Tau_passj.hip.add.l     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Qskj_nsc(jointi.hip_add.l,j+1),...
            Qdotskj_nsc(jointi.hip_add.l,j+1));
        Tau_passj.hip.add.r     = f_PassiveMoments(k_pass.hip.add,...
            theta.pass.hip.add,Qskj_nsc(jointi.hip_add.r,j+1),...
            Qdotskj_nsc(jointi.hip_add.r,j+1));
        Tau_passj.hip.rot.l     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Qskj_nsc(jointi.hip_rot.l,j+1),...
            Qdotskj_nsc(jointi.hip_rot.l,j+1));
        Tau_passj.hip.rot.r     = f_PassiveMoments(k_pass.hip.rot,...
            theta.pass.hip.rot,Qskj_nsc(jointi.hip_rot.r,j+1),...
            Qdotskj_nsc(jointi.hip_rot.r,j+1)); 
        if Options.dampingInKneeSec  
            Tau_passj.knee_add.r    = f_passiveMoments_kneeintmom(Qskj_nsc(jointi.knee_add.r,j+1),Qdotskj_nsc(jointi.knee_add.r,j+1));
            Tau_passj.knee_rot.r    = f_passiveMoments_kneeintmom(Qskj_nsc(jointi.knee_rot.r,j+1),Qdotskj_nsc(jointi.knee_rot.r,j+1));
        else
            Tau_passj.knee_add.r    = f_passiveMoments_kneeintmom(Qskj_nsc(jointi.knee_add.r,j+1));
            Tau_passj.knee_rot.r    = f_passiveMoments_kneeintmom(Qskj_nsc(jointi.knee_rot.r,j+1));
        end
        Tau_passj.knee_flex.l   = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Qskj_nsc(jointi.knee_flex.l,j+1),...
            Qdotskj_nsc(jointi.knee_flex.l,j+1));
        if Options.dampingInKneeSec
            Tau_passj.knee_tx.r     = f_passiveForce_kneeintf(Qskj_nsc(jointi.knee_tx.r,j+1),Qdotskj_nsc(jointi.knee_tx.r,j+1));
            Tau_passj.knee_ty.r     = f_passiveForce_kneeintf(                             0,Qdotskj_nsc(jointi.knee_ty.r,j+1));
            Tau_passj.knee_tz.r     = f_passiveForce_kneeintf(Qskj_nsc(jointi.knee_tz.r,j+1),Qdotskj_nsc(jointi.knee_tz.r,j+1));
        else
            Tau_passj.knee_tx.r     = f_passiveForce_kneeintf(Qskj_nsc(jointi.knee_tx.r,j+1));
            Tau_passj.knee_ty.r     = 0;
            Tau_passj.knee_tz.r     = f_passiveForce_kneeintf(Qskj_nsc(jointi.knee_tz.r,j+1));
        end
        Tau_passj.knee_flex.r        = f_PassiveMoments(k_pass.knee,...
            theta.pass.knee,Qskj_nsc(jointi.knee_flex.r,j+1),...
            Qdotskj_nsc(jointi.knee_flex.r,j+1));
        Tau_passj.ankle.l       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Qskj_nsc(jointi.ankle.l,j+1),...
            Qdotskj_nsc(jointi.ankle.l,j+1));
        Tau_passj.ankle.r       = f_PassiveMoments(k_pass.ankle,...
            theta.pass.ankle,Qskj_nsc(jointi.ankle.r,j+1),...
            Qdotskj_nsc(jointi.ankle.r,j+1));        
        Tau_passj.subt.l       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Qskj_nsc(jointi.subt.l,j+1),...
            Qdotskj_nsc(jointi.subt.l,j+1));
        Tau_passj.subt.r       = f_PassiveMoments(k_pass.subt,...
            theta.pass.subt,Qskj_nsc(jointi.subt.r,j+1),...
            Qdotskj_nsc(jointi.subt.r,j+1));        
        Tau_passj.trunk.ext     = f_PassiveMoments(k_pass.trunk.ext,...
            theta.pass.trunk.ext,Qskj_nsc(jointi.trunk.ext,j+1),...
            Qdotskj_nsc(jointi.trunk.ext,j+1));
        Tau_passj.trunk.ben     = f_PassiveMoments(k_pass.trunk.ben,...
            theta.pass.trunk.ben,Qskj_nsc(jointi.trunk.ben,j+1),...
            Qdotskj_nsc(jointi.trunk.ben,j+1));
        Tau_passj.trunk.rot     = f_PassiveMoments(k_pass.trunk.rot,...
            theta.pass.trunk.rot,Qskj_nsc(jointi.trunk.rot,j+1),...
            Qdotskj_nsc(jointi.trunk.rot,j+1)); 
        Tau_passj_all = [Tau_passj.hip.flex.l,Tau_passj.hip.flex.r,...
            Tau_passj.hip.add.l,Tau_passj.hip.add.r,...
            Tau_passj.hip.rot.l,Tau_passj.hip.rot.r,...
            Tau_passj.knee_add.r, Tau_passj.knee_rot.r, ...
            Tau_passj.knee_flex.l,Tau_passj.knee_flex.r,...
            Tau_passj.knee_tx.r, Tau_passj.knee_ty.r,Tau_passj.knee_tz.r, ...     
            Tau_passj.ankle.l,Tau_passj.ankle.r,...
            Tau_passj.subt.l,Tau_passj.subt.r,...
            Tau_passj.trunk.ext,Tau_passj.trunk.ben,...
            Tau_passj.trunk.rot]';
        Tau_passj_all_d(:,j)=Tau_passj_all;

        % Expression for the state derivatives at the collocation points
        Qsp_nsc          = Qskj_nsc*C(:,j+1);
        Qdotsp_nsc       = Qdotskj_nsc*C(:,j+1);            
        FTtildep_nsc    = FTtildekj_nsc*C(:,j+1);
        ap              = akj*C(:,j+1);
        a_ap            = a_akj*C(:,j+1);
        % Append collocation equations
        % Dynamic constraints are scaled using the same scale
        % factors as the ones used to scale the states
        % Activation dynamics (implicit formulation)
        eq_constr{end+1} = (h*vAk_nsc - ap)./scaling.a;
        g_names_coll = [g_names_coll; repmat({'vA_a_dyn'},NMuscle,1)];
        % Contraction dynamics (implicit formulation)     
        eq_constr{end+1} = (h*dFTtildej_nsc(:,j) - FTtildep_nsc)./...
            scaling.FTtilde';
        g_names_coll = [g_names_coll; repmat({'dFTtilde_FTtilde_dyn'},NMuscle,1)];
        % Skeleton dynamics (implicit formulation)               
        qdotj_nsc = Qdotskj_nsc(:,j+1); % velocity
        eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.QsQdots(1:2:end)';
        g_names_coll = [g_names_coll; repmat({'qdot_dyn'},length(qdotj_nsc),1)];
        eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./...
            scaling.QsQdots(2:2:end)';
        g_names_coll = [g_names_coll; repmat({'qdotdot_dyn'},length(Qdotsp_nsc),1)];
        % Arm activation dynamics (explicit formulation)  
        dadtj = f_ArmActivationDynamics(e_ak,a_akj(:,j+1)');
        eq_constr{end+1} = (h*dadtj - a_ap)./scaling.a_a;
        g_names_coll = [g_names_coll; repmat({'da_a_a_dyn'},length(a_ap),1)];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Add contribution to the cost function
        J = J + ...
            W.E          * B(j+1) *(f_J94(e_totj))/body_mass*h + ...        %min metabolic cost
            W.a          * B(j+1) *(f_J94(akj(:,j+1)'))*h + ...             %min muscle activations
            W.Ak         * B(j+1) *(f_J26(Aj(residuals_noarmsi,j)))*h + ... %min joint accelerations
            W.passMom    * B(j+1) *(f_J17(Tau_passj_all([1:10 14:end])))*h + ...     %min passive moments
            W.passFor    * B(j+1) *(f_J3(Tau_passj_all([11:13])))*h + ...           %min passive forces
            W.u          * B(j+1) *(f_J94(vAk))*h + ...                     %min slack control, derivative of muscle activation                  %min slack control
            W.u          * B(j+1) *(f_J94(dFTtildej(:,j)))*h;               %min slack control, derivative of norm. tendon force
    
        
        J = J + W.ArmE* B(j+1) *(f_J8(e_ak))*h;
      
        if nq.arms > 0
            J = J + W.AkArms  * B(j+1) *(f_J8(Aj(armsi,j)))*h;
        end
        
        J = J + W.Qskneesec * B(j+1) *(f_J4(Qskj([jointi.knee_add.r:jointi.knee_tx.r jointi.knee_tz.r],j+1)))*h; %minimize secondary dofs (knee add and rot, knee tx and tz)
        if Options.minBothContactForces
            J = J + W.minBothKCF * B(j+1)*...
                ((Tj(KCFi.L)/(body_weight*5)).^2+(Tj(KCFi.M)/(body_weight*5)).^2)*h;
        end
        if Options.mindiffKCFGRF
            KCF_totalj=Tj_debug(47)+Tj_debug(48);
            J = J + W.diffKCFGRF * B(j+1)*...
                (((Tj_debug(36)-KCF_totalj)/500).^2)*h;
        end
        
        %for debug
        JE=JE+W.E          * B(j+1) *(f_J94(e_totj))/body_mass*h;
        Ja=Ja+W.a          * B(j+1) *(f_J94(akj(:,j+1)'))*h;
        Jak=Jak+W.Ak       * B(j+1) *(f_J26(Aj(residuals_noarmsi,j)))*h;
        JpassMom=JpassMom+W.passMom*B(j+1) *(f_J17(Tau_passj_all([1:10 14:end])))*h;
        JpassFor=JpassFor+W.passFor*B(j+1) *(f_J3(Tau_passj_all([11:13])))*h;
        JvAk=JvAk+W.u      * B(j+1) *(f_J94(vAk))*h;
        JdFTtilde=JdFTtilde+W.u* B(j+1) *(f_J94(dFTtildej(:,j)))*h;
        JArmE=JArmE+W.ArmE* B(j+1) *(f_J8(e_ak))*h;
        JAjarms=JAjarms+W.AkArms  * B(j+1) *(f_J8(Aj(armsi,j)))*h;
        JQskneesec=JQskneesec+W.Qskneesec * B(j+1) *(f_J4(Qskj([jointi.knee_add.r:jointi.knee_tx.r jointi.knee_tz.r],j+1)))*h;
        if Options.minBothContactForces
            JminbothKCF = JminbothKCF + W.minBothKCF * B(j+1)*...
                ((Tj(KCFi.L)/(body_weight*5)).^2+(Tj(KCFi.M)/(body_weight*5)).^2)*h;
        end
        if Options.mindiffKCFGRF
            JdiffKCFGRF = JdiffKCFGRF + W.diffKCFGRF * B(j+1)*...
                (((Tj_debug(36)-KCF_totalj)/500).^2)*h;
        end

        % Add path constraints
        %Pelvis residuals
        eq_constr{end+1} = Tj(ground_pelvisi,1)./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_pelvis'},6,1)];

        % Muscle-driven joint torques for the lower limbs and the trunk
        % Hip flexion, left
        Ft_hip_flex_l   = FTj(mai(1).mus.l',1);
        T_hip_flex_l    = f_T28(MAj.hip_flex.l,Ft_hip_flex_l);
        eq_constr{end+1} = (Tj(jointi.hip_flex.l,1)-(T_hip_flex_l + ...
            Tau_passj.hip.flex.l))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_hip_flex_l'},1,1)];
        %Hip flexion, right
        Ft_hip_flex_r   = FTj(mai(1).mus.r',1);
        T_hip_flex_r    = f_T28(MAj.hip_flex.r,Ft_hip_flex_r);
        eq_constr{end+1} = (Tj(jointi.hip_flex.r,1)- (T_hip_flex_r + ...
            Tau_passj.hip.flex.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_hip_flex_r'},1,1)];
        % Hip adduction left
        Ft_hip_add_l    = FTj(mai(2).mus.l',1);
        T_hip_add_l     = f_T28(MAj.hip_add.l,Ft_hip_add_l);
        eq_constr{end+1} = (Tj(jointi.hip_add.l,1)-(T_hip_add_l + ...
            Tau_passj.hip.add.l))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_hip_add_l'},1,1)];
        % Hip adduction right
        Ft_hip_add_r    = FTj(mai(2).mus.r',1);
        T_hip_add_r     = f_T28(MAj.hip_add.r,Ft_hip_add_r);
        eq_constr{end+1} = (Tj(jointi.hip_add.r,1)-(T_hip_add_r + ...
            Tau_passj.hip.add.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_hip_add_r'},1,1)];
        %Hip rotation left
        Ft_hip_rot_l    = FTj(mai(3).mus.l',1);
        T_hip_rot_l     = f_T28(MAj.hip_rot.l,Ft_hip_rot_l);
        eq_constr{end+1} = (Tj(jointi.hip_rot.l,1)-(T_hip_rot_l + ...
            Tau_passj.hip.rot.l))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_hip_rot_l'},1,1)];
        %Hip rotation right
        Ft_hip_rot_r    = FTj(mai(3).mus.r',1);
        T_hip_rot_r     = f_T28(MAj.hip_rot.r,Ft_hip_rot_r);
        eq_constr{end+1} = (Tj(jointi.hip_rot.r,1)-(T_hip_rot_r + ...
            Tau_passj.hip.rot.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_hip_rot_r'},1,1)];
        

        %Knee add, right
        Ft_knee_add_r   = FTj(mai(4).mus.r',1);
        T_knee_add_r    = f_T13(MAj.knee_add.r,Ft_knee_add_r);
        eq_constr{end+1} = (Tj(jointi.knee_add.r,1)-(T_knee_add_r + ...
           Tau_passj.knee_add.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_knee_add_r'},1,1)];
        %Knee rot, right
        Ft_knee_rot_r   = FTj(mai(5).mus.r',1);
        T_knee_rot_r    = f_T13(MAj.knee_rot.r,Ft_knee_rot_r);
        eq_constr{end+1} = (Tj(jointi.knee_rot.r,1)-(T_knee_rot_r + ...
            Tau_passj.knee_rot.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_knee_rot_r'},1,1)];
        %Knee flexion, left
        Ft_knee_flex_l  = FTj(mai(6).mus.l',1);
        T_knee_flex_l   = f_T13(MAj.knee_flex.l,Ft_knee_flex_l);
        eq_constr{end+1} = (Tj(jointi.knee_flex.l,1)-(T_knee_flex_l + ...
            Tau_passj.knee_flex.l))./scaling.T(1); 
        g_names_coll = [g_names_coll; repmat({'T_knee_flex_l'},1,1)];
        %Knee flexion, right
        Ft_knee_flex_r= FTj(mai(6).mus.r',1);
        T_knee_flex_r   = f_T13(MAj.knee_flex.r,Ft_knee_flex_r);   
        eq_constr{end+1} = (Tj(jointi.knee_flex.r,1)-(T_knee_flex_r + ...
            Tau_passj.knee_flex.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_knee_flex_r'},1,1)];
        %Knee, tx
        Ft_knee_tx_r    = FTj(mai(7).mus.r',1);
        T_knee_tx_r     = f_T13(MAj.knee_tx.r,Ft_knee_tx_r);
        eq_constr{end+1} = (Tj(jointi.knee_tx.r,1) - (T_knee_tx_r + ...
            Tau_passj.knee_tx.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_knee_tx_r'},1,1)];
        % Knee ty, right
        Ft_knee_ty_r    = FTj(mai(8).mus.r',1);
        T_knee_ty_r     = f_T13(MAj.knee_ty.r,Ft_knee_ty_r);
        eq_constr{end+1} = (Tj(jointi.knee_ty.r,1) - (T_knee_ty_r+...
            Tau_passj.knee_ty.r))./scaling.T(1); 
        g_names_coll = [g_names_coll; repmat({'T_knee_ty_r'},1,1)];
        % Knee tz, right
        Ft_knee_tz_r    = FTj(mai(9).mus.r',1);
        T_knee_tz_r     = f_T13(MAj.knee_tz.r,Ft_knee_tz_r);
        eq_constr{end+1} = (Tj(jointi.knee_tz.r,1)- (T_knee_tz_r + ...
            Tau_passj.knee_tz.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_knee_tz_r'},1,1)];
        % Ankle, left
        Ft_ankle_l      = FTj(mai(10).mus.l',1);
        T_ankle_l       = f_T12(MAj.ankle.l,Ft_ankle_l);
        eq_constr{end+1} = (Tj(jointi.ankle.l,1)-(T_ankle_l + ...
            Tau_passj.ankle.l))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_ankle_l'},1,1)];
        % Ankle, right
        Ft_ankle_r      = FTj(mai(10).mus.r',1);
        T_ankle_r       = f_T12(MAj.ankle.r,Ft_ankle_r);
        eq_constr{end+1} = (Tj(jointi.ankle.r,1)-(T_ankle_r + ...
            Tau_passj.ankle.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_ankle_r'},1,1)];
        % Subtalar, left
        Ft_subt_l       = FTj(mai(11).mus.l',1);
        T_subt_l        = f_T12(MAj.subt.l,Ft_subt_l);
        eq_constr{end+1} = (Tj(jointi.subt.l,1)-(T_subt_l + ...
            Tau_passj.subt.l))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_subt_l'},1,1)];
        % Subtalar, right
        Ft_subt_r       = FTj(mai(11).mus.r',1);
        T_subt_r        = f_T12(MAj.subt.r,Ft_subt_r);
        eq_constr{end+1} = (Tj(jointi.subt.r,1)-(T_subt_r + ...
            Tau_passj.subt.r))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_subt_r'},1,1)];
        % Lumbar extension
        Ft_trunk_ext    = FTj([mai(12).mus.l,mai(12).mus.r]',1);
        T_trunk_ext     = f_T6(MAj.trunk_ext,Ft_trunk_ext);
        eq_constr{end+1} = (Tj(jointi.trunk.ext,1)-(T_trunk_ext +...
            Tau_passj.trunk.ext))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_lumb_ext'},1,1)];
        % Lumbar bending
        Ft_trunk_ben    = FTj([mai(13).mus.l,mai(13).mus.r]',1);
        T_trunk_ben     = f_T6(MAj.trunk_ben,Ft_trunk_ben);
        eq_constr{end+1} = (Tj(jointi.trunk.ben,1)-(T_trunk_ben + ...
            Tau_passj.trunk.ben))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_lumb_ben'},1,1)];
        % Lumbar rotation
        Ft_trunk_rot    = FTj([mai(14).mus.l,mai(14).mus.r]',1);
        T_trunk_rot     = f_T6(MAj.trunk_rot,Ft_trunk_rot);
        eq_constr{end+1} = (Tj(jointi.trunk.rot,1)-(T_trunk_rot +...
            Tau_passj.trunk.rot))./scaling.T(1);
        g_names_coll = [g_names_coll; repmat({'T_lumb_rot'},1,1)];
        % Torque-driven joint torques for the arms
        % Arms
        eq_constr{end+1} = Tj(armsi,1)/scaling.ArmTau - a_akj(:,j+1);
        g_names_coll = [g_names_coll; repmat({'T_arms'},8,1)];

        % Activation dynamics (implicit formulation)
        act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
        act2 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tact);
        ineq_constr1{end+1} = act1;
        ineq_constr2{end+1} = act2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Contraction dynamics (implicit formulation)
        eq_constr{end+1} = Hilldiffj;
        g_names_coll = [g_names_coll; repmat({'Hilldiff_eq'},NMuscle,1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Constraints to prevent parts of the skeleton to penetrate each
        % other.
        % Origins calcaneus (transv plane) at minimum 9 cm from each other.
        ineq_constr3{end+1} = f_J2(Tj(calcOr.r,1) - Tj(calcOr.l,1));          
        % Constraint to prevent the arms to penetrate the skeleton       
        % Origins femurs and ipsilateral hands (transv plane) at minimum 
        % 18 cm from each other.
        ineq_constr4{end+1} = f_J2(Tj(femurOr.r,1) - Tj(handOr.r,1));
        ineq_constr4{end+1} = f_J2(Tj(femurOr.l,1) - Tj(handOr.l,1));
        % Origins tibia (transv plane) at minimum 11 cm from each other.   
        ineq_constr5{end+1} = f_J2(Tj(tibiaOr.r,1) - Tj(tibiaOr.l,1));  

        if Options.defineGRFphases
            %Constraints to set GRF phases
            ineq_constr6{end+1}= GRFrj(2);
            ineq_constr7{end+1} = GRFlj(2);
        end
    end %End loop over collocation points
    eq_constr = vertcat(eq_constr{:});
    ineq_constr1 = vertcat(ineq_constr1{:});
    ineq_constr2 = vertcat(ineq_constr2{:});
    ineq_constr3 = vertcat(ineq_constr3{:});
    ineq_constr4 = vertcat(ineq_constr4{:});
    ineq_constr5 = vertcat(ineq_constr5{:});
    if Options.defineGRFphases
        ineq_constr6 = vertcat(ineq_constr6{:});
        ineq_constr7 = vertcat(ineq_constr7{:});
    else
        ineq_constr6=[];
        ineq_constr7=[];
    end
    
    f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
        Qdotsj,a_ak,a_aj,vAk,e_ak,dFTtildej,Aj},{eq_constr,ineq_constr1,...
        ineq_constr2,ineq_constr3,ineq_constr4,ineq_constr5,ineq_constr6,ineq_constr7,J});
    f_coll_map = f_coll.map(N,parallelMode,NThreads);
    [coll_eq_constr, coll_ineq_constr1, coll_ineq_constr2, coll_ineq_constr3,...
        coll_ineq_constr4, coll_ineq_constr5, coll_ineq_constr6, coll_ineq_constr7, Jall] = f_coll_map(tf,...
        a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
        Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, vA, ...
        e_a, dFTtilde_col, A_col);



    %%%%%%
    %to debug
    f_coll2=Function('f_coll2',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
        Qdotsj,a_ak,a_aj,vAk,e_ak,dFTtildej,Aj},{JE,Ja,Jak,JpassMom,JpassFor,JvAk,JdFTtilde,JArmE,JAjarms,JQskneesec,JminbothKCF,JdiffKCFGRF});
    f_coll_map2 = f_coll2.map(N,parallelMode,NThreads);
    [JEall,Jaall,Jakall,JpassMomall,JpassForall,JvAkall,JdFTtildeall,JArmEall,JAjarmsall,JQskneesecall,JminbothKCFall,JdiffKCFGRFall] = f_coll_map2(tf,...
        a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
        Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, vA, ...
        e_a, dFTtilde_col, A_col);

    f_coll3=Function('f_coll3',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
        Qdotsj,a_ak,a_aj,vAk,e_ak,dFTtildej,Aj},{lMtildej_all,lMTj_lr_all,Tj_all,Tau_passj_all_d,FTj_all,MAj_knee_ty_r_all,Qskj_nsc(:,2:end),Qdotskj_nsc(:,2:end)});
    f_coll_map3 = f_coll3.map(N,parallelMode,NThreads);
    [lMtildej_all_debug,lMTj_lr_all_debug,Tj_all_debug,Tau_passj_all_d_debug,FTj_all_debug,MAj_knee_ty_r_all_debug,Qsj_nsc_debug,Qdotsj_nsc_debug] = f_coll_map3(tf,...
        a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
        Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, vA, ...
        e_a, dFTtilde_col, A_col);
    %%%%%%

    opti.subject_to(coll_eq_constr == 0);
    opti.subject_to(coll_ineq_constr1(:) >= 0);
    opti.subject_to(coll_ineq_constr2(:) <= 1/tact);
    opti.subject_to(0.0081 < coll_ineq_constr3(:) < 4);
    opti.subject_to(0.0324 < coll_ineq_constr4(:) < 4);
    opti.subject_to(0.0121 < coll_ineq_constr5(:) < 4);  

    if Options.defineGRFphases
        Nstance_r=[1:7 25:50];
        Nswing_r=[8:24];
        Nstance_l=[1:32];
        Nswing_l=[33:50];
        coll_ineq_constr6_active=coll_ineq_constr6(:,Nstance_r);
        coll_ineq_constr7_active=coll_ineq_constr7(:,Nstance_l);
        opti.subject_to((coll_ineq_constr6_active(:)-50)/scaling.T(1)>=0);
        opti.subject_to((coll_ineq_constr7_active(:)-50)/scaling.T(1)>=0);
        coll_ineq_constr6_inactive=coll_ineq_constr6(:,Nswing_r);
        coll_ineq_constr7_inactive=coll_ineq_constr7(:,Nswing_l);
        opti.subject_to((coll_ineq_constr6_inactive(:)-50)/scaling.T(1)<=0);
        opti.subject_to((coll_ineq_constr7_inactive(:)-50)/scaling.T(1)<=0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %use unscaled X for continuity constraints
    Qs_nsc = Qs.*(scaling.QsQdots(1:2:end)'*ones(1,N+1));
    Qs_nsc(jointi.knee_ty.r,:) = Qs(jointi.knee_ty.r,:).*scaling.knee_ty.b+scaling.knee_ty.a;
    Qdots_nsc=Qdots.*(scaling.QsQdots(2:2:end)'*ones(1,N+1));
    Qs_col_nsc = Qs_col.*(scaling.QsQdots(1:2:end)'*ones(1,d*N));
    Qs_col_nsc(jointi.knee_ty.r,:)=Qs_col(jointi.knee_ty.r,:).*scaling.knee_ty.b+scaling.knee_ty.a;
    Qdots_col_nsc=Qdots_col.*(scaling.QsQdots(2:2:end)'*ones(1,d*N));
    % Loop over mesh points
    for k=1:N
        % Variables within current mesh interval
        % States      
        akj = [a(:,k), a_col(:,(k-1)*d+1:k*d)]; 
        FTtildekj = [FTtilde(:,k), FTtilde_col(:,(k-1)*d+1:k*d)];
        Qskj_nsc = [Qs_nsc(:,k), Qs_col_nsc(:,(k-1)*d+1:k*d)];
        Qdotskj_nsc = [Qdots_nsc(:,k), Qdots_col_nsc(:,(k-1)*d+1:k*d)];
        a_akj = [a_a(:,k), a_a_col(:,(k-1)*d+1:k*d)];
        % Add equality constraints (next interval starts with end values of 
        % states from previous interval)
        opti.subject_to(a(:,k+1) == akj*D);
        opti.subject_to(FTtilde(:,k+1) == FTtildekj*D); % scaled
        opti.subject_to(Qs_nsc(:,k+1) == Qskj_nsc*D); % scaled
        opti.subject_to(Qdots_nsc(:,k+1) == Qdotskj_nsc*D); % scaled
        opti.subject_to(a_a(:,k+1) == a_akj*D);
    end % End loop over mesh points  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Additional path constraints
    % Periodicity of the states
    % Qs and Qdots
    opti.subject_to(Qs([jointi.pelvis.tilt:jointi.pelvis.rot jointi.pelvis.ty:end],end) - Qs([jointi.pelvis.tilt:jointi.pelvis.rot jointi.pelvis.ty:end],1) == 0);
    opti.subject_to(Qdots(:,end) - Qdots(:,1) == 0);
    opti.subject_to(a(:,end) - a(:,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(:,end) - FTtilde(:,1) == 0);
    % Torque actuator activations
    opti.subject_to(a_a(:,end) - a_a(:,1) == 0);

    % Average speed
    % Provide expression for the distance traveled
    dist_trav_tot = Qs_nsc(jointi.pelvis.tx,end) -  Qs_nsc(jointi.pelvis.tx,1);
    vel_aver_tot = dist_trav_tot/tf; 
    opti.subject_to(vel_aver_tot - v_tgt == 0);
    opti.subject_to(dist_trav_tot >= bounds.dist_trav.lower);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %cost function terms related to pressure values
    if Options.minPressuredev||Options.minPressures
        q_knee=Qs_col_nsc([jointi.knee_flex.r jointi.knee_add.r jointi.knee_rot.r jointi.knee_tx.r jointi.knee_ty.r jointi.knee_tz.r],:);
        out=F2_debug(q_knee);
        pressures_1=full(out(9:8+Options.nfacestib1,:));
        pressures_2=full(out(9+Options.nfacestib1:8+Options.nfacestib1+Options.nfacestib2,:));
        
        
        if Options.minPressuredev
            M1=(sum((pressures_1(:)/1e6).^20))^(1/20);
            M2=(sum((pressures_2(:)/1e6).^20))^(1/20);
            Jall = Jall + (W.minpressuredev * ((log10(M1/(M2+1e-10))-log10(0.5))^2))/N;

            Jminpressdev=(W.minpressuredev * ((log10(M1/(M2+1e-10))-log10(0.5))^2))/N;
            Jminpress=[];
            M=[];
        elseif Options.minPressures
            M=(sum(([pressures_1(:);pressures_2(:)]/1e6).^20))^(1/20);
            Jall = Jall + (W.minpressures * M)/N;
    
            Jminpress=(W.minpressures * M)/N;
            Jminpressdev=[];
            M1=[];
            M2=[];
        end

        %%
    elseif ~Options.minPressuredev&&~Options.minPressures
        Jminpressdev=[];
        M1=[];
        M2=[];
        Jminpress=[];
        M=[];
        pressures_1=[];
        pressures_2=[];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Scale cost function      
    Jall_sc = sum(Jall)/dist_trav_tot;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    % Create NLP solver
    opti.minimize(Jall_sc);
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy      = 'adaptive'; %'adaptive';
    options.ipopt.max_iter = 4000;
    options.ipopt.tol = 1*10^(-tol_ipopt);
    options.ipopt.linear_solver = 'mumps';
    options.ipopt.print_level=5;
    options.ipopt.print_info_string='yes';
    opti.solver('ipopt', options);  
    
    % Create and save diary
    p = mfilename('fullpath');
    [~,namescript,~] = fileparts(p);
    pathresults = [pathRepo,'/Results'];
    if ~(exist([pathresults,'/',namescript],'dir')==7)
        mkdir(pathresults,namescript);
    end
    if (exist([pathresults,'/',namescript,'/D',savename],'file')==2)
        delete ([pathresults,'/',namescript,'/D',savename])
    end 
    diary([pathresults,'/',namescript,'/D',savename]);  
    % Data-informed (full solution at closest speed) initial guess    
    if IGm == 4   
        error('Initial guess not supported')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve problem
    % Opti does not use bounds on variables but constraints. This function
    % adjusts for that. 
    [w_opt,stats,g_opt,lambda_x,lambda_g] = solve_NLPSOL(opti,options);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    diary off
    % Extract results 
    % Create setup
    setup.tolerance.ipopt = tol_ipopt;
    setup.bounds = bounds;
    setup.scaling = scaling;
    setup.guess = guess;
    % Save results and set    
    save([pathresults,'/',namescript,'/w',savename],'w_opt');
    save([pathresults,'/',namescript,'/g',savename],'g_opt');
    save([pathresults,'/',namescript,'/lambda_x',savename],'lambda_x');
    save([pathresults,'/',namescript,'/lambda_g',savename],'lambda_g');
    save([pathresults,'/',namescript,'/stats',savename],'stats');

end

%% Analyze results
if analyseResults
    %% Load results
    if loadResults
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathRepo,'/Results'];
        load([pathresults,'/',namescript,'/w',savename]);
        load([pathresults,'/',namescript,'/stats',savename]);
    end 
    NParameters = 1;    
    tf_opt = w_opt(1:NParameters);
    starti = NParameters+1;
    a_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    a_col_opt = reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    FTtilde_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    FTtilde_col_opt =reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    Qs_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    Qs_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);
    Qdots_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    Qdots_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);    
    a_a_opt = reshape(w_opt(starti:starti+nq.arms*(N+1)-1),nq.arms,N+1)';
    starti = starti + nq.arms*(N+1);
    a_a_col_opt = reshape(w_opt(starti:starti+nq.arms*(d*N)-1),nq.arms,d*N)';
    starti = starti + nq.arms*(d*N);
    vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
    starti = starti + NMuscle*N;
    e_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
    starti = starti + nq.arms*N;   
    dFTtilde_col_opt=reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    qdotdot_col_opt =reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,(d*N))';
    starti = starti + nq.all*(d*N);
    if starti - 1 ~= length(w_opt)
        error('error when extracting results')
    end
    % Combine results at mesh and collocation points
    a_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    a_mesh_col_opt(1:(d+1):end,:)= a_opt;
    FTtilde_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    FTtilde_mesh_col_opt(1:(d+1):end,:)= FTtilde_opt;
    Qs_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    Qs_mesh_col_opt(1:(d+1):end,:)= Qs_opt;
    Qdots_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    Qdots_mesh_col_opt(1:(d+1):end,:)= Qdots_opt;
    a_a_mesh_col_opt=zeros(N*(d+1)+1,nq.arms);
    a_a_mesh_col_opt(1:(d+1):end,:)= a_a_opt;
    for k=1:N
        rangei = k*(d+1)-(d-1):k*(d+1);
        rangebi = (k-1)*d+1:k*d;
        a_mesh_col_opt(rangei,:) = a_col_opt(rangebi,:);
        FTtilde_mesh_col_opt(rangei,:) = FTtilde_col_opt(rangebi,:);
        Qs_mesh_col_opt(rangei,:) = Qs_col_opt(rangebi,:);
        Qdots_mesh_col_opt(rangei,:) = Qdots_col_opt(rangebi,:);
        a_a_mesh_col_opt(rangei,:) = a_a_col_opt(rangebi,:);
    end

    %% Unscale results
    % States at mesh points
    % Qs (1:N-1)
    q_opt_unsc.rad = Qs_opt(1:end-1,:).*repmat(...
        scaling.Qs,size(Qs_opt(1:end-1,:),1),1); 
    q_opt_unsc.rad(:,jointi.knee_ty.r) = Qs_opt(1:end-1,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
    % Convert in degrees
    q_opt_unsc.deg = q_opt_unsc.rad;
    q_opt_unsc.deg(:,[1:3,7:16 20:end]) = q_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qs (1:N)
    q_opt_unsc_all.rad = Qs_opt(1:end,:).*repmat(...
        scaling.Qs,size(Qs_opt(1:end,:),1),1); 
    q_opt_unsc_all.rad(:,jointi.knee_ty.r)=Qs_opt(1:end,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
    % Convert in degrees
    q_opt_unsc_all.deg = q_opt_unsc_all.rad;
    q_opt_unsc_all.deg(:,[1:3,7:16 20:end]) = ...
        q_opt_unsc_all.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qdots (1:N-1)
    qdot_opt_unsc.rad = Qdots_opt(1:end,:).*repmat(...
        scaling.Qdots,size(Qdots_opt(1:end,:),1),1);
    % Convert in degrees
    qdot_opt_unsc.deg = qdot_opt_unsc.rad;
    qdot_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
    qdot_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qdots (1:N)
    qdot_opt_unsc_all.rad =Qdots_opt.*repmat(scaling.Qdots,size(Qdots_opt,1),1); 
    % Muscle activations (1:N-1)
    a_opt_unsc = a_opt(1:end-1,:).*repmat(...
        scaling.a,size(a_opt(1:end-1,:),1),size(a_opt,2));
    a_opt_unsc_all = a_opt(1:end,:).*repmat(...
        scaling.a,size(a_opt(1:end,:),1),size(a_opt,2));
    % Muscle-tendon forces (1:N-1)
    FTtilde_opt_unsc = FTtilde_opt(1:end-1,:).*repmat(...
        scaling.FTtilde,size(FTtilde_opt(1:end-1,:),1),1);
    % Muscle-tendon forces
    FTtilde_opt_unsc_all = FTtilde_opt(1:end,:).*repmat(...
        scaling.FTtilde,size(FTtilde_opt(1:end,:),1),1);
    % Arm activations (1:N-1)
    a_a_opt_unsc = a_a_opt(1:end-1,:).*repmat(...
        scaling.a_a,size(a_a_opt(1:end-1,:),1),size(a_a_opt,2));
    % Arm activations (1:N)
    a_a_opt_unsc_all = a_a_opt(1:end,:).*repmat(...
        scaling.a_a,size(a_a_opt(1:end,:),1),size(a_a_opt,2));
    
    % Controls at mesh points
    % Time derivative of muscle activations (states)
    vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
    tact = 0.015;
    tdeact = 0.06;
    % Get muscle excitations from time derivative of muscle activations
    e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
    % Arm excitations
    e_a_opt_unsc = e_a_opt.*repmat(scaling.e_a,size(e_a_opt,1),...
        size(e_a_opt,2));
    % State and Controls at collocation points
    % Qs a
    q_col_opt=Qs_col_opt;
    qdot_col_opt=Qdots_col_opt;
    q_col_opt_unsc.rad = q_col_opt(1:end,:).*repmat(...
    scaling.Qs,size(q_col_opt(1:end,:),1),1); 
    q_col_opt_unsc.rad(:,jointi.knee_ty.r) = q_col_opt(1:end,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
    % Convert in degrees
    q_col_opt_unsc.deg = q_col_opt_unsc.rad;
    q_col_opt_unsc.deg(:,[1:3,7:16 20:end]) = q_col_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qdots
    qdot_col_opt_unsc.rad = qdot_col_opt(1:end,:).*repmat(...
        scaling.Qdots,size(qdot_col_opt(1:end,:),1),1);
    % Convert in degrees
    qdot_col_opt_unsc.deg = qdot_col_opt_unsc.rad;
    qdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
        qdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Muscle activations
    a_col_opt_unsc = a_col_opt.*repmat(...
        scaling.a,size(a_col_opt,1),size(a_col_opt,2));
    %Muscle-tendon forces
    FTtilde_col_opt_unsc= FTtilde_col_opt.*repmat(...
        scaling.FTtilde,size(FTtilde_col_opt,1),1);
    % Arm activations
    a_a_col_opt_unsc=a_a_col_opt.*repmat(...
        scaling.a_a,size(a_a_col_opt,1),size(a_a_col_opt,2));
    % "Slack" controls at collocation points   
    % Time derivative of Qdots
    qdotdot_col_opt_unsc.rad = ...
        qdotdot_col_opt.*repmat(scaling.Qdotdots,size(qdotdot_col_opt,1),1);
    % Convert in degrees
    qdotdot_col_opt_unsc.deg = qdotdot_col_opt_unsc.rad;
    qdotdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
        qdotdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Time derivative of muscle-tendon forces
    dFTtilde_col_opt_unsc = dFTtilde_col_opt.*repmat(...
        scaling.dFTtilde,size(dFTtilde_col_opt,1),size(dFTtilde_col_opt,2));

    %% Time grid    
    % Mesh points
    tgrid = linspace(0,tf_opt,N+1);
    dtime = zeros(1,d+1);
    for i=1:4
        dtime(i)=tau_root(i)*((tf_opt-0)/N);
    end
    % Mesh points and collocation points
    tgrid_ext = zeros(1,(d+1)*N+1);
    for i=1:N
        tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
    end
    tgrid_ext(end)=tf_opt(end); 

    %% Joint torques, ground reaction forces and moments at opt solution
    Xk_Qs_Qdots_opt                 = zeros(N,2*size(q_opt_unsc.rad,2));
    Xk_Qs_Qdots_opt(:,1:2:end)      = q_opt_unsc_all.rad(2:end,:);
    Xk_Qs_Qdots_opt(:,2:2:end)      = qdot_opt_unsc_all.rad(2:end,:);
    Xk_Qdotdots_opt                 = qdotdot_col_opt_unsc.rad(3:3:end,:); %to be removed!! 

    Xk_Qs_Qdots_col_opt             = zeros(N*d,2*size(q_opt_unsc.rad,2));
    Xk_Qs_Qdots_col_opt(:,1:2:end)  = q_col_opt_unsc.rad(1:N*d,:);
    Xk_Qs_Qdots_col_opt(:,2:2:end)  = qdot_col_opt_unsc.rad(1:N*d,:);
    Xk_Qdotdots_col_opt             = qdotdot_col_opt_unsc.rad;  

    for i = 1:N
        if Options.KCFasinputstoExternalFunction
            [~, KCF_opt_unsc_k(i,:), GRF_opt_r_k(i,:), GRF_opt_l_k(i,:), calcn_opt_r_k(i,:), calcn_opt_l_k(i,:),pressures_opt_1_k(i,:),pressures_opt_2_k(i,:)]=...
                ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_opt(i,:)',Xk_Qdotdots_opt(i,:)',{F2_skeletal_debug F2_debug},deri,residualsi,KCFi,armsi,scaling,tol_ipopt,a_a_opt_unsc(i,:),Options.nfacestib1,Options.nfacestib2,Options,jointi);
        else
            [~, KCF_opt_unsc_k(i,:), GRF_opt_r_k(i,:), GRF_opt_l_k(i,:), calcn_opt_r_k(i,:), calcn_opt_l_k(i,:),pressures_opt_1_k(i,:),pressures_opt_2_k(i,:)]=...
                ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_opt(i,:)',Xk_Qdotdots_opt(i,:)',F2_debug,deri,residualsi,KCFi,armsi,scaling,tol_ipopt,a_a_opt_unsc(i,:),Options.nfacestib1,Options.nfacestib2,Options,jointi);
        end


        for j=1:d
            if Options.KCFasinputstoExternalFunction
                [Tauk_out((i-1)*d+j,:), KCF_opt_unsc((i-1)*d+j,:), GRF_opt_r((i-1)*d+j,:), GRF_opt_l((i-1)*d+j,:), calcn_opt_r((i-1)*d+j,:), calcn_opt_l((i-1)*d+j,:),pressures_opt_1((i-1)*d+j,:),pressures_opt_2((i-1)*d+j,:)]=...
                    ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_col_opt((i-1)*d+j,:)',Xk_Qdotdots_col_opt((i-1)*d+j,:)',...
                    {F2_skeletal_debug F2_debug},deri,residualsi,KCFi,armsi,scaling,tol_ipopt,a_a_col_opt_unsc((i-1)*d+j,:),Options.nfacestib1,Options.nfacestib2,Options,jointi);
            else
                [Tauk_out((i-1)*d+j,:), KCF_opt_unsc((i-1)*d+j,:), GRF_opt_r((i-1)*d+j,:), GRF_opt_l((i-1)*d+j,:), calcn_opt_r((i-1)*d+j,:), calcn_opt_l((i-1)*d+j,:),pressures_opt_1((i-1)*d+j,:),pressures_opt_2((i-1)*d+j,:)]=...
                    ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_col_opt((i-1)*d+j,:)',Xk_Qdotdots_col_opt((i-1)*d+j,:)',...
                    F2_debug,deri,residualsi,KCFi,armsi,scaling,tol_ipopt,a_a_col_opt_unsc((i-1)*d+j,:),Options.nfacestib1,Options.nfacestib2,Options,jointi);
            end
        end
    end

    %% Muscle data / Hill's model variables
    for i=1:N
        for j=1:d
            % Get muscle-tendon lengths, velocities, and moment arms
            % Left leg
            qin_l = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2-1,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2-1,1), ...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2-1,1), 0,0,...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2-1,1),0,0.042,0 ...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2-1,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2-1,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2-1,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2-1,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2-1,1)];  
            qdotin_l = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2,1),0,0,...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2,1),0,0,0,...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2,1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2,1)];  
            [lMTk_l,vMTk_l,MA_l] = f_lMT_vMT_dM(qin_l,qdotin_l);    
            MA.hip_flex.l   =  MA_l(mai(1).mus.l',1);
            MA.hip_add.l    =  MA_l(mai(2).mus.l',2);
            MA.hip_rot.l    =  MA_l(mai(3).mus.l',3);
            MA.knee_flex.l  =  MA_l(mai(6).mus.l',6);
            MA.ankle.l      =  MA_l(mai(10).mus.l',10);  
            MA.subt.l       =  MA_l(mai(11).mus.l',11); 
            % For the back muscles, we want left and right together: left
            % first, right second. In MuscleInfo, we first have the right
            % muscles (45:47) and then the left muscles (48:50). Since the back
            % muscles only depend on back dofs, we do not care if we extract
            % them "from the left or right leg" so here we just picked left.
            MA.trunk_ext    =  MA_l([48:50,mai(12).mus.l]',12);
            MA.trunk_ben    =  MA_l([48:50,mai(13).mus.l]',13);
            MA.trunk_rot    =  MA_l([48:50,mai(14).mus.l]',14);
            % Right leg
            qin_r = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_ty.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2-1)];  
            qdotin_r = [Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_ty.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2)];      
            [lMTk_r,vMTk_r,MA_r] = f_lMT_vMT_dM(qin_r,qdotin_r);
            % Here we take the indices from left since the vector is 1:47 (is
            % 47?, are right instead of left...)
            MA.hip_flex.r   =  MA_r(mai(1).mus.l',1);
            MA.hip_add.r    =  MA_r(mai(2).mus.l',2);
            MA.hip_rot.r    =  MA_r(mai(3).mus.l',3);
            MA.knee_add.r   =  MA_r(mai(4).mus.l',4);
            MA.knee_rot.r   =  MA_r(mai(5).mus.l',5);
            MA.knee_flex.r  =  MA_r(mai(6).mus.l',6);
            MA.knee_tx.r    =  MA_r(mai(7).mus.l',7);
            MA.knee_ty.r    =  MA_r(mai(8).mus.l',8);
            MA.knee_tz.r    =  MA_r(mai(9).mus.l',9);
            MA.ankle.r      =  MA_r(mai(10).mus.l',10);
            MA.subt.r       =  MA_r(mai(11).mus.l',11);
            MA_mat((i-1)*d+j).MA=MA; %store matrix of MA for later calculations
            % Both legs
            % In MuscleInfo, we first have the right back muscles (45:47) and 
            % then the left back muscles (48:50). Here we re-organize so that
            % we have first the left muscles and then the right muscles.
            lMTk_lr     = [lMTk_l([1:44,48:50],1);lMTk_r(1:47,1)];
            vMTk_lr     = [vMTk_l([1:44,48:50],1);vMTk_r(1:47,1)]; 
    
            [Hilldiffk_opt_aux,FTk_opt_aux,~,~,~] =  f_forceEquilibrium_FtildeState_all_tendon(...
                        a_col_opt((i-1)*d+j,:),FTtilde_col_opt_unsc((i-1)*d+j,:),...
                        dFTtilde_col_opt_unsc((i-1)*d+j,:), lMTk_lr,vMTk_lr,tensions,aTendon,shift);
                    Hilldiffk_opt((i-1)*d+j,:)=full(Hilldiffk_opt_aux);
                    FTk_opt((i-1)*d+j,:)=full(FTk_opt_aux);
                    


            [Hilldiffk_opt_aux2(:,(i-1)*d+j), FTk_opt_aux2(:,(i-1)*d+j), Fce_opt_aux2(:,(i-1)*d+j), Fiso, vMmax, Fpe_opt_aux2(:,(i-1)*d+j), lMtilde_opt_aux2(:,(i-1)*d+j), FMvtilde_aux2(:,(i-1)*d+j), vMtilde_aux2(:,(i-1)*d+j)] = ...
            ForceEquilibrium_FtildeState_GC_test(a_col_opt((i-1)*d+j,:),FTtilde_col_opt_unsc((i-1)*d+j,:),dFTtilde_col_opt_unsc((i-1)*d+j,:),full(lMTk_lr)',full(vMTk_lr)',MTparameters_m,Fvparam,...
                Fpparam,Faparam);    
            
                MAr_all((i-1)*d+j,:,:)=full(MA_r);
                MAl_all((i-1)*d+j,:,:)=full(MA_l);
        
            % Get optimal passive torques
            Tau_pass_opt.hip.flex.l((i-1)*d+j,:)    = full(f_PassiveMoments(k_pass.hip.flex,...
                theta.pass.hip.flex,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.l*2)));
            Tau_pass_opt.hip.flex.r((i-1)*d+j,:)    = full(f_PassiveMoments(k_pass.hip.flex,...
                theta.pass.hip.flex,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_flex.r*2)));
            Tau_pass_opt.hip.add.l((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.add,...
                theta.pass.hip.add,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.l*2)));
            Tau_pass_opt.hip.add.r((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.add,...
                theta.pass.hip.add,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_add.r*2)));
            Tau_pass_opt.hip.rot.l((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.rot,...
                theta.pass.hip.rot,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.l*2)));
            Tau_pass_opt.hip.rot.r((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.hip.rot,...
                theta.pass.hip.rot,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.hip_rot.r*2))); 
            if Options.dampingInKneeSec
                Tau_pass_opt.knee_add.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2)));
                Tau_pass_opt.knee_rot.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2)));
            else
                Tau_pass_opt.knee_add.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_add.r*2-1)));
                Tau_pass_opt.knee_rot.r((i-1)*d+j,:)    = full(f_passiveMoments_kneeintmom(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_rot.r*2-1)));
            
            end
            Tau_pass_opt.knee_flex.l((i-1)*d+j,:)   = full(f_PassiveMoments(k_pass.knee,...
                theta.pass.knee,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.l*2)));
            if Options.dampingInKneeSec
                Tau_pass_opt.knee_tx.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2)));
                Tau_pass_opt.knee_ty.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(                                                  0,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_ty.r*2)));
                Tau_pass_opt.knee_tz.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2-1),Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2)));
            else
                Tau_pass_opt.knee_tx.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tx.r*2-1)));
                Tau_pass_opt.knee_ty.r((i-1)*d+j,:)     = 0;
                Tau_pass_opt.knee_tz.r((i-1)*d+j,:)     = full(f_passiveForce_kneeintf(Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_tz.r*2-1)));
            end
            Tau_pass_opt.knee_flex.r((i-1)*d+j,:)        = full(f_PassiveMoments(k_pass.knee,...
                theta.pass.knee,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.knee_flex.r*2)));
            Tau_pass_opt.ankle.l((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.ankle,...
                theta.pass.ankle,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.l*2)));
            Tau_pass_opt.ankle.r((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.ankle,...
                theta.pass.ankle,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.ankle.r*2)));        
            Tau_pass_opt.subt.l((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.subt,...
                theta.pass.subt,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.l*2)));
            Tau_pass_opt.subt.r((i-1)*d+j,:)       = full(f_PassiveMoments(k_pass.subt,...
                theta.pass.subt,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.subt.r*2)));        
            Tau_pass_opt.trunk.ext((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.trunk.ext,...
                theta.pass.trunk.ext,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ext*2)));
            Tau_pass_opt.trunk.ben((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.trunk.ben,...
                theta.pass.trunk.ben,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.ben*2)));
            Tau_pass_opt.trunk.rot((i-1)*d+j,:)     = full(f_PassiveMoments(k_pass.trunk.rot,...
                theta.pass.trunk.rot,Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2-1),...
                Xk_Qs_Qdots_col_opt((i-1)*d+j,jointi.trunk.rot*2))); 


            %% Recalculate path constraints
            % % Pelvis residuals (we want them to be zero at the optimal solution)
            eq_constr_opt((i-1)*d+j, 1:6) = (Tauk_out((i-1)*d+j,ground_pelvisi));

            % Hip flexion right
            Ft_hip_flex_r_opt  = FTk_opt((i-1)*d+j,mai(1).mus.r');
            T_hip_flex_r_opt   = full(f_T28(MA.hip_flex.r,Ft_hip_flex_r_opt));
            if Options.useJointResMom
                eq_constr_opt((i-1)*d+j, jointi.hip_flex.r) = Tauk_out((i-1)*d+j, jointi.hip_flex.r)-(T_hip_flex_r_opt + ...
                    Tau_pass_opt.hip.flex.r((i-1)*d+j,:))-JointMom_res_opt(:,jointi.knee_flex.r)*scaling.JRM_res;
            else
                eq_constr_opt((i-1)*d+j, jointi.hip_flex.r) = Tauk_out((i-1)*d+j, jointi.hip_flex.r)-(T_hip_flex_r_opt + ...
                    Tau_pass_opt.hip.flex.r((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.hip_flex.r)=T_hip_flex_r_opt;

            % Hip flexion left
            Ft_hip_flex_l_opt  = FTk_opt((i-1)*d+j,mai(1).mus.l');
            T_hip_flex_l_opt   = full(f_T28(MA.hip_flex.l,Ft_hip_flex_l_opt));
            if Options.useJointResMom
                eq_constr_opt((i-1)*d+j, jointi.hip_flex.l) = Tauk_out((i-1)*d+j, jointi.hip_flex.l)-(T_hip_flex_l_opt + ...
                    Tau_pass_opt.hip.flex.l((i-1)*d+j,:))-JointMom_res_opt(jointi.knee_flex.l,1)*scaling.JRM_res;
            else
                eq_constr_opt((i-1)*d+j, jointi.hip_flex.l) = Tauk_out((i-1)*d+j, jointi.hip_flex.l)-(T_hip_flex_l_opt + ...
                    Tau_pass_opt.hip.flex.l((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.hip_flex.l)=T_hip_flex_l_opt;

            % 
            % % Hip adduction right
            % 
            % Ft_hip_add_r_opt  = FTk_opt(i,mai(2).mus.r');
            % T_hip_add_r_opt   = full(f_T28(MAk.hip_add.r,Ft_hip_add_r_opt));
            % 
            % 
            % if Options.useJointResMom
            %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
            %         Tau_pass_opt.hip.add.r(i,:))-JointMom_res_opt(jointi.knee_add.r,1)*scaling.JRM_res;
            % else
            %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
            %         Tau_pass_opt.hip.add.r(i,:));
            % end
            % 
            % 
            % 
            % % Hip adduction left
            % 
            % Ft_hip_add_l_opt  = FTk_opt(i,mai(2).mus.l');
            % T_hip_add_l_opt   = full(f_T28(MAk.hip_add.l,Ft_hip_add_l_opt));
            % 
            % 
            % if Options.useJointResMom
            %     eq_constr_opt(i, jointi.hip_add.l) = Tauk_out(i, jointi.hip_add.l)-(T_hip_add_l_opt + ...
            %         Tau_pass_opt.hip.add.l(i,:))-JointMom_res_opt(jointi.knee_add.l,1)*scaling.JRM_res(1);
            % else
            %     eq_constr_opt(i, jointi.hip_add.l) = Tauk_out(i, jointi.hip_add.l)-(T_hip_add_l_opt + ...
            %         Tau_pass_opt.hip.add.l(i,:));
            % end
            % 
            % 
            % 
            % % Hip rotation right
            % 
            % Ft_hip_add_r_opt  = FTk_opt(i,mai(3).mus.r');
            % T_hip_add_r_opt   = full(f_T28(MAk.hip_add.r,Ft_hip_add_r_opt));
            % 
            % 
            % if Options.useJointResMom
            %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
            %         Tau_pass_opt.hip.add.r(i,:))-JointMom_res_opt(jointi.knee_add.r,1)*scaling.JRM_res(2);
            % else
            %     eq_constr_opt(i, jointi.hip_add.r) = Tauk_out(i, jointi.hip_add.r)-(T_hip_add_r_opt + ...
            %         Tau_pass_opt.hip.add.r(i,:));
            % end
            % 
            % 
            % 
            % % Hip rotation left
            % 
            % Ft_hip_rot_l_opt  = FTk_opt(i,mai(3).mus.l');
            % T_hip_rot_l_opt   = full(f_T28(MAk.hip_rot.l,Ft_hip_rot_l_opt));
            % 
            % 
            % if Options.useJointResMom
            %     eq_constr_opt(i, jointi.hip_rot.l) = Tauk_out(i, jointi.hip_rot.l)-(T_hip_rot_l_opt + ...
            %         Tau_pass_opt.hip.rot.l(i,:))-JointMom_res_opt(jointi.knee_rot.l,1)*scaling.JRM_res;
            % else
            %     eq_constr_opt(i, jointi.hip_rot.l) = Tauk_out(i, jointi.hip_rot.l)-(T_hip_rot_l_opt + ...
            %         Tau_pass_opt.hip.rot.l(i,:));
            % end
        
        
            %Knee add, right
            Ft_knee_add_r_opt   = FTk_opt((i-1)*d+j,mai(4).mus.r');
            T_knee_add_r_opt    = full(f_T13(MA.knee_add.r,Ft_knee_add_r_opt));
            if Options.useKCFresiduals
                eq_constr_opt((i-1)*d+j,jointi.knee_add.r) = Tauk_out((i-1)*d+j,jointi.knee_add.r)-(T_knee_add_r_opt + ...
                    Tau_pass_opt.knee_add.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,1).*scaling.KCF_res(1);
            else
                eq_constr_opt((i-1)*d+j,jointi.knee_add.r) = Tauk_out((i-1)*d+j,jointi.knee_add.r,1)-(T_knee_add_r_opt + ...
                   Tau_pass_opt.knee_add.r((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.knee_add.r)=T_knee_add_r_opt;

            %Knee rot, right
            Ft_knee_rot_r_opt   = FTk_opt((i-1)*d+j,mai(5).mus.r');
            T_knee_rot_r_opt    = full(f_T13(MA.knee_rot.r,Ft_knee_rot_r_opt));
            if Options.useKCFresiduals
                eq_constr_opt((i-1)*d+j,jointi.knee_rot.r) = Tauk_out((i-1)*d+j,jointi.knee_rot.r)-(T_knee_rot_r_opt + ...
                    Tau_pass_opt.knee_rot.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,2).*scaling.KCF_res(2);
            else
                eq_constr_opt((i-1)*d+j,jointi.knee_rot.r) = Tauk_out((i-1)*d+j,jointi.knee_rot.r)-(T_knee_rot_r_opt + ...
                    Tau_pass_opt.knee_rot.r((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.knee_rot.r)=T_knee_rot_r_opt;

            %Knee flexion, left
            Ft_knee_flex_l_opt  = FTk_opt((i-1)*d+j,mai(6).mus.l');
            T_knee_flex_l_opt   = full(f_T13(MA.knee_flex.l,Ft_knee_flex_l_opt));
            if Options.useJointResMom
                eq_constr_opt((i-1)*d+j,jointi.knee_flex.l) = Tauk_out((i-1)*d+j,jointi.knee_flex.l)-(T_knee_flex_l_opt + ...
                    Tau_pass_opt.knee_flex.l((i-1)*d+j,:))-JointMom_res_opt(jointi.knee_flex.l,1)*scaling.JRM_res;
            else
                eq_constr_opt((i-1)*d+j,jointi.knee_flex.l)= Tauk_out((i-1)*d+j,jointi.knee_flex.l)-(T_knee_flex_l_opt + ...
                    Tau_pass_opt.knee_flex.l((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.knee_flex.l)=T_knee_flex_l_opt;

            %Knee flexion, right
            Ft_knee_flex_r_opt= FTk_opt((i-1)*d+j,mai(6).mus.r');
            T_knee_flex_r_opt   = full(f_T13(MA.knee_flex.r,Ft_knee_flex_r_opt));
            if Options.useJointResMom
                eq_constr_opt((i-1)*d+j,jointi.knee_flex.r) = Tauk_out((i-1)*d+j,jointi.knee_flex.r)-(T_knee_flex_r_opt + ...
                    Tau_pass_opt.knee_flex.r((i-1)*d+j,:))- JointMom_res_opt(jointi.knee_flex.r,1)*scaling.JRM_res;
            else
                eq_constr_opt((i-1)*d+j,jointi.knee_flex.r) = Tauk_out((i-1)*d+j,jointi.knee_flex.r)-(T_knee_flex_r_opt + ...
                    Tau_pass_opt.knee_flex.r((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.knee_flex.r)=T_knee_flex_r_opt;

            %Knee, tx
            Ft_knee_tx_r_opt    = FTk_opt((i-1)*d+j,mai(7).mus.r');
            T_knee_tx_r_opt     = full(f_T13(MA.knee_tx.r,Ft_knee_tx_r_opt));
            if Options.useKCFresiduals
                eq_constr_opt((i-1)*d+j,jointi.knee_tx.r) = Tauk_out((i-1)*d+j,jointi.knee_tx.r) - (T_knee_tx_r_opt + ...
                    Tau_pass_opt.knee_tx.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,3).*scaling.KCF_res(3);
            else
                eq_constr_opt((i-1)*d+j,jointi.knee_tx.r) = Tauk_out((i-1)*d+j,jointi.knee_tx.r) - (T_knee_tx_r_opt + ...
                    Tau_pass_opt.knee_tx.r((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.knee_tx.r)=T_knee_tx_r_opt;

            % Knee ty, right
            Ft_knee_ty_r_opt    = FTk_opt((i-1)*d+j,mai(8).mus.r');
            T_knee_ty_r_opt((i-1)*d+j,:)     = full(f_T13(MA.knee_ty.r,Ft_knee_ty_r_opt));
            if Options.useKCFresiduals
                eq_constr_opt((i-1)*d+j,jointi.knee_ty.r) = Tauk_out((i-1)*d+j,jointi.knee_ty.r) - (T_knee_ty_r_opt((i-1)*d+j,:)+...
                    Tau_pass_opt.knee_ty.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,4).*scaling.KCF_res(4); %no passive force for knee ty for now
            else
                eq_constr_opt((i-1)*d+j,jointi.knee_ty.r) = Tauk_out((i-1)*d+j,jointi.knee_ty.r) - (T_knee_ty_r_opt((i-1)*d+j,:)+...
                    Tau_pass_opt.knee_ty.r((i-1)*d+j,:)); %no passive force for knee ty for now
            end
            MuscleMoments((i-1)*d+j, jointi.knee_ty.r)=T_knee_ty_r_opt((i-1)*d+j,:);

            % Knee tz, right
            Ft_knee_tz_r_opt    = FTk_opt((i-1)*d+j,mai(9).mus.r');
            T_knee_tz_r_opt     = full(f_T13(MA.knee_tz.r,Ft_knee_tz_r_opt));
            if Options.useKCFresiduals
                eq_constr_opt((i-1)*d+j,jointi.knee_tz.r) = Tauk_out((i-1)*d+j,jointi.knee_tz.r)- (T_knee_tz_r_opt + ...
                    Tau_pass_opt.knee_tz.r((i-1)*d+j,:))+KCF_res_opt((i-1)*d+j,5).*scaling.KCF_res(5);
            else
                eq_constr_opt((i-1)*d+j,jointi.knee_tz.r) = Tauk_out((i-1)*d+j,jointi.knee_tz.r)- (T_knee_tz_r_opt + ...
                    Tau_pass_opt.knee_tz.r((i-1)*d+j,:));
            end
            MuscleMoments((i-1)*d+j, jointi.knee_tz.r)=T_knee_tz_r_opt;

            
            % Torque-driven joint torques for the arms
            % Arms
            eq_constr_opt_arms((i-1)*d+j,:) = Tauk_out((i-1)*d+j,armsi)/scaling.ArmTau - a_a_col_opt_unsc((i-1)*d+j,:);

            %Reconstruct constraints to prevent parts of the skeleton to 
            % penetrate each other.
            % Origins calcaneus (transv plane) at minimum 9 cm from each other.
            ineq_constr3_opt((i-1)*d+j,:) = full(f_J2(calcn_opt_r((i-1)*d+j,[1 3]) - calcn_opt_l((i-1)*d+j,[1 3])));          
            % % Constraint to prevent the arms to penetrate the skeleton       
            % % Origins femurs and ipsilateral hands (transv plane) at minimum 
            % % 18 cm from each other.
            % ineq_constr4_opt((i-1)*d+j,:) = f_J2(Tj(femurOr.r,1) - Tj(handOr.r,1));
            % ineq_constr4{end+1} = f_J2(Tj(femurOr.l,1) - Tj(handOr.l,1));
            % % Origins tibia (transv plane) at minimum 11 cm from each other.   
            % ineq_constr5{end+1} = f_J2(Tj(tibiaOr.r,1) - Tj(tibiaOr.l,1)); 

        end
    end
    eq_constr_opt=eq_constr_opt./scaling.T(1);

    Tau_pass_opt_all = [Tau_pass_opt.hip.flex.l,Tau_pass_opt.hip.flex.r,...
    Tau_pass_opt.hip.add.l,Tau_pass_opt.hip.add.r,...
    Tau_pass_opt.hip.rot.l,Tau_pass_opt.hip.rot.r,...
    Tau_pass_opt.knee_add.r, Tau_pass_opt.knee_rot.r,...
    Tau_pass_opt.knee_flex.l,Tau_pass_opt.knee_flex.r,...
    Tau_pass_opt.knee_tx.r, Tau_pass_opt.knee_ty.r,Tau_pass_opt.knee_tz.r,...
    Tau_pass_opt.ankle.l,Tau_pass_opt.ankle.r,...
    Tau_pass_opt.subt.l,Tau_pass_opt.subt.r,...
    Tau_pass_opt.trunk.ext,Tau_pass_opt.trunk.ben,...
    Tau_pass_opt.trunk.rot]';

    %% Stride length and width  
    % For the stride length we also need the values at the end of the
    % interval so N+1 where states but not controls are defined
    Xk_Qs_Qdots_opt_all = zeros(N+1,2*size(q_opt_unsc_all.rad,2));
    Xk_Qs_Qdots_opt_all(:,1:2:end)  = q_opt_unsc_all.rad;
    Xk_Qs_Qdots_opt_all(:,2:2:end)  = qdot_opt_unsc_all.rad;
    % The stride length is the distance covered by the calcaneus origin
    % Right leg
    dist_r = sqrt(3*(f_J3(calcn_opt_r_k(end,:)-calcn_opt_r_k(1,:))));
    % Left leg
    dist_l = sqrt(3*(f_J3(calcn_opt_l_k(end,:)-calcn_opt_l_k(1,:))));
    % The total stride length is the sum of the right and left stride 
    % lengths after a half gait cycle, since we assume symmetry
    % StrideLength_opt = full(dist_r + dist_l);    
    % The stride width is the medial distance between the calcaneus origins
    StepWidth_opt = full(abs(calcn_opt_r_k(:,3)-calcn_opt_l_k(:,3)));
    stride_width_mean = mean(StepWidth_opt);
    stride_width_std = std(StepWidth_opt);

    %% Assert average speed    
    dist_trav_opt = q_opt_unsc_all.rad(end,jointi.pelvis.tx) - ...
        q_opt_unsc_all.rad(1,jointi.pelvis.tx); % distance traveled
    time_elaps_opt = tf_opt; % time elapsed
    vel_aver_opt = dist_trav_opt/time_elaps_opt; 
    % assert_v_tg should be 0
    assert_v_tg = abs(vel_aver_opt-v_tgt);
    if assert_v_tg > 1*10^(-tol_ipopt)
        error('Issue when reconstructing average speed')
    end 

    %% Decompose optimal cost
    J_opt           = 0;
    E_cost          = 0;
    A_cost          = 0;
    Arm_cost        = 0;
    Qdotdot_cost    = 0;
    Pass_cost       = 0;
    PassFor_cost    = 0;
    GRF_cost        = 0;
    vA_cost         = 0;
    dFTtilde_cost   = 0;
    QdotdotArm_cost = 0;
    Qskneesec_cost  = 0;
    minBothKCF_cost = 0;
    mindiffKCFGRF_cost=0;
    count           = 1;
    h_opt           = tf_opt/N;
    for k=1:N     
        for j=1:d 
            % Get muscle-tendon lengths, velocities, moment arms
            % Left leg
            qin_l_opt_all = [Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_flex.l*2-1,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_add.l*2-1,1), ...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_rot.l*2-1,1), 0,0,...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_flex.l*2-1,1),0,0.042,0 ...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.ankle.l*2-1,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.subt.l*2-1,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ext*2-1,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ben*2-1,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.rot*2-1,1)];  
            qdotin_l_opt_all = [Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_flex.l*2,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_add.l*2,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_rot.l*2,1),0,0,...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_flex.l*2,1),0,0,0,...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.ankle.l*2,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.subt.l*2,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ext*2,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ben*2,1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.rot*2,1)];  
            [lMTk_l_opt_all,vMTk_l_opt_all,~] = ...
                f_lMT_vMT_dM(qin_l_opt_all,qdotin_l_opt_all);    
            % Right leg
            qin_r_opt_all = [Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_add.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_rot.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_flex.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_tx.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_ty.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_tz.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.subt.r*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ext*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ben*2-1),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.rot*2-1)];  
            qdotin_r_opt_all = [Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_flex.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_add.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.hip_rot.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_add.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_rot.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_flex.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_tx.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_ty.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.knee_tz.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.ankle.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.subt.r*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ext*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.ben*2),...
                Xk_Qs_Qdots_col_opt((k-1)*d+j,jointi.trunk.rot*2)];      
            [lMTk_r_opt_all,vMTk_r_opt_all,~] = ...
                f_lMT_vMT_dM(qin_r_opt_all,qdotin_r_opt_all);
            % Both legs
            lMTk_lr_opt_all     = ...
                [lMTk_l_opt_all([1:44,48:50],1);lMTk_r_opt_all(1:47,1)];
            vMTk_lr_opt_all     = ...
                [vMTk_l_opt_all([1:44,48:50],1);vMTk_r_opt_all(1:47,1)]; 
       
            % Metabolic energy rate
            [~,~,Fce_opt_all,Fpass_opt_all,Fiso_opt_all,vMmax_opt_all,...
                massM_opt_all] = f_forceEquilibrium_FtildeState_all_tendon(...
                    a_col_opt_unsc(count,:)',FTtilde_col_opt_unsc(count,:)',...
                    dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
                    full(vMTk_lr_opt_all),tensions,aTendon,shift);                  
            [~,lMtilde_opt_all] = f_FiberLength_TendonForce_tendon(...
                FTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),aTendon,...
                shift);                
            [vM_opt_all,~] = ...
                f_FiberVelocity_TendonForce_tendon(...
                FTtilde_col_opt_unsc(count,:)',...
                dFTtilde_col_opt_unsc(count,:)',full(lMTk_lr_opt_all),...
                full(vMTk_lr_opt_all),aTendon,shift);   

            if mE == 0 % Bhargava et al. (2004)
                [e_tot_all,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    full(vM_opt_all),full(Fce_opt_all)',full(Fpass_opt_all)',...
                    full(massM_opt_all)',pctsts,full(Fiso_opt_all)',...
                    MTparameters_m(1,:)',body_mass,10); 
                %to debug
                [e_tot_all_debug,~,~,~,~,~] = ...
                getMetabolicEnergySmooth2004all(a_col_opt_unsc(count,:),...
                a_col_opt_unsc(count,:),full(lMtilde_opt_all'),...
                full(vM_opt_all'),full(Fce_opt_all)',full(Fpass_opt_all)',...
                    full(massM_opt_all)',pctsts',full(Fiso_opt_all)',...
                    MTparameters_m(1,:),body_mass,10);
            elseif mE == 1 % Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2003all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);                
            elseif mE == 2 % Umberger (2010)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);            
            elseif mE == 3 % Uchida et al. (2016)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2016all(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);      
            elseif mE == 4 % Umberger (2010) treating muscle lengthening 
                % heat rate as Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all_hl(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);  
            elseif mE == 5 % Umberger (2010) treating negative mechanical 
                % work as Umberger et al. (2003)
                % vMtilde defined for this model as vM/lMopt
                vMtildeUmbk_opt_all = full(vM_opt_all)./(MTparameters_m(2,:)');
                [e_tot_all,~,~,~,~] = fgetMetabolicEnergySmooth2010all_neg(...
                    a_col_opt_unsc(count,:)',a_col_opt_unsc(count,:)',...
                    full(lMtilde_opt_all),...
                    vMtildeUmbk_opt_all,full(vM_opt_all),full(Fce_opt_all)',...
                    full(massM_opt_all)',pctsts,full(vMmax_opt_all)',...
                    full(Fiso_opt_all)',body_mass,10);  
            end
            e_tot_opt_all((k-1)*d+j,:) = full(e_tot_all)';
            e_tot_opt_all_debug((k-1)*d+j,:) = full(e_tot_all_debug)';
                            
            J_opt = J_opt + (1/dist_trav_opt)*(...
                W.E*B(j+1)      *...
                    (f_J94(e_tot_opt_all((k-1)*d+j,:)))/body_mass*h_opt + ...                   
                W.a*B(j+1)      *(f_J94(a_col_opt(count,:)))*h_opt +...      
                W.ArmE*B(j+1)   *(f_J8(e_a_opt(k,:)))*h_opt +...                    
                W.Ak*B(j+1)     *...
                    (f_J26(qdotdot_col_opt(count,residuals_noarmsi)))*h_opt +...                     
                W.passMom*B(j+1)*(f_J17(Tau_pass_opt_all([1:10 14:end],count)))*h_opt + ... 
                W.passFor*B(j+1)*(f_J3(Tau_pass_opt_all([11:13],count)))*h_opt + ... 
                W.u*B(j+1)      *(f_J94(vA_opt(k,:)))*h_opt + ...        
                W.u*B(j+1)      *(f_J94(dFTtilde_col_opt(count,:)))*h_opt + ...                           
                W.AkArms*B(j+1) *(f_J8(qdotdot_col_opt(count,armsi)))*h_opt+...
                W.Qskneesec*B(j+1)*(f_J4(q_col_opt(count,[jointi.knee_add.r:jointi.knee_tx.r jointi.knee_tz.r])))*h_opt+...
                W.minBothKCF*B(j+1)*(sum((KCF_opt_unsc(count,:)/(body_weight*5)).^2))*h_opt+...
                W.diffKCFGRF *B(j+1)*(((GRF_opt_r(count,2)-(sum(KCF_opt_unsc(count,:))))/500).^2)*h_opt);


                E_cost = E_cost + W.E*B(j+1)*...
                    (f_J94(e_tot_opt_all((k-1)*d+j,:)))/body_mass*h_opt;
                A_cost = A_cost + W.a*B(j+1)*...
                    (f_J94(a_col_opt(count,:)))*h_opt;
                Arm_cost = Arm_cost + W.ArmE*B(j+1)*...
                    (f_J8(e_a_opt(k,:)))*h_opt;
                Qdotdot_cost = Qdotdot_cost + W.Ak*B(j+1)*...
                    (f_J26(qdotdot_col_opt(count,residuals_noarmsi)))*h_opt;
                Pass_cost = Pass_cost + W.passMom*B(j+1)*...
                    (f_J17(Tau_pass_opt_all([1:10 14:end],count)))*h_opt;
                PassFor_cost = PassFor_cost + W.passFor*B(j+1)*...
                    (f_J3(Tau_pass_opt_all([11:13],count)))*h_opt;
                vA_cost = vA_cost + W.u*B(j+1)*...
                    (f_J94(vA_opt(k,:)))*h_opt;
                dFTtilde_cost = dFTtilde_cost + W.u*B(j+1)*...
                    (f_J94(dFTtilde_col_opt(count,:)))*h_opt;
                QdotdotArm_cost = QdotdotArm_cost + W.AkArms*B(j+1)*...
                    (f_J8(qdotdot_col_opt(count,armsi)))*h_opt;   
                Qskneesec_cost = Qskneesec_cost + ...
                    W.Qskneesec*B(j+1)*(f_J4(q_col_opt(...
                    count,[jointi.knee_add.r:jointi.knee_tx.r jointi.knee_tz.r])))*h_opt;
                minBothKCF_cost= minBothKCF_cost + ...
                    W.minBothKCF*B(j+1)*(sum((KCF_opt_unsc(count,:)/(body_weight*5)).^2))*h_opt;
                mindiffKCFGRF_cost = mindiffKCFGRF_cost + ...
                    W.diffKCFGRF *B(j+1)*(((GRF_opt_r(count,2)-(sum(KCF_opt_unsc(count,:))))/500).^2)*h_opt;
                count = count + 1;  

                %to debug
                lMTkj_lr_opt_all((k-1)*d+j,:)=full(lMTk_lr_opt_all);
                lMtildej_opt_all((k-1)*d+j,:)=full(lMtilde_opt_all);
        end
    end   

    if Options.minPressuredev||Options.minPressures
        q_knee=q_col_opt_unsc.rad(:,[jointi.knee_flex.r jointi.knee_add.r jointi.knee_rot.r jointi.knee_tx.r jointi.knee_ty.r jointi.knee_tz.r]);
        for i=1:size(q_knee,1)
            out_opt(i,:)=full(F2_debug(q_knee(i,:)));
        end
        pressures_1_opt=full(out_opt(:,9:8+Options.nfacestib1));
        pressures_2_opt=full(out_opt(:,9+Options.nfacestib1:8+Options.nfacestib1+Options.nfacestib2));

        if Options.minPressuredev
            M1_opt=(sum((pressures_1_opt(:)/1e6).^20))^(1/20);
            M2_opt=(sum((pressures_2_opt(:)/1e6).^20))^(1/20);
            J_opt = J_opt + (1/dist_trav_opt)* W.minpressuredev * ((log10(M1_opt/(M2_opt+1e-10))-log10(0.5))^2);
            
            Jminpressdev_cost = W.minpressuredev * ((log10(M1_opt/(M2_opt+1e-10))-log10(0.5))^2);
            Jminpress_cost=0;
        elseif Options.minPressures
            M_opt=(sum(([pressures_1_opt(:);pressures_2_opt(:)]/1e6).^20))^(1/20);
            J_opt = J_opt + (1/dist_trav_opt)*(W.minpressures * M_opt);
    
            Jminpress_cost=W.minpressures * M_opt;
            Jminpressdev_cost=0;
        end
    elseif ~Options.minPressuredev&&~Options.minPressures
        Jminpressdev_cost=0;
        Jminpress_cost=0;
    end

    J_optf = full(J_opt);     
    E_costf = full(E_cost);
    A_costf = full(A_cost);
    Qdotdot_costf = full(Qdotdot_cost);
    Pass_costf = full(Pass_cost);
    PassFor_costf = full(PassFor_cost);
    vA_costf = full(vA_cost);
    dFTtilde_costf = full(dFTtilde_cost);
    Arm_costf = full(Arm_cost);
    QdotdotArm_costf = full(QdotdotArm_cost);
    Qskneesec_costf = full(Qskneesec_cost);
    minBothKCF_costf = full(minBothKCF_cost);
    mindiffKCFGRF_costf=full(mindiffKCFGRF_cost);
    minpressdev_costf=full(Jminpressdev_cost);
    minpressures_costf=full(Jminpress_cost);
    assertCost = J_optf - 1/(dist_trav_opt)*(E_costf+A_costf+Arm_costf+...
        Qdotdot_costf+Pass_costf+PassFor_costf+vA_costf+dFTtilde_costf+...
        QdotdotArm_costf+Qskneesec_costf+minBothKCF_costf+mindiffKCFGRF_costf+Jminpressdev_cost+Jminpress_cost);
    bar(1/(dist_trav_opt)*[E_costf A_costf Arm_costf ...
        Qdotdot_costf Pass_costf PassFor_costf vA_costf dFTtilde_costf ...
        QdotdotArm_costf Qskneesec_costf minBothKCF_costf mindiffKCFGRF_costf Jminpressdev_cost Jminpress_cost]);
    set(gca,'XTickLabels',{'E_costf' 'A_costf' 'Arm_costf' ...
        'Qdotdot_costf' 'Pass_costf' 'PassFor_costf' 'vA_costf' 'dFTtilde_costf' ...
        'QdotdotArm_costf' 'Qskneesec_costf' 'minBothKCF_costf' 'mindiffKCFGRF_costf' 'minpressdev' 'minpress'})
    assertCost2 = abs(stats.iterations.obj(end) - J_optf);
    if assertCost > 1*10^(-tol_ipopt)
        error('Issue when reconstructing optimal cost wrt sum of terms')
    end 
    if assertCost2 > 1*10^(-tol_ipopt)
        error('Issue when reconstructing optimal cost wrt stats')
    end
    X_aux = opti.x;
    f_cost=Function('f_cost',{X_aux},{Jall,Jall_sc, JEall,Jaall,Jakall,JpassMomall,JpassForall,JvAkall,JdFTtildeall,JArmEall,JAjarmsall,JQskneesecall,JminbothKCFall,JdiffKCFGRFall,Jminpressdev,Jminpress,pressures_1,pressures_2,M1,M2,M,dist_trav_tot});
    [Jall_val,Jall_sc_val,JEall_val,Jaall_val,Jakall_val,JpassMomall_val,JpassForall_val,JvAkall_val,JdFTtildeall_val,JArmEall_val,JAjarmsall_val,JQskneesecall_val,JminbothKCFall_val,JdiffKCFGRFall_val, Jminpressdev_val,Jminpress_val,pressures_1_val,pressures_2_val,M1_val,M2_val,M_val,dist_trav_tot_val]=f_cost(w_opt);
    f_lMT_debug=Function('f_lMT_debug',{X_aux},{lMtildej_all_debug,lMTj_lr_all_debug});
    [lMtildej_all_debug_val,lMTj_lr_all_debug_val]=f_lMT_debug(w_opt);
    f_coll_eq_constr=Function('f_coll_eq_constr',{X_aux},{coll_eq_constr});
    [coll_eq_constr_val]=f_coll_eq_constr(w_opt);
    f_ineq_constr=Function('f_ineq_constr',{X_aux},{coll_ineq_constr1,coll_ineq_constr2,coll_ineq_constr3,coll_ineq_constr4,coll_ineq_constr5});
    [coll_ineq_constr1_val,coll_ineq_constr2_val,coll_ineq_constr3_val,coll_ineq_constr4_val,coll_ineq_constr5_val]=f_ineq_constr(w_opt);
    f_Tj_Tau=Function('f_Tj_Tau',{X_aux},{Tj_all_debug,Tau_passj_all_d_debug});
    [Tj_all_debug_val,Tau_passj_all_d_debug_val]=f_Tj_Tau(w_opt);
    f_FT_MA=Function('f_FT_MA',{X_aux},{FTj_all_debug,MAj_knee_ty_r_all_debug});
    [FTj_all_debug_val,MAj_knee_ty_r_all_debug_val]=f_FT_MA(w_opt);
    f_Qs_Qdots_j=Function('f_Qs_Qdots_j',{X_aux},{Qsj_nsc_debug,Qdotsj_nsc_debug});
    [Qsj_nsc_debug_val,Qdotsj_nsc_debug_val]=f_Qs_Qdots_j(w_opt);
    X_Qs_Qdots_j_debug(1:2:34*2,:)=full(Qsj_nsc_debug_val);
    X_Qs_Qdots_j_debug(2:2:34*2,:)=full(Qdotsj_nsc_debug_val);

    %% Reconstruct full gait cycle

    % joint accelerations controls on mesh points (2:N)
    qddot_opt_unsc.deg = qdotdot_col_opt_unsc.deg(d:d:end,:);
    qddot_opt_unsc.rad = qdotdot_col_opt_unsc.rad(d:d:end,:);

    % express slack controls on mesh points 1:N to be consistent
    qddot_opt_unsc.deg = [qddot_opt_unsc.deg(end,:); qddot_opt_unsc.deg(1:end-1,:)];
    qddot_opt_unsc.rad = [qddot_opt_unsc.rad(end,:); qddot_opt_unsc.rad(1:end-1,:)];
    dFTtilde_opt_unsc = [dFTtilde_col_opt_unsc(end,:); dFTtilde_col_opt_unsc(1:end-1,:)];

    %% Gait cycle starts at right side initial contact

    % % Ground reaction forces at mesh points (1:N-1)
    % Foutk_opt                   = zeros(size(q_opt_unsc.rad,1),F2_debug.nnz_out);
    % for i = 1:size(q_opt_unsc.rad,1)
    %     % Create zero input vector for external function
    %     F_ext_input = zeros(F2_debug.nnz_in,1);
    %     % Assign Qs, Qdots and Qdotdots
    %     F_ext_input=[Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)'];
    % 
    %     % Evaluate external function
    %     res = F2_debug(F_ext_input);
    %     Foutk_opt(i,:) = full(res);
    % end
    % GRFk_opt = Foutk_opt(:,[35:37 38:40]); %first right and then left GRFs
    GRFk_opt=[GRF_opt_r_k GRF_opt_l_k]; 
    [idx_GC,idx_GC_base_forward_offset,HS1,HS_threshold] = getStancePhaseSimulation(GRFk_opt,body_mass/3);
    GRFj_opt=[GRF_opt_r GRF_opt_l];
    [idx_GCj,idx_GC_base_forward_offsetj,HS1j,HS_thresholdj] = getStancePhaseSimulation(GRFj_opt,body_mass/3);

    Qs_GC = q_opt_unsc.deg(idx_GC,:);
    Qs_GC_j=q_col_opt_unsc.deg(idx_GCj,:);
    Qdots_GC = qdot_opt_unsc.deg(idx_GC,:);
    Qdots_GC_j =qdot_col_opt_unsc.deg(idx_GCj,:);
    Qdotdots_GC = qddot_opt_unsc.deg(idx_GC,:);
    Acts_GC = a_opt_unsc(idx_GC,:);
    Acts_GCj= a_col_opt_unsc(idx_GCj,:);
    dActs_GC = vA_opt_unsc(idx_GC,:);
    FTtilde_GC = FTtilde_opt_unsc(idx_GC,:);
    FTtilde_GCj= FTtilde_col_opt_unsc(idx_GCj,:);
    dFTtilde_GC = dFTtilde_opt_unsc(idx_GCj,:);
    a_a_GC = a_a_opt_unsc(idx_GC,:);
    e_a_GC = e_a_opt_unsc(idx_GC,:);
    KCF_GCj = KCF_opt_unsc(idx_GCj,:);
    GRF_GC = GRFk_opt(idx_GC,:);
    GRF_GCj = GRFj_opt(idx_GCj,:);
    if Options.minPressuredev||Options.minPressures
        pressures_1_GCj=pressures_1_opt(idx_GCj,:);
        pressures_2_GCj=pressures_2_opt(idx_GCj,:);
    else
        pressures_1_GCj=pressures_opt_1(idx_GCj,:);
        pressures_2_GCj=pressures_opt_2(idx_GCj,:);
    end
    MuscleMoments_GCj=MuscleMoments(idx_GCj,:);
    Tauk_out_GCj= Tauk_out(idx_GCj,:);

    % adjust forward position to be continuous and start at 0
    Qs_GC(idx_GC_base_forward_offset,jointi.pelvis.tx) = Qs_GC(idx_GC_base_forward_offset,jointi.pelvis.tx) + dist_trav_opt;
    Qs_GC(:,jointi.pelvis.tx) = Qs_GC(:,jointi.pelvis.tx) - Qs_GC(1,jointi.pelvis.tx);
    
    %% Unscale actuator torques
    T_a_GC = zeros(size(a_a_GC));
    for i=1:size(a_a_GC,2)
        T_a_GC(:,i) = a_a_GC(:,i).*scaling.ArmTau;
    end

    Options.W=W;

    %% Save Results
    Results.Qs_opt=Qs_GC;
    Results.Qs_optj=Qs_GC_j;
    Results.Qdots_opt=Qdots_GC;
    Results.Qdots_optj=Qdots_GC_j;
    Results.GRFs_opt=GRF_GC;
    Results.GRFs_optj=GRF_GCj;
    Results.objective.E_cost=E_costf;
    Results.objective.A_cost=A_costf;
    Results.objective.Arm_cost=Arm_costf;
    Results.objective.Qdotdot_cost=Qdotdot_costf;
    Results.objective.Pass_cost=Pass_costf;
    Results.objective.vA_cost=vA_costf;
    Results.objective.dFTtilde_cost=dFTtilde_costf;
    Results.objective.QdotdotArm_cost=QdotdotArm_costf;
    Results.objective.Qskneesec_cost=Qskneesec_costf;
    Results.objective.total=J_opt;
    Results.tgrid=tgrid;
    Results.tgrid_ext=tgrid_ext;
    Results.muscle_names=muscleNames;
    Results.muscles.a=Acts_GC;
    Results.muscles.a_j=Acts_GCj;
    Results.muscles.da=dActs_GC;
    Results.muscles.FTtilde=FTtilde_GC;
    Results.muscles.FTtilde_j=FTtilde_GCj;
    Results.muscles.dFTtilde=dFTtilde_GC;
    Results.muscles.MuscleMoments=MuscleMoments_GCj;
    Results.torque_actuators.a=a_a_GC;
    Results.torque_actuators.e=e_a_GC;
    Results.torque_actuators.T=T_a_GC;
    Results.GRFs.threshold=HS_threshold;
    Results.GRFs.initial_contact_side=HS1;
    Results.GRFs.idx_GC=idx_GC;
    Results.KCFj=KCF_GCj;
    Results.spatiotemp.dist_trav=dist_trav_opt;
    Results.Options=Options;
    Results.W=W;
    if Options.minPressuredev
        Results.M1 = M1_opt;
        Results.M2 = M2_opt;
    elseif Options.minPressures
        Results.M=M_opt;
    end
    Results.pressures.compartment1 = pressures_1_GCj;
    Results.pressures.compartment2 = pressures_2_GCj;

    Results.Tauk_out_GCj=Tauk_out_GCj;
    Results.stride_r=dist_r;
    Results.stride_l=dist_l;
    Results.stride_width=StepWidth_opt;
    Results.calcn_r=calcn_opt_r_k;
    Results.calcn_l=calcn_opt_l_k;

    %save results
Outname = ['..\Results\PredSim_3D_GC_v2\Results_PredSim_a' num2str(W.a) ...
           '_E' strrep(num2str(W.E),'.','p') ...
           '_P' strrep(num2str(W.minpressuredev),'.','p') '.mat'];
disp(['Saving results as: ' Outname]);
save(Outname,'w_opt','stats','setup','Results','Options');
            

    keyboard;

    % Two gait cycles
    t_mesh = [tgrid(1:N) tgrid(N+1)+tgrid(1:N)-tgrid(1)] ;
    % Joint angles
    q_opt_GUI_GC_1 = Qs_GC;
    q_opt_GUI_GC_2 = q_opt_GUI_GC_1;
    q_opt_GUI_GC_2(:,jointi.pelvis.tx) =...
        q_opt_GUI_GC_2(:,jointi.pelvis.tx) +...
        dist_trav_opt;
    JointAngle.labels = {'time','pelvis_rz','pelvis_rx','pelvis_ry','pelvis_tx',...
                'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
                'hip_rotation_l','hip_flexion','hip_adduction','hip_rotation',...
                'knee_flexion_l','knee_flexion','knee_adduction','knee_rotation',...
                'knee_tx','knee_ty','knee_tz','ankle_angle_l','ankle_angle',...
                'subtalar_angle_l','subtalar_angle',...
                'lumbar_extension','lumbar_bending','lumbar_rotation',...
                'arm_flex_l','arm_add_l','arm_rot_l',...
                'arm_flex_r','arm_add_r','arm_rot_r',...
                'elbow_flex_l','elbow_flex_r','pro_sup_l','pro_sup_r'};
    
    q_opt_GUI_GC = [t_mesh',[q_opt_GUI_GC_1;q_opt_GUI_GC_2]];
    % Muscle activations (to have muscles turning red when activated).
    Acts_GC_GUI = [Acts_GC;Acts_GC];
    % Combine data joint angles and muscle activations
    JointAngleMuscleAct.data = [q_opt_GUI_GC,Acts_GC_GUI];
    % Combine labels joint angles and muscle activations
    JointAngleMuscleAct.labels = JointAngle.labels;
    for i = 1:NMuscle/2
            Muscle_names_all{i} = ...
                [muscleNames{i}(1:end-2),'_l'];
            Muscle_names_all{i+NMuscle/2} = ...
                [muscleNames{i}(1:end-2),'_r'];
    end
    for i = 1:NMuscle
        JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC,2)} = ...
            [Muscle_names_all{i},'/activation'];
    end
    JointAngleMuscleAct.inDeg = 'yes';
    filenameJointAngles = fullfile('..\Results\PredSim_3D_GC_v2\',...
        ['out_motion.mot']);
    write_motionFile(JointAngleMuscleAct, filenameJointAngles);

     %% Recalculate constraint values for joint moments
    for i=1:N
        for j=1:d
            


            %% Recalculate dynamic constraints
            % Skeleton dynamics (implicit formulation)      
            
            
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
            %         Qsp_nsc          = Xkj_nsc(1:2:end,:)*C(:,j+1);
            %         Qdotsp_nsc       = Xkj_nsc(2:2:end,:)*C(:,j+1);    
            %         % Skeleton dynamics (implicit formulation)               
            %         qdotj_nsc = Xkj_nsc(2:2:end,j+1); % velocity
            %         eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.QsQdots(1:2:end)';
            %         eq_constr{end+1} = (h*Ak_nsc - Qdotsp_nsc)./...
            %                scaling.QsQdots(2:2:end)';
                    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
            X_col_opt_nsc = Xk_Qs_Qdots_col_opt;
            
            q_colkj=[q_opt_unsc.rad(i,:); X_col_opt_nsc((i-1)*d+1:i*d,1:2:end)]; %coordinates for all 4 collocation points of interval i unscaled
            qdot_colkj =[qdot_opt_unsc.rad(i,:); X_col_opt_nsc((i-1)*d+1:i*d,2:2:end)]; % velocity
    
            Qsp_opt_nsc = q_colkj'*C(:,j+1);
            Qdotsp_opt_nsc = qdot_colkj'*C(:,j+1);
            
            FTtildekj_opt_unsc=[FTtilde_opt_unsc(i,:); FTtilde_col_opt_unsc(((i-1)*d+1):((i-1)*d+3),:)];
            FTtildep_opt_unsc=FTtildekj_opt_unsc'*C(:,j+1);
    
            akj_opt_unsc=[a_opt_unsc(i,:); a_col_opt(((i-1)*d+1):((i-1)*d+3),:)];
            ap_opt_unsc=akj_opt_unsc'*C(:,j+1);
    
            eq_constr_dynkin((i-1)*d+j,:)= (h_opt*qdot_colkj(j+1,:) - Qsp_opt_nsc')./scaling.QsQdots(1:2:end);
            
            eq_constr_dynkinvel((i-1)*d+j,:)= (h_opt*qdotdot_col_opt_unsc.rad((i-1)*d+j,:) - Qdotsp_opt_nsc')./scaling.QsQdots(2:2:end);
           
            eq_constr_dynFT((i-1)*d+j,:)=(h_opt*dFTtilde_col_opt_unsc((i-1)*d+j,:)-FTtildep_opt_unsc')./scaling.FTtilde;
    
            eq_constr_act((i-1)*d+j,:)=(h_opt*vA_opt_unsc(i,:)-ap_opt_unsc')./scaling.a;

            %to debug
            MA_knee_ty_r_opt((i-1)*d+j,:)=full(MA.knee_ty.r);

            %ineq activation dynamics
            act1_opt(:,(i-1)*d+j) = vA_opt_unsc(i,:) + akj_opt_unsc(j+1,:)./(ones(1,size(akj(:,j+1),1))*tdeact); %needs to be >0
            act2_opt(:,(i-1)*d+j) = vA_opt_unsc(1,:) + akj_opt_unsc(j+1,:)./(ones(1,size(akj(:,j+1),1))*tact); %needs to be < 1/tact

        end

        %% check continuity constraints
         % Variables within current mesh interval
        % States      
        akj_opt = [a_opt_unsc(i,:); a_col_opt_unsc((i-1)*d+1:i*d,:)]; 
        FTtildekj_opt = [FTtilde_opt_unsc_all(i,:); FTtilde_col_opt_unsc((i-1)*d+1:i*d,:)];
        Qskj_nsc_opt = [q_opt_unsc_all.rad(i,:); q_col_opt_unsc.rad((i-1)*d+1:i*d,:)];
        Qdotskj_nsc_opt = [qdot_opt_unsc_all.rad(i,:); qdot_col_opt_unsc.rad((i-1)*d+1:i*d,:)];
        a_akj_opt = [a_a_opt_unsc_all(i,:); a_a_col_opt_unsc((i-1)*d+1:i*d,:)];
        % Add equality constraints (next interval starts with end values of 
        % states from previous interval)
        eq_constr_cont_a(i,:)=(a_opt_unsc_all(i+1,:)' - akj_opt'*D);
        eq_constr_cont_FTtilde(i,:)=(FTtilde_opt_unsc_all(i+1,:)' - FTtildekj_opt'*D); % scaled
        eq_constr_cont_q(i,:)=(q_opt_unsc_all.rad(i+1,:)' - Qskj_nsc_opt'*D); % scaled
        eq_constr_cont_qdot(i,:)= qdot_opt_unsc_all.rad(i+1,:)' - Qdotskj_nsc_opt'*D; % scaled
        eq_constr_cont_a_a(i,:)=a_a_opt_unsc_all(i+1,:)' - a_akj_opt'*D;
    end
    %% Check peridoicity of the states
    eq_constr_per_qs=Qs_opt(end,[jointi.pelvis.tilt:jointi.pelvis.rot jointi.pelvis.ty:end])-Qs_opt(1,[jointi.pelvis.tilt:jointi.pelvis.rot jointi.pelvis.ty:end]); 
    eq_constr_per_qdots=Qdots_opt(end,:)-Qdots_opt(1,:);
    eq_constr_per_a=a_opt(end,:)-a_opt(1,:);
    % Muscle-tendon forces
    eq_constr_per_FTtilde=FTtilde_opt(end,:)-FTtilde_opt(1,:);
    % Torque actuator activations
    eq_constr_per_a_a=a_a_opt(end,:)-a_a_opt(1,:);
    
    if Options.defineGRFphases
        %GRF stance and swing
        ineq_constr_GRFr_stance=GRF_opt_r_k(Nstance_r,2);
        ineq_constr_GRFr_swing=GRF_opt_r_k(Nswing_r,2);
        ineq_constr_GRFl_stance=GRF_opt_l_k(Nstance_l,2);
        ineq_constr_GRFl_swing=GRF_opt_l_k(Nswing_l,2);
    end

    

end

function [Tauk_out, KCF_opt_unsc, GRF_opt_r, GRF_opt_l, calcn_opt_r, ...
    calcn_opt_l,pressures_opt_1,pressures_opt_2]=ExtractOptDataFromExternalFuncs(Xk_Qs_Qdots_opt,Xk_Qdotdots_opt,F2_in,deri,residualsi,KCFi,armsi,scaling,tol_ipopt,a_a_opt_unsc,nfacestib1,nfacestib2,Options,jointi)

    if Options.KCFasinputstoExternalFunction
        F2_skeletal_debug=F2_in{1};
        F2=F2_in{2};
    else
        F2=F2_in;
    end  

    if Options.KCFasinputstoExternalFunction
        q_knee=Xk_Qs_Qdots_opt([[jointi.knee_flex.r jointi.knee_add.r jointi.knee_rot.r jointi.knee_tx.r jointi.knee_ty.r jointi.knee_tz.r]*2-1]);
        out=F2(q_knee);
                SumForces=full(out(3:5));
                SumMoments=full(out(6:8));
                out2_res = ...
                        F2_skeletal_debug([Xk_Qs_Qdots_opt;Xk_Qdotdots_opt;SumForces;SumMoments]); 
    else
        if deri == 2
            out2_res = F2(Xk_Qs_Qdots_opt,Xk_Qdotdots_opt); 
        else
            out2_res = ...
                F2([Xk_Qs_Qdots_opt;Xk_Qdotdots_opt]);     
        end
    end
        out2_res_opt(:,:) = full(out2_res); 

        % Optimal joint torques, ground reaction forces and moments
        Tauk_out        = out2_res_opt(residualsi);
        if Options.KCFasinputstoExternalFunction
            KCF_opt_unsc  = full(out(1:2));
            GRF_opt_r=out2_res_opt(35:37); 
            GRF_opt_l=out2_res_opt(38:40);
            calcn_opt_r=out2_res_opt(41:43);
            calcn_opt_l=out2_res_opt(44:46);
            pressures_opt_1=full(out(9:8+nfacestib1));
            pressures_opt_2=full(out(9+nfacestib1:8+nfacestib1+nfacestib2));
        else
            keyboard; %%to check!
            GRF_opt_r=out2_res_opt(35:37); 
            GRF_opt_l=out2_res_opt(38:40);
            calcn_opt_r=out2_res_opt(41:43);
            calcn_opt_l=out2_res_opt(44:46);
            KCF_opt_unsc  = out2_res_opt(47:48);
            pressures_opt_1=out2_res_opt(49:48+nfacestib1);
            pressures_opt_2=out2_res_opt(49+nfacestib1:48+nfacestib1+nfacestib2);
        end

        % assertArmTmax should be 0
        assertArmTmax = max(max(abs(Tauk_out(armsi)-(a_a_opt_unsc')*...
            scaling.ArmTau))); 
        if assertArmTmax > 1*10^(-tol_ipopt)
            disp('Issue when reconstructing residual forces')
        end 

end