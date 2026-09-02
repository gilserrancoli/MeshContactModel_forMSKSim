solAD=load('C:\Gil\MeshesInAD\gaitWithKneeProsthesis\trackingSimulations_3D\Results\TrackSim_3D_GC\Optsol_withbettertrack_badkneeflex_withJRM_v2_openaccbounds_negslope\Results_3D.mat');
t0=solAD.tgrid_ext(1);
tf=solAD.tgrid_ext(end);
q_opt=solAD.Results_3D.Derivative_AD_Recorder.Hessian_Approximated.LinearSolver_mumps.InitialGuess_2.Qs_opt;
q_opt(:,[1:3 7:16 20:34])=q_opt(:,[1:3 7:16 20:34])*pi/180;

out=MatlabFunction_CalculatePressures(q_opt(1,:))