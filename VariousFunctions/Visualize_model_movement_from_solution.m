% fem=stlread('Femoral Component simpler.stl');
% tib=stlread('Tibial Insert simpler.stl');
close all
clear all
clc

options.debug=0; %if 0, Results_3D.mat is grabbed from the results folder in the Github directory
if options.debug==0
    solpath='..\Results\TrackSim_3D_GC\Results_3D.mat';
elseif options.debug==1
    solpath='C:\Users\jbadi\Documents\TFG\Resultats(postJbugfix)(14.4)\Resultats(postJbugfix)\WKCF=10(noKCFcab)\TrackSim_3D_GC\Results_3D.mat';
elseif options.debug==2
    solpath='C:\Users\jbadi\Documents\TFG\WQs100_WKCF10_WKCFres1_WJRMres1_withJRM_JRMboundsascons\WQs100_WKCF10_WKCFres1_WJRMres1_withJRM_JRMboundsascons\Results_3D.mat';
end
    

fem=stlread('..\..\..\contactsKneeProsthesis\Femoral Component simpler.stl');
tib=stlread('..\..\..\contactsKneeProsthesis\Tibial Insert simpler_superior.stl');

centerfem=[0 0.042 0];
patch('Faces',fem.ConnectivityList,'Vertices',fem.Points+centerfem,'Facecolor',[0.8 0.8 1]);
hold all;
patch('Faces',tib.ConnectivityList,'Vertices',tib.Points,'Facecolor',[0.2 0.6 0.2]);

plot3(0,0,0,'ko','Color',[0.2 0.6 0.2],'MarkerFaceColor','k');
plot3(centerfem(1),centerfem(2),centerfem(3),'ko','Color',[0.2 0.6 0.2],'MarkerFaceColor','k');

IKdata=importdata('..\..\..\contactsKneeProsthesis\IK_ngait_og1_Extmuscle_withfluoro.mot');
knee_rot=IKdata.data(89:220,11:13)*pi/180;
knee_tra=IKdata.data(89:220,14:16);

import casadi.*
%Function for rotation 3x3
psi=MX.sym('psi',1);
theta=MX.sym('theta',1);
phi=MX.sym('phi',1);
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];

R=R1*R2*R3;
fR3x3=Function('fR',{psi,theta,phi},{R});

%function for translation 4x4
x=MX.sym('x',1);
y=MX.sym('y',1);
z=MX.sym('z',1);
Rtrans=[1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1];
ftrans=Function('ftrans',{x,y,z},{Rtrans});

%function for translation 0.042 4x4
Rtrans0042=[1 0 0 0; 0 1 0 0.042; 0 0 1 0; 0 0 0 1];

%function for rotation 4x4
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
R_aux=R1*R2*R3;
Rrot=MX(4,4);
Rrot(1:3,1:3)=R_aux;
Rrot(4,4)=1;
fR=Function('fR',{psi,theta,phi},{Rrot});

% Rtib=Rtrans*(Rrot*inv(Rtrans));
% fTransf=Function('fTransf',{x,y,z,psi,theta,phi},{Rtib}); %error defining translations (wrong origin and coordinate axes)

aux=R*[-x;-y;-z];
Rtranstib_4x4=MX(4,4);
Rtranstib_4x4(1,1)=1;
Rtranstib_4x4(2,2)=1;
Rtranstib_4x4(3,3)=1;
Rtranstib_4x4(4,4)=1;
Rtranstib_4x4(1:3,4)=aux;

Rtib=Rtranstib_4x4*Rtrans0042*Rrot;
fTransf=Function('fTransf_v2',{x,y,z,psi,theta,phi},{Rtib});


%% Write forces to visualize in blender
solAD=load(solpath); %Results_3D.mat file to visualise here
load('..\Blender\tgrid_ext.mat'); %time grid
solAD.tgrid_ext=tgrid_ext;
solAD.tgrid=solAD.Results_3D.NMesh_50.Qs_toTrack;
t0=solAD.tgrid_ext(1);
tf=solAD.tgrid_ext(end);

q_opt=solAD.Results_3D.Derivative_AD_Recorder.Hessian_Approximated.LinearSolver_mumps.InitialGuess_2.Qs_opt;
q_opt(:,[1:16 20:34])=q_opt(:,[1:16 20:34])*pi/180;

fContact=external('f','..\..\..\contactsKneeProsthesis\CalculateKneeContact_k10000_sepML_forTesting.dll');
for i=1:size(q_opt,1)
    [Sum_Force,Sum_Moments,Sum_Force_part1,Sum_Force_part2,Sum_Moments_part1,Sum_Moments_part2,p_all1,p_all2]=fContact(q_opt(i,17),q_opt(i,18),q_opt(i,19),-q_opt(i,14),-q_opt(i,15),q_opt(i,16));
    Sum_Force_full(i,:)=full(Sum_Force);
    Sum_Moments_full(i,:)=full(Sum_Moments);
    Sum_Force_part1_full(i,:)=full(Sum_Force_part1);
    Sum_Force_part2_full(i,:)=full(Sum_Force_part2);
    Sum_Moments_part1_full(i,:)=full(Sum_Moments_part1);
    Sum_Moments_part2_full(i,:)=full(Sum_Moments_part2);
    p_all1_full(i,:)=full(p_all1);
    p_all2_full(i,:)=full(p_all2);
end
csvwrite('pressures1.csv',p_all1_full);
csvwrite('pressures2.csv',p_all2_full);

%% check visualization with movement

% solAD_interp=interp1(solAD.tgrid_ext,solAD. x_opt_ext,solAD.tgrid_ext(1):(tf-t0)/100:solAD.tgrid_ext(end));
IK_exp_interp=interp1(IKdata.data(:,1),IKdata.data,solAD.tgrid);
knee_rot=IK_exp_interp(:,11:13)*pi/180;
knee_tra=IK_exp_interp(:,14:16);

Options.writeVideo=0;
if Options.writeVideo
    v=VideoWriter('video_model_states_from_solution.avi');
    v.FrameRate=20/(tf-t0);
    open(v);
end

fem=stlread('..\..\..\contactsKneeProsthesis\Femoral Component simpler_superior_94.stl');
tib=stlread('..\..\..\contactsKneeProsthesis\Tibial Insert simpler_superior_24.stl');

faces1=csvread('..\..\..\contactsKneeProsthesis\pairs1tibia_forblender.csv');
faces2=csvread('..\..\..\contactsKneeProsthesis\pairs2tibia_forblender.csv');
load('..\..\..\contactsKneeProsthesis\pairs_mesh94x24.mat');
%faces1=pairs(:,1);
%faces2=pairs(:,2);
pressures1=csvread('pressures1.csv');
pressures2=csvread('pressures2.csv');
for i=1:size(q_opt,1)
    patch('Faces',fem.ConnectivityList,'Vertices',fem.Points+centerfem,'Facecolor',[0.8 0.8 1],'FaceAlpha',0.5);
    hold all;
%     Mtransformtib=full(fTransf(knee_tra(i,1),knee_tra(i,2),knee_tra(i,3),knee_rot(i,3),knee_rot(i,1),knee_rot(i,2)));
%     for j=1:size(tib.Points,1)
%         pointsTib(j,:)=Mtransformtib*[tib.Points(j,:) 1]';
%     end
%     hold all;
%     patch('Faces',tib.ConnectivityList,'Vertices',pointsTib(:,1:3),'Facecolor',[0.2 0.6 0.2],'FaceAlpha',0.5);

    %get color of faces
    for j=1:size(tib.ConnectivityList,1)
        if any(j==faces1)
            pressurej=pressures1(i,find(j==faces1));
        elseif any(j==faces2)
            pressurej=pressures2(i,find(j==faces2));
        else
        end
        if pressurej<1.1e6
            pressurej=1.1e6;
        elseif pressurej>5.5e7
            pressurej=5.5e7;
        else
        end
        %C(j,:)=[(pressurej-8.88e9)/(1.41e10-8.88e9) 0 (pressurej-1.41e10)/(8.88e9-1.41e10)];    
        C(j,:)=(pressurej-1.1e6)/(5.5e7-1.1e6);
    end

    %get kinematics transformation
    Mtransformtib2=full(fTransf(q_opt(i,17),q_opt(i,18),q_opt(i,19),-q_opt(i,14),-q_opt(i,15),q_opt(i,16)));
    for j=1:size(tib.Points,1)
        pointsTib2(j,:)=Mtransformtib2*[tib.Points(j,:) 1]';
    end
    hold all;
    Cmap=load('CustomColorMap.mat');
    patch('Faces',tib.ConnectivityList,'Vertices',pointsTib2(:,1:3),'FaceVertexCData',C,'FaceColor','flat');
    colormap(Cmap.CustomColormap);
    set(gca,'View',[1.268e+02,5.372],'CameraPosition',[0.5315,0.3845,0.0764],'XTick',[],'YTick',[],'ZTick',[],'XColor','none','YColor','none','ZColor','none','Color','none');
    xlim([-0.04 0.04]);
    ylim([0 0.0644]);
    zlim([-0.04 0.04]);
%    xlim([-0.08 0.08]);
%    ylim([0 0.0644*2]); %debug purposes, ignore
%    zlim([-0.08 0.08]);
    [leg_hand icons]=legend({'','model data tib comp'},'Location','southeast','box','off');
%     set(icons(1),'rotation',90)
%     title(num2str(i));
%     kin2blender(i,4:6)=Mtransformtib2(1:3,4);
%     kin2blender(i,1:3)=q_opt(i,17:19);

    agtheta = asin(Mtransformtib2(3,2));
    agphi = asin(Mtransformtib2(3,1)/(-cos(agtheta)));
    agpsi = asin(Mtransformtib2(1,2)/(-cos(agtheta)));


    kin2blender(i,4:6)= [agtheta agphi agpsi]; %Blender uses rad for Euler rotations
    kin2blender(i,1:3)=Mtransformtib2(1:3,4);

    if Options.writeVideo
        F=getframe(gcf);
        writeVideo(v,F)
    end
    pause(0.1);
%    if i==43
%        keyboard; %debug purposes, ignore
%    end
   clf;
end
if Options.writeVideo
    close(v);
end

%% Write kinematics data to visualize in Blender
csvwrite('kin2blender.csv',kin2blender);

%% Write pressure map to paint each face in blender
pressuremat=[];
for i=1:length(faces1)
    pressuremat(:,faces1(i))=pressures1(:,i);
end
for i=1:length(faces2)
    pressuremat(:,faces2(i))=pressures2(:,i);
end

writematrix(pressuremat,'pressure2blender.csv');

%% Write stats file for Blender
minptot=min(min(pressuremat));
maxptot=max(max(pressuremat)); %absolute minima and maxima of the sample

 minpinst=[];
 maxpinst=[]; %minima and maxima at each instant

for t=1:length(pressuremat)
    minpinst(t)=min(pressuremat(t,:));
    maxpinst(t)=max(pressuremat(t,:));
end

statsblender=minpinst;
statsblender(2,:)=maxpinst;
statsblender(3,1)=minptot;
statsblender(3,2)=maxptot;

writematrix(statsblender,'statsblender.csv');

%% Adapt kin2blender to particular python program

kin2btransposed=kin2blender';
writematrix(kin2btransposed,'kin2btransposed.csv')