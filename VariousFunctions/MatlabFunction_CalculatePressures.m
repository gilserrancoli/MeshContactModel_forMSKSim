function out=MatlabFunction_CalculatePressures(q)


import casadi.*

%Load input data
% fem=stlread('Femoral Component simpler_superior.stl');
% tib=stlread('Tibial Insert simpler_superior.stl');
load('pairs_mesh94x24.mat');

% centerfem=[0 0.042 0];
Options.sepML=1; %whether to separate into medial and lateral forces (1) or not (0)
Options.optE=0;

h=0.006;
poisson=0.46; %0.45 in ISB 2021
if Options.optE
    E=MX.sym('E',1);
else
%     E=5*10e6; % in ISB 2021
    E=400*1e6; 
end
centerfem=[0 0.042 0];

% %% check visualization
fem=stlread('C:\Gil\MeshesInAD\contactsKneeProsthesis\Femoral Component simpler_superior_94.stl');
tib=stlread('C:\Gil\MeshesInAD\contactsKneeProsthesis\Tibial Insert simpler_superior_24.stl');
% patch('Faces',fem.ConnectivityList,'Vertices',fem.Points+centerfem,'Facecolor',[0.8 0.8 1]);
% hold all;
% patch('Faces',tib.ConnectivityList,'Vertices',tib.Points(:,1:3),'Facecolor',[0.2 0.6 0.2]);

%% Discard too far cells  
femPoints=fem.Points+centerfem;
pointsfem=find(((femPoints(:,2))>0.00375)&((femPoints(:,2))<0.05));
pointstib=find((tib.Points(:,2)>0.00375)&(tib.Points(:,2)<0.04));
nfacesFem=size(fem.ConnectivityList,1);
nfacesTib=size(tib.ConnectivityList,1);

[m,Ifem]=sort(fem.ConnectivityList(:,1));
facesFem=fem.ConnectivityList(Ifem,:);
[m,Itib]=sort(tib.ConnectivityList(:,1));
facesTib=tib.ConnectivityList(Itib,:);

facesFem_1=[];
facesTib_1=[];
k=1;
for i=1:nfacesFem
    j=1;
    found=false;
    while ~found & (j<size(pointsfem,1))
        if any(facesFem(i,:)==pointsfem(j))
            facesFem_1(k,:)=facesFem(i,:);
            k=k+1;
            found=true;
        end
        j=j+1;
    end
end
k=1;
for i=1:nfacesTib
    j=1;
    found=false;
    while ~found & (j<size(pointstib,1))
        if any(facesTib(i,:)==pointstib(j))
            facesTib_1(k,:)=facesTib(i,:);
            k=k+1;
            found=true;
        end
        j=j+1;
    end
end


%% Define casadi functions
%Function for rotation 3x3
psi=q(14); %q needs to be in radians
theta=q(15);
phi=q(16);
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];

R=R1*R2*R3;
% fR3x3=Function('fR',{psi,theta,phi},{R});

%function for translation 4x4
x=q(17);
y=q(18);
z=q(19);
Rtrans=[1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1];
% ftrans=Function('ftrans',{x,y,z},{Rtrans});

%function for translation 0.042 4x4
Rtrans0042=[1 0 0 0; 0 1 0 0.042; 0 0 1 0; 0 0 0 1];

%function for translation and change of orientation 4x4
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
R_aux=R1*R2*R3;
% Rrottrans=MX(4,4);
Rrottrans(1:3,1:3)=R_aux;
Rrottrans(4,4)=1;
Rrottrans(2,4)=0;
% frottrans=Function('frottrans',{psi,theta,phi},{Rrottrans});

%function for rotation 4x4
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
R_aux=R1*R2*R3;
% Rrot=MX(4,4);
Rrot(1:3,1:3)=R_aux;
Rrot(4,4)=1;
% fR=Function('fR',{psi,theta,phi},{Rrot}); %4x4 matrix
% f_Raux=Function('f_Raux',{x,y,z,psi,theta,phi},{R_aux}); %only the 3x3 matrix

%translation in tibia coord system
aux=R*[-x;-y;-z];
% Rtranstib_4x4=MX(4,4);
Rtranstib_4x4(1,1)=1;
Rtranstib_4x4(2,2)=1;
Rtranstib_4x4(3,3)=1;
Rtranstib_4x4(4,4)=1;
Rtranstib_4x4(1:3,4)=aux;

Rtib=Rtranstib_4x4*Rtrans0042*Rrot;
% fTransf=Function('fTransf',{x,y,z,psi,theta,phi},{Rtib});


if Options.sepML==1
    facesTib_1_part1=[];
    facesTib_1_part2=[];
    tib_part1.ConnectivityList=[];
    tib_part2.ConnectivityList=[];
    tib_part1.Points=tib.Points;
    tib_part2.Points=tib.Points;
    pairs1=[];
    pairs2=[];
    i1=1;
    i2=1;
    faces_1_part1_in_global=[];
    faces_1_part2_in_global=[];
    %Separate into Part1 Lateral, and Part2 Medial
    for i=1:size(facesTib_1,1)
            if (tib.Points(facesTib_1(i,1),3)>0)&(tib.Points(facesTib_1(i,2),3)>0)&(tib.Points(facesTib_1(i,3),3)>0)
                facesTib_1_part1=[facesTib_1_part1; facesTib_1(i,:)];
                tib_part1.ConnectivityList=[tib_part1.ConnectivityList; tib.ConnectivityList(i,:)];
                Iloc1=find(pairs(:,1)==i);
                pairs1=[pairs1; repmat(i1,length(Iloc1),1) pairs(Iloc1,2)];
                i1=i1+1;
                faces_1_part1_in_global=[faces_1_part1_in_global; i];
            elseif (tib.Points(facesTib_1(i,1),3)<0)&(tib.Points(facesTib_1(i,2),3)<0)&(tib.Points(facesTib_1(i,3),3)<0)
                facesTib_1_part2=[facesTib_1_part2; facesTib_1(i,:)];
                tib_part2.ConnectivityList=[tib_part2.ConnectivityList; tib.ConnectivityList(i,:)];
                Iloc2=find(pairs(:,1)==i);
                pairs2=[pairs2; repmat(i2,length(Iloc2),1) pairs(Iloc2,2)];
                i2=i2+1;
                faces_1_part2_in_global=[faces_1_part2_in_global; i];
            else
                keyboard;
            end
    end
    [Sum_Force_part1, Sum_Moments_part1, p_all1]=CalculateForces(x,y,z,psi,theta,phi, fem, tib_part1, pairs1, h, poisson, E, Rtib, R_aux, facesFem_1, facesTib_1_part1, centerfem);
    [Sum_Force_part2, Sum_Moments_part2, p_all2]=CalculateForces(x,y,z,psi,theta,phi, fem, tib_part2, pairs2, h, poisson, E, Rtib, R_aux, facesFem_1, facesTib_1_part2, centerfem);
    Sum_Force=Sum_Force_part1+Sum_Force_part2;
    Sum_Moments=Sum_Moments_part1+Sum_Moments_part2;
    
else
    [Sum_Force, Sum_Moments, p_all]=CalculateForces(x,y,z,psi,theta,phi, fem, tib, pairs, h, poisson, E, Rtib, R_aux, facesFem_1, facesTib_1, centerfem);
end

% if Options.optE
%     f=Function('f',{x,y,z,psi,theta,phi,E},{Sum_Force,Sum_Moments,p_all});
% else
%     if Options.sepML
%         f=Function('f',{x,y,z,psi,theta,phi},{Sum_Force,Sum_Moments,Sum_Force_part1,Sum_Force_part2,Sum_Moments_part1,Sum_Moments_part2,p_all1,p_all2});
%     else
%         f=Function('f',{x,y,z,psi,theta,phi},{Sum_Force,Sum_Moments,p_all});
%     end
% end
% if Options.optE
%     if Options.sepML
%         cg=CodeGenerator('CalculateKneeContact_varE_sepML_forTesting.c');
%     else
%         cg=CodeGenerator('CalculateKneeContact_varE_forTesting.c');
%     end
% else
%     if Options.sepML
%         cg=CodeGenerator('CalculateKneeContact_sepML_forTesting.c');
%     else
%         cg=CodeGenerator('CalculateKneeContact_forTesting.c');
%     end
% end
% cg.add(f);
% cg.add(f.jacobian);
% cg.generate;
out.Sum_Force_part1=Sum_Force_part1;
out.Sum_Force_part2=Sum_Force_part2;
out.Sum_Moments_part1=Sum_Moments_part1;
out.Sum_Moments_part2=Sum_Moments_part2;
out.p_all1=p_all1;
out.p_all2=p_all2;
keyboard;

end

function [Sum_Force, Sum_Moments, p_all]=CalculateForces(x,y,z,psi,theta,phi, fem, tib, pairs, h, poisson, E, Rtib, R_aux, facesFem, facesTib, centerfem)
import casadi.*

fem_Points=fem.Points+centerfem;
% Mtransformtib=fTransf(x,y,z,psi,theta,phi);
Mtransformtib=Rtib;

% pointsTib=MX(size(tib.Points,1),4);
for j=1:size(tib.Points,1)
    pointsTib(j,:)=Mtransformtib*[tib.Points(j,:) 1]';
end

originTib_G=Mtransformtib*[0 0 0 1]'; %origin of the tibia in ground coordinates
originTib_G=originTib_G(1:3)';

k=1; %number of contacting elements in femur
l=1; %count number of elements contacting element k of femur
d=zeros(size(pairs,1),3);
nt=zeros(size(pairs,1),3);
force_s_l=zeros(48,3); %for 189x50 -->97
Ctib=zeros(48,3);
mom_O=zeros(48,3);
p_all=zeros(48,1);
for i=1:size(pairs,1)

    try
        [d_aux,As,At,ns_aux,Cfem,Ctib_aux,nt_aux]=CalculateIntersection(fem_Points(facesFem(pairs(i,2),:),:),pointsTib(facesTib(pairs(i,1),:),1:3)); %first element in pairs is tibia, second femur
    catch
        keyboard;
    end
    d(i,:)=d_aux; %d is the distance between face tibia and face femur d=Cfem-Ctib
    nt(i,:)=nt_aux; %nt is the normal to the tibia face
    
    if i>1
        if ((pairs(i,1)==(pairs(i-1,1))))&&(i~=size(pairs,1))
            l=l+1;
        elseif ((pairs(i,1)~=(pairs(i-1,1))))||(i==size(pairs,1))
            try
            [mindist,nt_l]=CalculateMinimumDistance(d(i-l+1:i-1,:),nt(i-l+1:i-1,:));
            catch
                keyboard;
            end
            p=CalculatePressure(poisson,E,mindist,h);
            try
            force_s_l(k,:)=p*At*nt_l;
            p_all(k)=p;
            Ctib(k,:)=Ctib_aux;
            catch
                keyboard;
            end
            face_s(k)=pairs(i-1,1);
            mom_O(k,:)=cross(Ctib_aux-originTib_G,force_s_l(k,:));
            k=k+1;
            l=2;
        end
    else
        l=l+1;
    end
end
for k=k:48
    p_all(k)=0;
end

Sum_Force_G=sum(force_s_l);
Sum_Moments_G=sum(mom_O,1);

% Mrottib=full(f_Raux(x,y,z,psi,theta,phi));
Mrottib=R_aux;
Sum_Force=[(Mrottib')*Sum_Force_G']'; %in tibia frame
Sum_Moments=[(Mrottib')*Sum_Moments_G']'; %in tibia frame

end
function [d,As,Atib,ns,Cfem,Ctib,nt]=CalculateIntersection(fem_points,tib_points)

Ctib=mean(tib_points); %center triangle tibia
%edges floor/prism
edge1_p=tib_points(2,:)-tib_points(1,:);
edge2_p=tib_points(3,:)-tib_points(1,:);
%unit normal vector to triangle of tibia
nt=cross(edge1_p,edge2_p);
nt=nt/norm(nt); 
anglec=acos(sum(edge1_p.*edge2_p)/(norm(edge1_p)*norm(edge2_p)));
Atib=(1/2)*norm(edge1_p)*norm(edge2_p)*sin(anglec);

Cfem=mean(fem_points);
%edges femur
edge1_s=fem_points(2,:)-fem_points(1,:);
edge2_s=fem_points(3,:)-fem_points(1,:);
%unit normal vector to triangle of femur
ns=cross(edge1_s,edge2_s);
ns=ns/norm(ns);
% if ns(2)>0
%     keyboard;
% end
angles=acos(sum(edge1_s.*edge2_s)/(norm(edge1_s)*norm(edge2_s)));
As=(1/2)*norm(edge1_s)*norm(edge2_s)*sin(angles);

% hold all;
% plot3([Cfem(1) Cfem(1)+ns(1)*0.1],[Cfem(2) Cfem(2)+ns(2)*0.1],[Cfem(3) Cfem(3)+ns(3)*0.1]); %plot ray

dist_betweenCenters=Cfem-Ctib;

% d=sum(dist_betweenCenters.*nc)/sum(ns.*nc); % similar to Smith et al 2018
d=dist_betweenCenters;

end
function p=CalculatePressure(poisson,E,d,h)
    k=10000;
    pen=-d;
    %% Used for ISB wrong (natural knees)
%     dh=((1+tanh(k*(pen-h)))/2)*(pen/h); 
%     p=((1-poisson)*E/((1+poisson)*(1-2*poisson)))*log(1+dh); %from Bei and Fregly 2004 smoothed
    %% Used for TGCS? (artificial knees)
    p_init=((1-poisson)*E/((1+poisson)*(1-2*poisson)))*pen/h; %from Bei and Fregly 2004 smoothed
    p=p_init*(1+tanh(k*pen))/2;
end
function [mindist,nt_l]=CalculateMinimumDistance(d_v,nt_v)
proj=sum(d_v.*nt_v,2);

% dist=-sqrt(sum(d_v.^2,2)).*(dp./abs(dp)); %min distance, taking into account the sign, negative dp --> in contact
dist_v=-proj.*nt_v;
dist=-proj;
k=10000;
mindist=-[log(sum(exp(k*(dist))))/k]; %calculate minimum distance
nt_l=nt_v(1,:);

end
