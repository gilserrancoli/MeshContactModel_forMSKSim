function fTransf = create_fTransf()

import casadi.*

%Function for rotation 3x3
psi=MX.sym('psi',1);
theta=MX.sym('theta',1);
phi=MX.sym('phi',1);

%function for translation 4x4
x=MX.sym('x',1);
y=MX.sym('y',1);
z=MX.sym('z',1);

%function for translation 0.042 4x4
Rtrans0042=[1 0 0 0; 0 1 0 0.042; 0 0 1 0; 0 0 0 1];

%function for rotation 4x4
R1=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R2=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
R3=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
R=R1*R2*R3;

R_aux=R1*R2*R3;
Rrot=MX(4,4);
Rrot(1:3,1:3)=R_aux;
Rrot(4,4)=1;

%translation in tibia coord system
aux=R*[-x;-y;-z];
Rtranstib_4x4=MX(4,4);
Rtranstib_4x4(1,1)=1;
Rtranstib_4x4(2,2)=1;
Rtranstib_4x4(3,3)=1;
Rtranstib_4x4(4,4)=1;
Rtranstib_4x4(1:3,4)=aux;

Rtib=Rtranstib_4x4*Rtrans0042*Rrot;
fTransf=Function('fTransf',{x,y,z,psi,theta,phi},{Rtib});

end