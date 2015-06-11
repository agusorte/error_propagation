%symbolic error propgation for the pose object.

syms p11 p12 p13 p14 p21 p22 p23 p24 p31 p32 p33 p34

P=[p11 p12 p13 p14;
    p21 p22 p23 p24;
    p31 p32 p33 p34];


v_x = [p12 p22 p32, p13, p23, p33,p14,p24,p34];
v_y = [p11 p21 p31, p13, p23, p33,p14,p24,p34];
v_z = [p11 p21 p31, p12, p22, p32,p14,p24,p34];
v_t = [p11 p21 p31, p12, p22, p32,p13,p23,p33];

P1=[p11;...
    p21; ...
    p31];

P2=[p12;...
    p22; ...
    p32 ];

P3=[p13 ;...
    p23; ...
    p33];

P4=[p14 ;...
    p24; ...
    p34 ];

X=det([P2 P3 P4]);
Y=-det([P1 P3 P4]);
Z=det([P1 P2 P4]);
T=-det([P1 P2 P3]);

C=[X/T; Y/T; Z/T; 1];
C2=[X; Y; Z; T];
    

J_x=jacobian(X,v_x);
J_y=jacobian(X,v_y);
J_z=jacobian(X,v_z);
J_t=jacobian(X,v_t);

J_C=jacobian(C,[p11, p12, p13, p14, p21,p22,p23,p24, p31,p32,p33,p34])
J_C2=jacobian(C2,[p11, p12, p13, p14, p21,p22,p23,p24, p31,p32,p33,p34])

%mirar esta parte del proceso y calcular bien el jacobiano
% primero R1=Qz*P(1:3,1:3)
%luego R2=Qy(R1)R1
%R3=Qx(R2)R2

Qz=[p11/sqrt(p11^2+p12^2) -p21/sqrt(p11^2+p12^2) 0;
    p21/sqrt(p11^2+p12^2) p11/sqrt(p11^2+p12^2)  0
     0 0 1];
Qy=[p11/sqrt(p11^2+p12^2)  0 -p21/sqrt(p11^2+p12^2);
    0 1 0;
    p11/sqrt(p11^2+p12^2)  0 -p21/sqrt(p11^2+p12^2);];

Qx=[1 0 0;
    0 p11/sqrt(p11^2+p12^2) -p21/sqrt(p11^2+p12^2);
    0 p21/sqrt(p11^2+p12^2) p11/sqrt(p11^2+p12^2)];

Q=Qx*Qy*Qz

Q_vec=reshape(Q,1,9);

J_Q=jacobian(Q,[p11 p12 p13 p14 p21 p22 p23 p24 p31 p32 p33 p34])

R=Q*P(1:3,1:3)

R_vec=reshape(R,1,9);

J_R=jacobian(R_vec,[p11 p12 p13 p14 p21 p22 p23 p24 p31 p32 p33 p34]);


%QR error propagation
