%simulation of error propagation of projection matrix
%using points without radial distortion.
%By Agustin Ortega
%oct 2012

clc,clear all,close all;

%% add libraries

addpath('/home/aortega/proyectos/calibracion/IST_Actual/matlab.my/util');
%%
 R = [ -0.0567   -0.9692    0.2397;
        -0.7991    0.1879    0.5710;
        -0.5985   -0.1592   -0.7852];


T=[0.5; -0.1; 8];


alpha = 760;         % horizontal focal length
beta = 760;         % vertical focal length

a = 320;                  % principal point u coord
b = 240;                  % principal point v coord

f=[alpha beta]';

K= [ alpha 0 a ;...
    0 beta b ; ...
    0 0 1 ];

kd=zeros(1,5);%distortion radial

P=K*[R T];


%% First generate data
cube_points= [0,0,0;0,0,1;1,0,0;1,0,1;1,1,0;1,1,1;0,1,0;0,1,1;0,0,0;...
              0,1,0;1,0,0;1,1,0;1,0,1;1,1,1;0,0,1;0,1,1;0,0,0;1,0,0;%
              0,0,1;1,0,1;0,1,1;1,1,1;0,1,0;1,1,0;]';       
          
M=[cube_points;ones(1,size(cube_points,2))];
%3d points we think to generate to cube structure

N=size(M,2);

m=P*M;

%% deshomogenizar m

for i=1:size(m,2),
    m(:,i)=m(:,i)/m(3,i);
end;




%% find the solution of projection matrix using lines

P=CalibNormDLT(m, M);

%% compute covariance error of proyection matrix
Cov_m=eye(N*2)*0.01;
Cov_M=eye(N*3)*0.0001;

CovP=Compute_Cov_P(P,m,M,Cov_m,Cov_M);

%% descompose P
[K_est, R_est, C_est, pp, pv] = decomposecamera(P);

K_est=K_est/K_est(3,3);
%%covariance of position

[Cov_t]=Cov_position(P,CovP)


%visualize :)
hFOV = 40;%2*atan(cam.cc(1) / cam.fc(1)) * 180/pi;  % not correct...
vFOV = 40;%2*atan(

%3d points
figure(33);

hold on;
for i=1:2:N,
    
plot3([M(1,i) M(1,i+1)],[M(2,i) M(2,i+1)],[M(3,i) M(3,i+1)],'b-.');

end;
options.hFOV = hFOV;
options.scale = 5;%camera scale

[K, r0, t0] = draw_camera(P,options);
plotGMM3D([ C_est(1),  C_est(2),  C_est(3)]',...
    Cov_t(1:3,1:3), [0 0.01 1], 0.1,3)

axis equal;

figure(222)
hold on;
for i=1:2:N,
plot([m(1,i), m(1,i+1)],[m(2,i), m(2,i+1)],'r-.');

end;


axis equal;
