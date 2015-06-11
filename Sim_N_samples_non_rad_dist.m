%Simulation of samples N experiments
%we dont consider radial distortion
% for computing the projection matrix
% Agustin Ortega
% Oct 2012

clc,clear all, close all;

%% dir 

addpath('/home/aortega/proyectos/calibracion/IST_Actual/matlab.my/util');


%% data of calibration parameters
R = [ -0.0567   -0.9692    0.2397;
        -0.7991    0.1879    0.5710;
        -0.5985   -0.1592   -0.7852];

T=[0.5; -0.1; 8];


alpha = 760;         % horizontal focal length
beta = 760;         % vertical focal length

a = 320;                  % principal point u coord
b = 240;                  % principal point v coord

f=[alpha beta]'; %focal length

K= [ alpha 0 a ;...
    0 beta b ; ...
    0 0 1 ];

kd=zeros(1,5);%distortion radial

P_real=K*[R T];%projection matrix


% number of samples
N_sim=200;


% 3d lines that belong to a cube
cube_points= [0,0,0;0,0,1;1,0,0;1,0,1;1,1,0;1,1,1;0,1,0;0,1,1;0,0,0;...
              0,1,0;1,0,0;1,1,0;1,0,1;1,1,1;0,0,1;0,1,1;0,0,0;1,0,0;%
              0,0,1;1,0,1;0,1,1;1,1,1;0,1,0;1,1,0;]';       
          
M_real=[cube_points;ones(1,size(cube_points,2))];%3d lines

%3d points we think to generate to cube structure

N=size(M_real,2);

m_real=P_real*M_real;

m_real=deshomogenizar(m_real);




%%
% noise adding to lines 3d, image points and projection matrix
% the noise is distibuted as Gaussian
% the idea is generate n samples

sigma_m=2.0;
sigma_M=0.001;
sigma_T=0.01;

Cov_m=eye(N*2)*sigma_m;
Cov_M=eye(N*3)*sigma_M;

%%
% error propgation
CovP=Compute_Cov_P(P_real,m_real,M_real,Cov_m,Cov_M);

[K_, R_, C_method, pp_, pv_] = decomposecamera(P_real);
 

[Cov_t]=Cov_position(P_real,CovP);

%%
% principal bucle
    
m_noisy=[];
M_noisy=[];
cov_P_sample=[];
C_sal=[];
NEES=[];
for i=1:N_sim,
    
    cube_noise=cube_points+randn(3,N)*sigma_M;
    M=[cube_noise;ones(1,size(cube_points,2))];%3d lines
    
    m=P_real*M_real;%real o noisy ?
    
    m=deshomogenizar(m);
    
    m(1:2,:)=m(1:2,:)+randn(2,N)*sigma_m;
    
    
    
    m_noisy=[m_noisy  m];
    M_noisy=[M_noisy  M];
    
    P_est=CalibNormDLT(m, M);
    

    
    [K_est, R_est, C_est, pp, pv] = decomposecamera(P_est);
    
    C_sal=[C_sal  C_est];

    %here we compute the cov per each point
    %[Cov_t]=Cov_position(P_est,CovP);
    
    %NEES error
    e=C_est- C_method;
    
    NESS_sample=e'*inv(Cov_t(1:3,1:3))*e;
    NEES=[NEES NESS_sample];
    
end;


%visualize data real and samples
figure(333);
title('3D camera pose covariance');
hold on;

for i=1:2:N,
    
plot3([M_real(1,i) M_real(1,i+1)],...
    [M_real(2,i) M_real(2,i+1)],...
    [M_real(3,i) M_real(3,i+1)],'b-.');

end;

plot3(M_noisy(1,:), M_noisy(2,:),M_noisy(3,:),'r.');
plot3(C_sal(1,:), C_sal(2,:),C_sal(3,:),'g.');

hFOV = 40;%2*atan(cam.cc(1) / cam.fc(1)) * 180/pi;  % not correct...

options.hFOV = hFOV;
options.scale = 5;%camera scale

%compute aproximate covariance montecalor 
cov_montecarlo=cov(C_sal');
C_montecarlo=mean(C_sal')';

NEES2=[];
for i=1:size(C_sal,2)
    e=C_sal(1:3,i) - C_montecarlo;
    
    NESS_sample=e'*inv(cov_montecarlo(1:3,1:3))*e;
    NEES2=[NEES2 NESS_sample];
    
end



[K, r0, t0] = draw_camera(P_real,options);



plotGMM3D([ C_method(1),  C_method(2),  C_method(3)]',...
    Cov_t(1:3,1:3), [1 0.01 0], 0.3,1)

plotGMM3D([ C_montecarlo(1),  C_montecarlo(2),  C_montecarlo(3)]',...
    cov_montecarlo, [0 0.01 1], 0.3,1)

axis equal;

figure(444);

hold on;
title('Image points');
for i=1:2:N,
plot([m_real(1,i), m_real(1,i+1)],...
    [m_real(2,i), m_real(2,i+1)],'b-.');

end;

plot(m_noisy(1,:), m_noisy(2,:),'r.');


axis equal;

%% NEES results

figure(332);
hold on;
title('NEES Linear error propagation');
plot(1:length(NEES),NEES,'b', ...
    1:length(NEES),ones(1,length(NEES))*7.82,'r--' );
xlabel('t');
ylabel('NEES');
nees = sum(NEES<7.82)/length(NEES)

figure(334);
hold on;
title('NEES Montecarlo');
plot(1:length(NEES2),NEES2,'b', ...
    1:length(NEES2),ones(1,length(NEES2))*7.82,'r--' );
xlabel('t');
ylabel('NEES');
nees2 = sum(NEES2<7.82)/length(NEES2)


%
%measure book first order prediction

measure=norm(Cov_t(1:3,1:3)-cov_montecarlo)/norm(cov_montecarlo)


