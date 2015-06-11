%Simulation of samples N experiments
%Here we consider raial distortion
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


alpha = 800;         % horizontal focal length
beta = 800;         % vertical focal length

a = 320;                  % principal point u coord
b = 240;                  % principal point v coord

f=[alpha beta]'; %focal length

K= [ alpha 0 a ;...
    0 beta b ; ...
    0 0 1 ];

kd=zeros(1,5);%distortion radial

P_real=K*[R T];%projection matrix

r_lambda=0;%this is the radial distortion

% number of samples
N_sim=200;


% 3d lines that belong to a cube
cube_points= [0,0,0;0,0,1;1,0,0;1,0,1;1,1,0;1,1,1;0,1,0;0,1,1;0,0,0;...
              0,1,0;1,0,0;1,1,0;1,0,1;1,1,1;0,0,1;0,1,1;0,0,0;1,0,0;%
              0,0,1;1,0,1;0,1,1;1,1,1;0,1,0;1,1,0;]';       
          
M_real=[cube_points;ones(1,size(cube_points,2))];%3d lines

%3d points we think to generate to cube structure

N=size(M_real,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%here change the projection now we introduce radial distortion 
%in our model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_real=P_real*M_real; %points without radial distortion


m_real=deshomogenizar(m_real);


%adding radial distortion

md=m_real+r_lambda*[zeros(1,N);...
                  zeros(1,N);...
                  ((m_real(1,:)-320).^2+(m_real(2,:)-240).^2)]; %center of projection
% 

md=deshomogenizar(md);
md(1:2,:)=md(1:2,:)+m_real(1:2,:);

% [P_est r_lamda_est]=Compute_P_lambda(M_real,mu)

%%
% noise adding to lines 3d, image points and projection matrix
% the noise is distibuted as Gaussian
% the idea is generate n samples

sigma_m=1.5;
sigma_M=0.001;


Cov_m=eye(N*2)*sigma_m;
Cov_M=eye(N*3)*sigma_M;

%%
% error propgation

%this covariance includes the radial distortion
CovP=Compute_Cov_P_lambda(P_real,md,M_real,r_lambda,Cov_m,Cov_M,320,240);
[K_, R_, C_method, pp_, pv_] = decomposecamera(P_real);
[Cov_t]=Cov_position(P_real,CovP(1:12,1:12));

%%
% principal bucle
    K
m_noisy=[];
M_noisy=[];
cov_P_sample=[];
C_sal=[];
C_sal2=[];
C_sal3=[];
C_sal4=[];
lamda_sal=[];
NEES=[];
for i=1:N_sim,
    
    
    %%ading noise to 3d lines
    cube_noise=cube_points+randn(3,N)*sigma_M;
    M=[cube_noise;ones(1,size(cube_points,2))];%3d lines
    
    m=P_real*M_real;%real o noisy ?
    
    m=deshomogenizar(m);
    
    
%     m(1:2,:)=m(1:2,:)+randn(2,N)*sigma_m;

    
    md_=m+...
        r_lambda*[zeros(1,N); ...
                  zeros(1,N);...
                  (((m(1,:)-320).^2+(m(2,:)-240).^2))];
              
    md_=m   ;       
              
    md_=deshomogenizar(md_);          
    
    md_(1:2,:)=md_(1:2,:)+m(1:2,:);
    
    md_(1:2,:)=md(1:2,:)+randn(2,N)*sigma_m;
    
    m_noisy=[m_noisy  md_];
    M_noisy=[M_noisy  M];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %here compute the estimation method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %P_est=CalibNormDLT(m, M);
    [P_est1,P_est2,P_est3,P_est4, r_lambda_]=Compute_P_lambda(M,md_);
    
    lamda_sal=[lamda_sal 1/r_lambda_];%here we'll see the behavior of lambda
    [K_est, R_est, C_est, pp, pv] = decomposecamera(P_est1);
    [K_est, R_est, C_est2, pp, pv] = decomposecamera(P_est2);
    [K_est, R_est, C_est3, pp, pv] = decomposecamera(P_est3);
    [K_est, R_est, C_est4, pp, pv] = decomposecamera(P_est4);
    
    C_sal=[C_sal  C_est];
    C_sal2=[C_sal2  C_est2];
    C_sal3=[C_sal3  C_est3];
    C_sal4=[C_sal4  C_est4];

    
    
    %NEES error
    e=C_est- C_method;
    
    NESS_sample=e'*inv(Cov_t(1:3,1:3))*e;
    NEES=[NEES NESS_sample];
    
end;


%visualize data real and samples
figure(333);
hold on;

for i=1:2:N,
    
plot3([M_real(1,i) M_real(1,i+1)],...
    [M_real(2,i) M_real(2,i+1)],...
    [M_real(3,i) M_real(3,i+1)],'b-.');

end;

plot3(M_noisy(1,:), M_noisy(2,:),M_noisy(3,:),'r.');
plot3(C_sal(1,:), C_sal(2,:),C_sal(3,:),'g.');
plot3(C_sal2(1,:), C_sal2(2,:),C_sal2(3,:),'m.');
plot3(C_sal3(1,:), C_sal3(2,:),C_sal3(3,:),'y.');
plot3(C_sal4(1,:), C_sal4(2,:),C_sal4(3,:),'k.');

hFOV = 40;%2*atan(cam.cc(1) / cam.fc(1)) * 180/pi;  % not correct...

options.hFOV = hFOV;
options.scale = 5;%camera scale

%compute aproximate covariance montecalor 
cov_montecarlo=cov(C_sal');
C_montecarlo=mean(C_sal');

[K, r0, t0] = draw_camera(P_real,options);



plotGMM3D([ C_method(1),  C_method(2),  C_method(3)]',...
    Cov_t(1:3,1:3), [0 0.01 1], 0.3,1)

plotGMM3D([ C_montecarlo(1),  C_montecarlo(2),  C_montecarlo(3)]',...
    cov_montecarlo, [1 0.01 0], 0.3,1)

axis equal;

figure(444);

hold on;
plotting_line_sampled(20,md,r_lambda)

for i=1:2:N,
% plot([md(1,i), md(1,i+1)],...
%     [md(2,i), md(2,i+1)],'b-.');

plot([m_real(1,i), m_real(1,i+1)],...
    [m_real(2,i), m_real(2,i+1)],'b-.');

end;


plot(m_noisy(1,:), m_noisy(2,:),'r.');

% 
% axis equal;
% 
% NEES2=[];
% for i=1:size(C_sal,2)
%     e=C_sal(1:3,i) - C_montecarlo;
%     
%     NESS_sample=e'*inv(cov_montecarlo(1:3,1:3))*e;
%     NEES2=[NEES2 NESS_sample];
%     
% end



%% NEES results

figure(332);
title('NEES graphic');
plot(1:length(NEES),NEES,'b', ...
    1:length(NEES),ones(1,length(NEES))*5.99,'r--' );
xlabel('t');
ylabel('NEES');


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


figure(444)
title('estimation radial distortion');

