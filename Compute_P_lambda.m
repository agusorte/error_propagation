%function that computes P and radial distortion
%Gaspar Paper
%by Agustin Ortega
%oct 2012

function [P_p1,P_p2,P_p3,P_p4,lambda_]=Compute_P_lambda(M,md)


if nargin<=0,
    close all;
    N=6;
    P = [-236.88,-826.304,-59.504,2960;
        -782.92,112.112,268.352,1840;
        -0.5985,-0.1592,-0.7852,8;]
    
%     M=rand(4,N);
%     m=P*M;
    
    cube_points= [0,0,0;0,0,1;1,0,0;1,0,1;1,1,0;1,1,1;0,1,0;0,1,1;0,0,0;...
        0,1,0;1,0,0;1,1,0;1,0,1;1,1,1;0,0,1;0,1,1;0,0,0;1,0,0;%
        0,0,1;1,0,1;0,1,1;1,1,1;0,1,0;1,1,0;]';
    
    M=[cube_points;ones(1,size(cube_points,2))];%3d lines
    
    %3d points we think to generate to cube structure
    
    N=size(M,2);
    
    r_lambda=0.05;
    
    mu=P*M;
    
    [K_est, R_est, C_est, pp, pv] =decomposecamera(P);
    
    mu=deshomogenizar(mu);%undistorted point
    
    mu_centered(1,:)=mu(1,:)-320;
    mu_centered(2,:)=mu(2,:)-240;

    md(1,:)=(mu(1,:))./(1+r_lambda*(mu_centered(1,:).^2+mu_centered(2,:).^2))+mu(1,:);
    md(2,:)=(mu(2,:))./(1+r_lambda*(mu_centered(1,:).^2+mu_centered(2,:).^2))+mu(2,:);
    md(3,:)=ones(1,N);
    %plot for debbuging 
    
%     figure(444);
%     
%     hold on;
%     for i=1:2:N,
%         plot([mu(1,i), mu(1,i+1)],...
%             [mu(2,i), mu(2,i+1)],'b-.');
%         
%         plot([md(1,i), md(1,i+1)],...
%             [md(2,i), md(2,i+1)],'g-.');
%         
%     end;
    
    %DEBUG ONLY TESTING THE RADIAL DISTORTION GENERATE 1 LINES 
    %WITH NPOINTS 
    figure(2343); hold on;
    
      for i=1:2:N,
        plot([mu(1,i), mu(1,i+1)],...
            [mu(2,i), mu(2,i+1)],'b-.');
        
        plot([md(1,i), md(1,i+1)],...
            [md(2,i), md(2,i+1)],'g-.');
        
    end;
    
    for i=1:2:N,
        N_lines=10;
        
        line=cross(mu(:,i),mu(:,i+1));
        x=mu(1,i):1/N_lines:mu(1,i+1);
        
        
        y=-(line(3)+line(1)*x)/line(2);
        
        
        xd=x./(1+r_lambda*((x-320).^2+(y-240).^2))+x;
        yd=y./(1+r_lambda*((x-320).^2+(y-240).^2))+y;
       
  
        
       % plot(x,y,'r-.');
      plot(xd,yd,'r.-');
        
    end;
   axis([260 420 140 320])
    
    
    
    
end;

% %% centroids of the points
centroid1 = mean(md(1:2,:)')';
centroid2 = mean(M(1:3,:)')';

%% Shift the origin of the points to the centroid
md(1,:) = md(1,:) - centroid1(1);
md(2,:) = md(2,:) - centroid1(2);

M(1,:) = M(1,:) - centroid2(1);
M(2,:) = M(2,:) - centroid2(2);
M(3,:) = M(3,:) - centroid2(3);

%% Normalize the points so that the average distance from the origin is equal to sqrt(2).
averagedist1 = mean(sqrt(md(1,:).^2 + md(2,:).^2));
averagedist2 = mean(sqrt(M(1,:).^2 + M(2,:).^2 + M(3,:).^2));

scale1 = sqrt(2)/averagedist1;
scale2 = sqrt(3)/averagedist2;

md(1:2,:) = scale1*md(1:2,:);
M(1:3,:) = scale2*M(1:3,:);

%% similarity transform 1
T1 = [scale1         0  -scale1*centroid1(1)
           0    scale1  -scale1*centroid1(2)
           0         0                     1];

%% similarity transform 2
T2 = [scale2       0      0   -scale2*centroid2(1)
           0  scale2      0   -scale2*centroid2(2)
           0       0  scale2  -scale2*centroid2(3)
           0       0       0                     1];

% [m, T1] = normalise2dpts(m);
% [M, T2] = normalise2dpts(M);
    
       
N=size(M,2); %M is 4xN
%now compute matrices
% matrix size is 3Nx12
Ai1=[];
Ai2=[];

for i=1:N,

  A_1_aux=kron(M(:,i)',vec2skew(md(:,i)));
  
  Ai1=[Ai1;  A_1_aux];
  
  A_2_aux=kron(M(:,i)',vec2skew([0 0 sqrt(md(1,i)^2+md(2,i)^2)]));
  
  Ai2=[Ai2;  A_2_aux];
   
end;

AiT=Ai1'*Ai1;
Ai1_A2T=Ai1'*Ai2;


[P_est,e,s]=polyeig(AiT,Ai1_A2T);

e_non_inf=find(~isinf(e));


if length(e_non_inf)>4, % we have 3 solutions
    
    % [val i]=min(e);
    [val i]=min(abs(e)); %find the closets value to 0
    
    lambda_=val;
    P_mat=reshape(P_est(:,e_non_inf(1)),3,4);
    P_mat2=reshape(P_est(:,e_non_inf(2)),3,4);
    P_mat3=reshape(P_est(:,e_non_inf(3)),3,4);
    P_mat4=reshape(P_est(:,e_non_inf(4)),3,4);
    
    %P_p=T2\P_mat*T1;
    P_p1=inv(T1)*P_mat*T2;
    P_p2=inv(T1)*P_mat2*T2;
    P_p3=inv(T1)*P_mat3*T2;
    P_p4=inv(T1)*P_mat4*T2;
    
    
    %P_p=P_mat;
    
    %debugging
    %  [K_p, R_p, C_p, pp_p, pv_p] = decomposecamera(P_p);
    %  [K, R, C, pp, pv] = decomposecamera(P);
    %  disp(K_p/K_p(3,3)); disp(K);
    %  disp(R_p); disp(R);
    %  disp(C_p); disp(C);
end;
disp('');

