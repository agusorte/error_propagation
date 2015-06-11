%Propagating error en projection matrix 
%considering radial distortion.
%the solution is Minimum Least Squares
%Agustin Ortega
%Oct 2012


function CovP=Compute_Cov_P_lambda(P,m,M,lamda,Cov_m,Cov_M,x0,y0)
%P is the projection Matrix
%m are the images points or lines
%M 3d Points or lines
%lamda is the radial distortion componet 3d Points or lines
%cov_m covariance of points
%cov_M Cov of 3d points
% the number of points must be mayor of 6
if nargin<=0,
    %for testing
    N=100;
    P=rand(3,4);%Projection Matrix
    m=ones(3,N);
    M=ones(4,N);
    
    m(1:2,:)=rand(2,N);%column vector [ui vi 1]
    M(1:3,:)=rand(3,N); %column vector [Xi Yi Zi 1]
    
    %covariances
    Cov_m=eye(2*N)*0.5;
    Cov_M=eye(3*N)*0.5;
    
    %luego probamos esto
    
    Cov_m_M=ones(3*N);%luego lo cambiamos
    Cov_m=ones(2*N);%luego lo cambiamos

    
    
    lamda=rand;
end;

N=size(M,2);

N2=size(m,2);
%validate number of points thes must be equal


if(N~=N2 || N<6)
    disp('The number points has to be more than 6')
    return;
end;

%compute Hessians

%matrix values
%matrix values
Sum_X=sum(M(1,:));
Sum_X_2=sum(M(1,:).^2);

Sum_Y=sum(M(2,:));
Sum_Y_2=sum(M(2,:).^2);

Sum_Z=sum(M(3,:));
Sum_Z_2=sum(M(3,:).^2);

Sum_X_Y=sum(M(1,:).*M(2,:));
Sum_X_Z=sum(M(1,:).*M(3,:));
Sum_Y_Z=sum(M(2,:).*M(3,:));

Sum_U_V_X=sum(-sqrt((m(1,:)-x0).^2+(m(2,:)-y0).^2).*M(1,:));

Sum_U_V_Y=sum(-sqrt((m(1,:)-x0).^2+(m(2,:)-y0).^2).*M(2,:));

Sum_U_V_Z=sum(-sqrt((m(1,:)-x0).^2+(m(2,:)-y0).^2).*M(3,:));

Sum_U_V=sum(-sqrt((m(1,:)-x0).^2+(m(2,:)-y0).^2));
Sum_U_V_2=sum(((m(1,:)-x0).^2+(m(2,:)-y0).^2));



%Hessian matrix dC^2/dP^2
HC_P=2*[Sum_X_2 Sum_X_Y Sum_X_Z Sum_X 0 0 0 0 0 0 0 0 0;
       Sum_X_Y Sum_Y_2 Sum_Y_Z Sum_Y 0 0 0 0 0 0 0 0 0;
       Sum_X_Z Sum_Y_Z Sum_Z_2 Sum_Z 0 0 0 0 0 0 0 0 0;
       Sum_X Sum_Y Sum_Z N 0 0 0 0 0 0 0 0 0;
       0 0 0 0 Sum_X_2 Sum_X_Y Sum_X_Z Sum_X 0 0 0 0 0;
       0 0 0 0 Sum_X_Y Sum_Y_2 Sum_Y_Z Sum_Y  0 0 0 0 0;
       0 0 0 0 Sum_X_Z Sum_Y_Z Sum_Z_2 Sum_Z 0 0 0 0 0;
       0 0 0 0 Sum_X Sum_Y Sum_Z N 0 0 0 0 0;
       0 0 0 0 0 0 0 0 Sum_X_2 Sum_X_Y Sum_X_Z Sum_X Sum_U_V_X;
       0 0 0 0 0 0 0 0 Sum_X_Y Sum_Y_2 Sum_Y_Z Sum_Y Sum_U_V_Y;
       0 0 0 0 0 0 0 0 Sum_X_Z Sum_Y_Z Sum_Z_2 Sum_Z Sum_U_V_Z;
       0 0 0 0 0 0 0 0 Sum_X Sum_Y Sum_Z N Sum_U_V;
       0 0 0 0 0 0 0 0 Sum_U_V_X Sum_U_V_Y Sum_U_V_Z Sum_U_V 2*Sum_U_V_2];
 
   
HC_P_Y=zeros(13,5*N);   
M_aux_i=[];
for i=1:N,
    X=M(1,i);
    Y=M(2,i);
    Z=M(3,i);
    u=m(1,i);
    v=m(2,i);
    
    
    %p11
    HC_P_Y(1,2*i-1)=-2*X;
    HC_P_Y(1,2*i)=0;
    
    HC_P_Y(1,2*N+3*i-2)=4*P(1,1)*X-2*u+2*P(1,2)*Y+2*P(1,3)*Z+2*P(1,4);
    HC_P_Y(1,2*N+3*i-1)=2*P(1,2)*X;
    HC_P_Y(1,2*N+3*i)=2*P(1,3)*X;
    
    %p12
    HC_P_Y(2,2*i-1)=-2*Y;
    HC_P_Y(2,2*i)=0;
    
    HC_P_Y(2,2*N+3*i-2)=2*P(1,1)*Y;
    HC_P_Y(2,2*N+3*i-1)=4*P(1,2)*Y-2*u+2*P(1,1)*X+2*P(1,2)*Z+2*P(1,4);
    HC_P_Y(2,2*N+3*i)=2*P(1,3)*Y;
    
    %p13
    
    HC_P_Y(3,2*i-1)=-2*Z;
    HC_P_Y(3,2*i)=0;
    
    HC_P_Y(3,2*N+3*i-2)=2*P(1,1)*Z;
    HC_P_Y(3,2*N+3*i-1)=2*P(1,2)*Z;
    HC_P_Y(3,2*N+3*i)=4*P(1,3)*Z-2*u+2*P(1,1)*X+2*P(1,2)*Y+2*P(1,4);
    
    %p14
    HC_P_Y(4,2*i-1)=-2;
    HC_P_Y(4,2*i)=0;
    
    HC_P_Y(4,2*N+3*i-2)=2*P(1,1);
    HC_P_Y(4,2*N+3*i-1)=2*P(1,2);
    HC_P_Y(4,2*N+3*i)=2*P(1,3);
    
    %p21
    HC_P_Y(5,2*i-1)=0;
    HC_P_Y(5,2*i)=-2*X;
    
    HC_P_Y(5,2*N+3*i-2)=4*P(2,1)*X-2*v+2*P(2,2)*Y+2*P(2,3)*Z+2*P(2,4);
    HC_P_Y(5,2*N+3*i-1)=2*P(2,2)*X;
    HC_P_Y(5,2*N+3*i)=2*P(2,3)*X;
    
   %p22
   HC_P_Y(6,2*i-1)=-2*X;
    HC_P_Y(6,2*i)=-2*Y;
    
    HC_P_Y(6,2*N+3*i-2)=2*P(2,1)*Y;
    HC_P_Y(6,2*N+3*i-1)=4*P(2,2)*Y-2*v+2*P(2,1)*X+2*P(2,3)*Z+2*P(2,4);
    HC_P_Y(6,2*N+3*i)=2*P(2,3)*Y;
    
   %p23
    HC_P_Y(7,2*i-1)=0;
    HC_P_Y(7,2*i)=-2*Z;
    
    HC_P_Y(7,2*N+3*i-2)=2*P(2,1)*Z;
    HC_P_Y(7,2*N+3*i-1)=2*P(2,2)*Z;
    HC_P_Y(7,2*N+3*i)=4*P(2,3)*Z-2*v+2*P(2,1)*X+2*P(2,2)*Y+2*P(2,4);

    %P24
    HC_P_Y(8,2*i-1)=0;
    HC_P_Y(8,2*i)=-2;
    
    HC_P_Y(8,2*N+3*i-2)=2*P(2,1);
    HC_P_Y(8,2*N+3*i-1)=2*P(2,2);
    HC_P_Y(8,2*N+3*i)=2*P(2,3);
    
    %P31
    HC_P_Y(9,2*i-1)=-2*lamda*u*X/sqrt((u-x0)^2+(v-y0)^2);
    HC_P_Y(9,2*i)=-2*lamda*v*X/sqrt((u-x0)^2+(v-y0)^2);
    
    HC_P_Y(9,2*N+3*i-2)=4*P(3,1)*X-2-2*lamda*sqrt((u-x0)^2+(v-y0)^2)+2*P(3,2)*Y+2*P(3,3)*Z+2*P(3,4);
    HC_P_Y(9,2*N+3*i-1)=2*P(3,2)*X;
    HC_P_Y(9,2*N+3*i)=2*P(3,3)*X;
    
    %P32
    HC_P_Y(10,2*i-1)=-2*lamda*u*Y/sqrt((u-x0)^2+(v-y0)^2);
    HC_P_Y(10,2*i)=-2*lamda*v*Y/sqrt((u-x0)^2+(v-y0)^2);
    
    HC_P_Y(10,2*N+3*i-2)=-2*P(3,1)*Y;
    HC_P_Y(10,2*N+3*i-1)=4*P(3,2)*Y-2-2*lamda*sqrt((u-x0)^2+(v-y0)^2)+2*P(3,1)*X+2*P(3,3)*Z+2*P(3,4);
    HC_P_Y(10,2*N+3*i)=2*P(3,3)*Y;
    
   %p33
    HC_P_Y(11,2*i-1)=-2*lamda*u*Z/sqrt((u-x0)^2+(v-y0)^2);
    HC_P_Y(11,2*i)=-2*lamda*v*Z/sqrt((u-x0)^2+(v-y0)^2);
    
    HC_P_Y(11,2*N+3*i-2)=2*P(3,1)*Z;
    HC_P_Y(11,2*N+3*i-1)=2*P(3,2)*Z;
    HC_P_Y(11,2*N+3*i)=4*P(3,3)*Z-2-2*lamda*sqrt((u-x0)^2+(v-y0)^2)+2*P(3,1)*X+...
        2*P(3,2)*Y+2*P(3,4);
   
    %p34
    HC_P_Y(12,2*i-1)=-2*lamda*u/sqrt((u-x0)^2+(v-y0)^2);
    HC_P_Y(12,2*i)=-2*lamda*v/sqrt((u-x0)^2+(v-y0)^2);
    
    HC_P_Y(12,2*N+3*i-2)=2*P(3,1);
    HC_P_Y(12,2*N+3*i-1)=2*P(3,2);
    HC_P_Y(12,2*N+3*i)=2*P(3,3);
    
    %lambdad
    HC_P_Y(13,2*i-1)=2*lamda*u+2*(1+lamda*sqrt((u-x0)^2+(v-y0)^2)-P(3,1)*X-...
        P(3,2)*Y-P(3,3)*Z-P(3,4))*u/sqrt((u-x0)^2+(v-y0)^2);
    HC_P_Y(13,2*i)=2*lamda*v+2*(1+lamda*sqrt((u-x0)^2+(v-y0)^2)-P(3,1)*X-...
        P(3,2)*Y-P(3,3)*Z-P(3,4))*v/sqrt((u-x0)^2+(v-y0)^2);
    
    HC_P_Y(13,2*N+3*i-2)=-2*P(3,1)*sqrt((u-x0)^2+(v-y0)^2);
    HC_P_Y(13,2*N+3*i-1)=-2*P(3,2)*sqrt((u-x0)^2+(v-y0)^2);
    HC_P_Y(13,2*N+3*i)=-2*P(3,3)*sqrt((u-x0)^2+(v-y0)^2);
    
    
end;

% fill covariance matrix is 5Nx5N

 Cov_m_M(1:2*N,1:2*N)= Cov_m;
 Cov_m_M(2*N+1:5*N,2*N+1:5*N)= Cov_M;
 

alphaF=-inv(HC_P)*HC_P_Y;


CovP=alphaF*Cov_m_M*alphaF';
%AlphaF=
disp('');

   

   