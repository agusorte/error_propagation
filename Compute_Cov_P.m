%this program propagate the error using the the matrix the proyection and
%points 2d and lines 3d
%the solution of Minumum Least Square is 
% in the program is not considered the radial distortion.
%By Agustin Ortega
%Oct 2012

function CovP=Compute_Cov_P(P,m,M,Cov_m,Cov_M)

%P is the projection Matrix
%m are the images points or lines
%M 3d Points or lines
%cov_m covariance of points
%cov_M Cov of 3d points
% the number of points must be mayor of 6
if nargin<=0,
    %for testing
    N=6;
    P=rand(3,4);%Projection Matrix
    m=ones(3,N);
    M=ones(4,N);
    
    m(1:2,:)=rand(2,N);%column vector [ui vi 1]
    M(1:3,:)=rand(3,N); %column vector [Xi Yi Zi 1]
    
    %covariances
    Cov_m=eye(2*N)*0.5;
    Cov_M=eye(3*N)*0.5;
    
    %luego probamos esto
    
    Cov_m_M=zeros(5*N);%luego lo cambiamos
end;

%validate number of points thes must be equal

%

N=size(M,2);

N2=size(m,2);

if(N~=N2 || N<6)
    disp('The number points has to be more than 6')
    return;
end;

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





%Hessian matrix dC^2/dP^2
HC_P=2*[Sum_X_2 Sum_X_Y Sum_X_Z Sum_X 0 0 0 0 0 0 0 0;
       Sum_X_Y Sum_Y_2 Sum_Y_Z Sum_Y 0 0 0 0 0 0 0 0;
       Sum_X_Z Sum_Y_Z Sum_Z_2 Sum_Z 0 0 0 0 0 0 0 0;
       Sum_X Sum_Y Sum_Z N 0 0 0 0 0 0 0 0;
       0 0 0 0 Sum_X_2 Sum_X_Y Sum_X_Z Sum_X 0 0 0 0;
       0 0 0 0 Sum_X_Y Sum_Y_2 Sum_Y_Z Sum_Y  0 0 0 0;
       0 0 0 0 Sum_X_Z Sum_Y_Z Sum_Z_2 Sum_Z 0 0 0 0;
       0 0 0 0 Sum_X Sum_Y Sum_Z N 0 0 0 0;
       0 0 0 0 0 0 0 0 Sum_X_2 Sum_X_Y Sum_X_Z Sum_X;
       0 0 0 0 0 0 0 0 Sum_X_Y Sum_Y_2 Sum_Y_Z Sum_Y;
       0 0 0 0 0 0 0 0 Sum_X_Z Sum_Y_Z Sum_Z_2 Sum_Z;
       0 0 0 0 0 0 0 0 Sum_X Sum_Y Sum_Z N];
   
   
   
HC_P_Y=zeros(12,5*N);   
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
    HC_P_Y(9,2*i-1)=0;
    HC_P_Y(9,2*i)=0;
    
    HC_P_Y(9,2*N+3*i-2)=4*P(3,1)*X-2+2*P(3,2)*Y+2*P(3,3)*Z+2*P(3,4);
    HC_P_Y(9,2*N+3*i-1)=2*P(3,2)*X;
    HC_P_Y(9,2*N+3*i)=2*P(3,3)*X;
    
    %P32
    HC_P_Y(10,2*i-1)=0;
    HC_P_Y(10,2*i)=0;
    
    HC_P_Y(10,2*N+3*i-2)=2*P(3,1)*Y;
    HC_P_Y(10,2*N+3*i-1)=4*P(3,2)*Y-2+2*P(3,1)*X+2*P(3,3)*Z+2*P(3,4);
    HC_P_Y(10,2*N+3*i)=2*P(3,3)*Y;
    
   %p33
    HC_P_Y(11,2*i-1)=0;
    HC_P_Y(11,2*i)=0;
    
    HC_P_Y(11,2*N+3*i-2)=2*P(3,1)*Z;
    HC_P_Y(11,2*N+3*i-1)=2*P(3,2)*Z;
    HC_P_Y(11,2*N+3*i)=4*P(3,3)*Z-2+2*P(3,1)*X+2*P(3,2)*Y+2*P(3,4);
   
    %p34
    HC_P_Y(12,2*i-1)=0;
    HC_P_Y(12,2*i)=0;
    
    HC_P_Y(12,2*N+3*i-2)=2*P(3,1);
    HC_P_Y(12,2*N+3*i-1)=2*P(3,2);
    HC_P_Y(12,2*N+3*i)=2*P(3,3);
    
end;


 Cov_m_M(1:2*N,1:2*N)= Cov_m;
 Cov_m_M(2*N+1:5*N,2*N+1:5*N)= Cov_M;

 
alphaF=-inv(HC_P)*HC_P_Y;


CovP=alphaF*Cov_m_M*alphaF';
%AlphaF=
disp('');

