%ploting lines by sampling
%Agustin Ortega
%oct 2012

function plotting_line_sampled(N,mu,r_lambda)
%N number of points by line
% m lines

for i=1:2:size(mu,2),

    line=cross(mu(:,i),mu(:,i+1));
    x=mu(1,i):1/N:mu(1,i+1);
    
    
    y=-(line(3)+line(1)*x)/line(2);
    
    
    xd=x./(1+r_lambda*((x-320).^2+(y-240).^2))+x;
    yd=y./(1+r_lambda*((x-320).^2+(y-240).^2))+y;
    
    
    

    plot(xd,yd,'r.-');
    
end;