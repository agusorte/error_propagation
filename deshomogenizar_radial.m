%deshomogenizar with radial distortion
%Agustin Ortega
%oct 2012

function m=deshomogenizar_radial(m,r)

for i=1:size(m,2),
    m(:,i)=m(:,i)/(1+r*(m(1,i)^2+m(2,i)^2));
end;
