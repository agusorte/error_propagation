%deshomogenizar
%Agustin Ortega
%oct 2012

function m=deshomogenizar(m)

for i=1:size(m,2),
    m(:,i)=m(:,i)/m(3,i);
end;
