%convert vector to skew matrix
%Agustin Ortega
%oct 2012

function skew=vec2skew(vec)

skew=[0    -vec(3) vec(2);
      vec(3) 0   -vec(1);
     -vec(2) vec(1) 0];