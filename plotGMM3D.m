function plotGMM3D(Mu, Sigma, col1, alpha1, w)

 %col color;
 %alpha transparencia
 %w scala w
nbData = size(Mu,2);
nbpoints = 10;

 
darkColor=col1-.3;
darkColor(find(darkColor<0)) = 0;

 
for n=1:nbData

 
  stdev = sqrtm(w.*Sigma(:,:,n));

 
  [x,y,z] = sphere(nbpoints);
  x2 = reshape(x,[size(x,1)*size(x,2) 1]);
  y2 = reshape(y,[size(y,1)*size(y,2) 1]);
  z2 = reshape(z,[size(z,1)*size(z,2) 1]);

 
  D = [x2 y2 z2] * real(stdev) + repmat(Mu(:,n)',size(x2,1),1);
  X = reshape(D(:,1), [size(x,1) size(x,2)]);
  Y = reshape(D(:,2), [size(y,1) size(y,2)]);
  Z = reshape(D(:,3), [size(z,1) size(z,2)]);

 
  h(1) = surface(X,Y,Z);

  
  set(h(1), 'FaceColor', col1, 'EdgeColor', darkColor, 'FaceAlpha', alpha1);

  
  set(h(1), 'FaceColor', col1, 'EdgeColor', 'none', 'FaceAlpha', alpha1);
  camlight('headlight');
  lighting('gouraud');
  material('metal');

 
  set(gca, 'AmbientLightColor', [0.75 0.75 0.75]);
  set(h(1), 'AmbientStrength', 0.7);
end