%Gruppe 18
%Members:
%Balkis Saidi (124785)
%Ahsan Qadeer (123892)

function assignment5

fh = fopen('bh.dat', 'r');
A = fscanf(fh, '%f%f%f%f', [4 inf]);
fclose(fh);
x1 = A(1:2, :);
x1(3, :)=1;
x2 = A(3:4, :);
x2(3, :)=1;

 %Conditioning

   t1=condition2(x1');
   c1 = t1*x1;

   t2=condition2(x2');
   c2 = t2*x2;

  %Af=0

   A = design_homo2(c1,c2);

   [U, D, V] = svd(A);

   f = reshape(V(:,9),3,3)';

  %Reverse conditioning

   F = t2' * f * t1;

 %singularity constraint

   F = singularity_constraint(F);


  [PN,P1] = projective_matrices(F);

  X = triangulation(PN,P1,x1,x2);

  figure; scatter3(X(1, :), X(2, :), X(3, :), 10, 'filled');
  axis square; view(32, 75);


  % task 2

  fh = fopen('pp.dat', 'r');
  A = fscanf(fh, '%f%f%f%f%f%f%f', [7 inf]);
  fclose(fh);
  x1 = A(1:2, :);
  x1(3, :)=1;
  x2 = A(3:4, :);
  x3(3, :)=1;
  x3 = A(5:7, :);
  x3(4, :)=1;

  X2 = triangulation(PN,P1,x1,x2);

  H = homography3D(X2,x3);

  XF = H*X;

  XF = XF ./ XF(4,:);

  figure; scatter3(XF(1, :), XF(2, :), XF(3, :), 10, 'filled');
  axis square; view(32, 75);



endfunction

function H = homography3D(x1,x2)

   %Conditioning

   t1=condition3(x1');
   c1 = t1*x1;

   t2=condition3(x2');
   c2 = t2*x2;

   A = design_homo3(c1,c2);

   [U, D, V] = svd(A);

   H = reshape(V(:,16),4,4)';

   H = inv(t2) * H * t1;

  endfunction

function X= triangulation(PN,P1,x1,x2)

  for i=1 : size(x1,2)
  A=[x1(1,i)*PN(3,:)-PN(1,:);
  x1(2,i)*PN(3,:)-PN(2,:);
  x2(1,i)*P1(3,:)-P1(1,:);
  x2(2,i)*P1(3,:)-P1(2,:)];

  [U, D, V] = svd (A);

   X(:,i) = V(:,end);

   X(:,i) = X(:,i) ./ X(4,i);

  endfor

 endfunction

function [PN,P1] = projective_matrices(F)
  [U1,D1,V1]= svd(F);
  e2=U1(:,end);
  a=[0,-e2(3),e2(2);
     e2(3),0,-e2(1);
     -e2(2),e2(1),0];


  PN=[1,0,0,0;
      0,1,0,0;
      0,0,1,0];

   P1= [a*F+[e2,e2,e2],e2];


  endfunction

function T = condition2(points)

  tx= mean(points(:,1));
  ty= mean(points(:,2));
  sx= mean(abs(points(:,1) - tx));
  sy= mean(abs(points(:,2) - ty));

  T = [1/sx 0 -tx/sx;
        0 1/sy -ty/sy;
        0  0  1];
endfunction

function T = condition3(points)

  tx= mean(points(:,1));
  ty= mean(points(:,2));
  tz = mean(points(:,3));
  sx= mean(abs(points(:,1) - tx));
  sy= mean(abs(points(:,2) - ty));
  sz= mean(abs(points(:,3) - tz));

  T = [1/sx 0 0 -tx/sx;
        0 1/sy  0 -ty/sy;
        0 0 1/sz -tz/sz;
        0  0  0 1];
endfunction

function A = design_homo3(c1,c2)
  A=[];
  for i=1 : size(c1,2)
    A = [ A; -c2(4,i)*c1(:,i)' 0 0 0 0 0 0 0 0 c2(1,i)*c1(:,i)';
          0 0 0 0 -c2(4,i)*c1(:,i)' 0 0 0 0 c2(2,i)*c1(:,i)';
          0 0 0 0 0 0 0 0 -c2(4,i)*c1(:,i)' c2(3,i)*c1(:,i)'];
  endfor
endfunction

function A = design_homo2(c1,c2)
  A=[];
  for i=1 : size(c1,2)
    A= [A; c1(1, i)*c2(1, i) c1(2,i)*c2(1,i) c2(1,i) c1(1,i)*c2(2,i) c1(2,i)*c2(2,i) c2(2,i) c1(1,i) c1(2,i) 1];
  endfor
endfunction

function F = singularity_constraint(F)

  [U1,D1,V1]= svd(F);
  D1(3,3)=0;
  F = U1*D1*V1';

endfunction


