%Group 04
%Balkis Saidi (124785)
%Ahsan Qadeer (123892)
%Final Project
%Essential Matrix Estimation and Non-linear Optimization

function ProjectPCV
pkg load optim
%Task1
%a)compute the calibration matrices K1 and K2

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

  %Calibration matrix K1 and K2

  K1=[1,0,0;
      0,1,0;
      0,0,1];
 %compute calibration matrix K2

P1 = P1/P1(3,4);

M = P1(:,1:3);

m3 = M(3,:);

if det(M) > 0
  lambda=1/norm(m3)
else
  lambda=-1/norm(m3)
endif

M=lambda*M;

[Q, R] = qr (inv(M));

R=inv(R)

K2 = R

%b)Estimate the essential matrix

E = K2'*F*K1;
E = singularity_constraint2(E)

%c)Resolve the fourfold ambiguity of the essential matrix

fh = fopen('pp.dat', 'r');
  A = fscanf(fh, '%f%f%f%f%f%f%f', [7 inf]);
  fclose(fh);
  x1 = A(1:2, :);
  x1(3, :)=1;
  x2 = A(3:4, :);
  x2(3, :)=1;
  x3 = A(5:7, :);
  x3(4, :)=1;

%the four pose configurations can be computed from E.

[U,D,V]= svd(E);

  W=[0,-1,0;
      1,0,0;
      0,0,1];

R1=U*W*V'; C1 = U(:,3);
R2=U*W*V'; C2 = -U(:,3);
R3=U*W'*V'; C3 = U(:,3);
R4=U*W'*V'; C4 = -U(:,3);

%note that the  det(R)=1. If det(R)=−1, the camera pose must be corrected i.e. C=−C and R=−R.

if det(R1) == -1
    R1 = -R1;
    C1 = -C1;
  endif

 if det(R2) == -1
    R2 = -R2;
    C2 = -C2;
  endif

if det(R3) == -1
    R3 = -R3;
    C3 = -C3;
  endif

if det(R4) == -1
    R4 = -R4;
    C4 = -C4;
  endif

%the camera pose can be written as: P=KR[I3×3 −C]

%Case 1
PN = K1*R1*[eye(3,3) -C1];
P1 = K2*R1*[eye(3,3) -C1];

X= triangulation(PN,P1,x1,x2)

%Case 2
PN = K1*R2*[eye(3,3) -C2];
P1 = K2*R2*[eye(3,3) -C2];

X= triangulation(PN,P1,x1,x2)

%Case 3
PN = K1*R3*[eye(3,3) -C3];
P1 = K2*R3*[eye(3,3) -C3];

X= triangulation(PN,P1,x1,x2)

%Case 4
PN = K1*R4*[eye(3,3) -C4];
P1 = K2*R4*[eye(3,3) -C4];

X= triangulation(PN,P1,x1,x2)

%Check the depth (z) of the 3D points computed, it should be positive for the plausible combination
%we have to choose either 1 or 3
%skew-symmetric matrix S associated with C
 S =[0,-C1(3),C1(2);
      C1(3),0,-C1(1);
      -C1(2),C1(1),0];
%New essential matrix
        EN = S*R1;

%d)Compute and plot the epipolar lines from the essential matrix
%New fundamental matrix
FN = inv(K2')*EN*inv(K1)

 %epipolar lines

   L1 = FN' * x2
   L2 = FN * x1

   figure
   plot (x1(1,:),x1(2,:),"x");
   hold on
   for i=1:5
     hline(L1(:,i));
   endfor

   figure,
   plot (x2(1,:),x2(2,:),"x");
   hold on
   for i=1:5
     hline(L2(:,i));
   endfor

%Task2
%a)the geometric error
 a = sum (x2.* L2).^2
 b = L2(1,:).^2+L2(2,:).^2 + L1(1,:).^2+L1(2,:).^2
 err = mean(a./b)

%b) Perform a non-linear optimization by means of the indirect optimization using a built-in function of your choice (Levenberg-Marquart, etc.)

[f, p, cvg, iter]  = leasqr(x1, x2, FN, @(F, x1, x2) symmetric_epipolar_distance(F, x1, x2))

% Non-linear optimization using Levenberg-Marquardt with Trust Region
%options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'TolX', 1e-6, 'TolFun', 1e-6);
%F_opt = lsqnonlin(@(F) symmetric_epipolar_distance(F, x1, x2), FN, [], [], options);

%d) Re-calculate the geometric error. Report and comment on the results.
%The re-calculated geometric error evaluate the accuracy of the optimized solution
endfunction


% Define the function to be optimized
function y = symmetric_epipolar_distance(F, x1, x2)
    L1 = F' * x2
    L2 = F * x1
    y = mean((sum (x2.* L2).^2)./(L2(1,:).^2+L2(2,:).^2 + L1(1,:).^2+L1(2,:).^2));
endfunction

function  h = hline(l, varargin)
%        ==================
	if abs(l(1)) < abs(l(2))                                  % More horizontal
		xlim = get(get(gcf, 'CurrentAxes'), 'Xlim');
		x1 = cross(l, [1; 0; -xlim(1)]);
		x2 = cross(l, [1; 0; -xlim(2)]);
	else                                                        % More vertical
		ylim = get(get(gcf, 'CurrentAxes'), 'Ylim');
		x1 = cross(l, [0; 1; -ylim(1)]);
		x2 = cross(l, [0; 1; -ylim(2)]);
	end
	x1 = x1 / x1(3);
	x2 = x2 / x2(3);
	h = line([x1(1) x2(1)], [x1(2) x2(2)], varargin{:});
end


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

function E = singularity_constraint2(E)

  [U1,D1,V1]= svd(E);
  d1 = D1(1,1);
  d2 = D1(2,2);
  D1(1,1) = (d1+d2)/2;
  D1(2,2) = (d1+d2)/2;
  D1(3,3)=0;
  E = U1*D1*V1';

endfunction
