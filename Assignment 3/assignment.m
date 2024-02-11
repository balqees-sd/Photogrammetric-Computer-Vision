% Asignment 3

function Assignment()

img = imread('img3.jpeg');
imshow(img);

% Read the example points

points1 = [317.61 232.03;813.09 373.6; 723.43 939.87;836.69 1156.9; 1072.6 609.54; 1369.9 538.76];

points1 = [points1,ones(6,1)]

points2 = [0 270 282;0 196 83; 82 0 207; 223 0 250; 74 78 0; 238 160  0];

points2 = [points2,ones(6,1)]

T1=condition2(points1);
T2=condition3(points2);

c1 = T1*points1';
c2 = T2*points2';

A = design_homo2(c1,c2);

[U,D,V] = svd(A);

h = V(:,end);

H = inv(T1) * reshape(h,4,3)' * T2

%Projection matrix H

H = H/H(3,4);

M = H(:,1:3);

m3 = M(3,:);

if det(M) > 0
  lambda=1/norm(m3)
else
  lambda=-1/norm(m3)
endif

M=lambda*M;

[Q, R] = qr (inv(M));
% Rotation matrix
Q=inv(Q);
Q([1,3], :)= Q([1,3], :)*(-1)
% Calibration matrix
R=inv(R);
R(:, [1,3])= R(:, [1,3])*(-1)




[Uc,Dc,Vc] = svd(H);

C = Vc(:,end)';
C = C/C(1,4)

endfunction

% Code taken from solution 2
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
  tz= mean(points(:,3));
  sx= mean(abs(points(:,1) - tx));
  sy= mean(abs(points(:,2) - ty));
  sz= mean(abs(points(:,3) - tz));

  T = [1/sx 0 0 -tx/sx;
        0 1/sy 0 -ty/sy;
        0 0 1/sz -tz/sz;
        0  0 0 1];
endfunction

function A = design_homo2(c1,c2)
  A=[];
  for i=1 : size(c1,2)
    A= [A;
       -c1(3, i)*c2(:, i)' 0 0 0 0 c1(1, i)*c2(:, i)';
       0 0 0 0 -c1(3, i)*c2(:, i)' c1(2, i)*c2(:, i)'];
  endfor
endfunction

