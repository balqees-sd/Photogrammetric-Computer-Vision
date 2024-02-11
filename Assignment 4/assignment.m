function assignment4

   img1 = imread('img1.jpg');
   img2 = imread('img2.jpg');

 %picking 8 homologous points in the image pair

   figure, imshow(img1), hold on
   [x1, y1] = ginput(8);
   P1 = [x1, y1, ones(8,1)]
   hold on;
   plot(x1,y1,'y+');

   figure, imshow(img2), hold on
   [x2, y2] = ginput(8);
   P2 = [x2,y2,ones(8,1)]
   hold on;
   plot(x2,y2,'y+');


 %Conditioning

   t1=condition2(P1);
   c1 = t1*P1'

   t2=condition2(P2);
   c2 = t2*P2'

 %Af=0

   A = design_homo2(c1,c2)

   [U, D, V] = svd(A);

   f = reshape(V(:,9),3,3)

 %Reverse conditioning

   F = t2' * f * t1

 %singularity constraint

   F = singularity_constraint(F)

 %epipolar lines

   L1 = F'* P2'
   L2 = F * P1'

   figure, imshow(img1), hold on;
   for i=1:8
     hline(L1(:,i));
   endfor

   figure, imshow(img2), hold on;
   for i=1:8
     hline(L2(:,i));
   endfor

 %geometric image error (sampson distance)

 a = sum ((P2'.* L2).^2)
 b = L2(1,:).^2+L2(2,:).^2 + L1(1,:).^2+L1(2,:).^2
 err = sum(a./b)/8

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

  [U1,D1,V1]= svd(F)
  D1(3,3)=0
  F = U1*D1*V1'

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




