%#Balkis Saidi (124785)
%#Ahsan Qadeer (123892)


function Assignment6

  img1 = double(imread('left.png'));
  img2 = double(imread('right.png'));

  %# parameters
  r = 5; %# Radio
  mind = 10; %# Minimum disparity
  maxd = 50; %# Maximum disparity
  best_c= 0; %Best correlation

  %# disparity matrix

  [heigh,width] = size(img1);

  disparity_map = zeros(heigh,width);

  % mean and mean_square

  % Left image

  for i = r+1 : heigh - r
     for j = r+1 : width - r
        mean_1(i,j) = mean2(img1(i-r:i+r,j-r:j+r));
        mean_square1(i,j) = mean2(img1(i-r:i+r,j-r:j+r).*img1(i-r:i+r,j-r:j+r));
     end
  end

  % Right image
  for i = r+1 : heigh - r
     for j = r+1 : width - r
        mean_2(i,j) = mean2(img2(i-r:i+r,j-r:j+r));
        mean_square2(i,j) = mean2(img2(i-r:i+r,j-r:j+r).*img2(i-r:i+r,j-r:j+r));
     end
  end

  for i = r+1 : heigh - r
    for j = r+1 : width - r
      best_c = 0;
      variance1 = mean_square1(i,j) - mean_1(i,j)^2; %# NCC is not defined for homogeneous image areas (variance is zero)
      x = img1(i-r:i+r,j-r:j+r);
      if (variance1 > 0)
        for m = min(j+mind, width-r): -1 : max(j - maxd, r+1)
            y = img2(i-r:i+r,m-r:m+r);
            variance2 = mean_square2(i,m) - mean_2(i,m)^2;
            if (variance2 > 0)
              covariancex_y = (mean2(x.*y) - mean_1(i,j)*mean_2(i,m));
              ncc = covariancex_y / sqrt(variance1 * variance2);
              if ncc > best_c
                best_c = ncc;
                disparity_map(i,j) = abs(m - j);
              end
            end
        end
      end
    end
  end


  %# Plot disparity map

  figure();
  disparity_map = medfilt2(disparity_map);
  imshow(disparity_map, []);

end
