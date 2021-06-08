close
clear all
clc
original_image = imread('boop.jpg');%Colormap is Empty
figure
imshow(original_image);
original_image = rgb2gray(original_image);
mask_size = 3; %Size of Gaussian Mask
sigma = 2; %This value was decided upon arbitrarily - decides how smooth
%the image should be
%Gaussian Mask
%The scaling is done by first generating all Gaussian coefficients, then
%summing up all the coefficients, and dividing the element by this sum -
%this ensures that the area underneath the gaussian is always 1;
index = -floor(mask_size/2) : floor(mask_size/2);%Rounding the NxN gaussian mask's range down
%toward infinity
[X Y] = meshgrid(index, index); %This function generates 2d coordinates
%The important aspect of this function is that it creates spatial
%coordinates that will be the same size as the mask desired. For odd values
%of N, the spatial coordinates will be symmetric.
% Gaussian is as below:
% 
% $$G(x,y) = \frac{1}{2 \pi \sigma ^{2}} e^{- \frac{x^{2} + y^{2}}{2 \sigma ^{2}}}$$
% 
%Mask below - 
generated_mask = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
%Normalizing the generated mask, so the total area is 1
generated_mask_normalized = generated_mask / sum(generated_mask(:));
generated_mask_column = generated_mask_normalized(:);%Method of creating a column vector
%Filtering should start here
%The image is converted to class type double, for higher precision
original_image_double_precision = im2double(original_image);
%We need to pad the array in order to properly filter the image, if the 
%border of the image is not padded, the spatial mask will not filter pixels
%beyond the image pixels.
padded_array = padarray(original_image_double_precision, [floor(mask_size/2) floor(mask_size/2)]);
%im2col will take distinct pixels around one pixel - the size of the
%gaussian mask since we have a column vector - converting the surrounding
%pixels into columns, and creating a matrix of pixel groups NxN
columns = im2col(padded_array, [mask_size mask_size], 'sliding');
%We have come to the point where we must sum over the rows for each column,
%this multiplication is done by a built in function bsxfun
pixel_group_filter = sum(bsxfun(@times, columns, generated_mask_column), 1);
%Converting all pixel groups back to their respective positions after
%applying the filter
output_image = col2im(pixel_group_filter, [mask_size mask_size], size(padded_array), 'sliding');
%# Filter it - The Filter can also being created using built in matlab functions
G = fspecial('gaussian',[3 3],2);
%Two ways to apply filter  - Method one using imfilter
image_meth1 = imfilter(original_image,G,'same');
figure
imshow(image_meth1)
title(['Gaussian Filter with \sigma = 2']);
%Method 2 - convolving
convolution = conv2(original_image,G, 'same');
imshow(convolution, []);
title('Gaussian Filter with \sigma = 2');
image_meth2 = convolution;
%%Calculating Gradient - Method 1
[Gx,Gy] = gradient(G);   
%%Method 2
%This can be calculated using the first order derivative of the gaussian
%The kernal must be the same size 
%https://www.mathworks.com/matlabcentral/fileexchange/8060-gradient-using-first-order-derivative-of-gaussian
epsilon=1e-2;
halfsize=ceil(sigma*sqrt(-2*log(sqrt(2*pi)*sigma*epsilon)));
size=2*halfsize+1;
for i=1:size
    for j=1:size
        u=[i-halfsize-1 j-halfsize-1];
        hx(i,j)=normrnd(u(1),sigma)*normrnd(u(2),sigma);
    end
end
hx=hx/sqrt(sum(sum(abs(hx).*abs(hx))));
hy=hx';

X = conv2(image_meth2, hx, 'same');
Y = conv2(image_meth2, hy, 'same');
wiki_eq_gradient = sqrt(Gy^2*Gx^2);
wiki_eq_direction = atan2 (X, Y);
wiki_eq_direction = wiki_eq_direction*180/pi;
%%Finding Direction
%Edge direction angle is rounded to one of four angles, for each direction
%0 45 90 135 as said in the wikipedia page. Edge detection within these
%regions are floored to zero
%Doing this for a max of a 1024 by 1024 image
directions = zeros(1024, 1024);
%Floor value to 0 if inbetween the 4 direction values
for i = 1  : 183
    for j = 1 : 276
        if ((wiki_eq_direction(i, j) > 0 ) && (wiki_eq_direction(i, j) < 22.5) || (wiki_eq_direction(i, j) > 157.5) && (wiki_eq_direction(i, j) < -157.5))
            directions(i, j) = 0;
        end
        
        if ((wiki_eq_direction(i, j) > 22.5) && (wiki_eq_direction(i, j) < 67.5) || (wiki_eq_direction(i, j) < -112.5) && (wiki_eq_direction(i, j) > -157.5))
            directions(i, j) = 45;
        end
        
        if ((wiki_eq_direction(i, j) > 67.5 && wiki_eq_direction(i, j) < 112.5) || (wiki_eq_direction(i, j) < -67.5 && wiki_eq_direction(i, j) > 112.5))
            directions(i, j) = 90;
        end
        
        if ((wiki_eq_direction(i, j) > 112.5 && wiki_eq_direction(i, j) <= 157.5) || (wiki_eq_direction(i, j) < -22.5 && wiki_eq_direction(i, j) > -67.5))
            directions(i, j) = 135;
        end
    end
end
%% Non Maximum Surpression Edge thinning