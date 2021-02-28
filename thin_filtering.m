% thin_filtering.
% Created on 11/12/2020.
% By @SpencerAndTheMatt.
% This programme is designed use the method of eigenvalue decompositiion 
% to filter out non-vessel like structures from an image, in a similar
% manner to Frangi et al 1998.
% This is a re-write of thin_vessel_filtering.m.
% This file uses the files provided by L. Beltrachini as an aid. 


%Begin function definition
function [vesselness] = thin_filtering(image, s, ps) 
% ------------------------------------------------------
% This function takes arguments of image, s and ps and 
% returns the image, filtered from objects with a low
% measure of 'vesselness'. The aim is to obtain a clear
% picture of only vessels
% 
% Arguments:
%     image: The image to be filtered. Default value is
%            image with vessel down centre and 15 by 15
%            blob in corner.
%              
%     s: scale factor. Default value is 1.
% 
%     ps: pixel size. Default value is 0.1.
% ------------------------------------------------------

%Use nargin with switch statement for default values
switch nargin
    %No input arguments
    case 0
        %Form image and make a blood vessel and blob
        image = zeros(500);
        image(:,[248:252]) = 1;
        for index = 1:15;
            for jndex = 1:15;
                image(index, jndex) = 1;
            end
        end
        
        %Set s and ps
        s = 1;
        ps = 0.1;
        
    %1 input argument
    case 1
        s = 1;
        ps = 0.1;
    
    %2 input arguments
    case 2
        ps = 0.1;
end

%Define meshgrid
[x, y] = meshgrid(-6*s:ps:6*s);
    
%Define Gaussian Kernel, equation (3)
gauss = exp(-(x.^2/2 + y.^2/2)/s^2);

%Find derivatives of Gauss kernel
gauss_dxx = (x.^2.*gauss)/(2*pi*s^6) - gauss/(2*pi*s^4);
gauss_dyy = - gauss / (2*pi*s^4) + y.^2.*gauss/(2*pi*s^6);
gauss_dxy = x.*y.*gauss/(2*pi*s^6);

%Find Gaussian derivatives of the image
gamma = 1; %Frangi 1998, just under equation (3)
image_xx = conv2(s^gamma*image, gauss_dxx, 'same');
image_yy = conv2(s^gamma*image, gauss_dyy, 'same');
image_xy = conv2(s^gamma*image, gauss_dxy, 'same');

%Find eigenvalues
e_vals(: ,: ,1) = image_xx/2 + image_yy/2 - (image_xx.^2 - 2*image_xx.*image_yy + 4*image_xy.^2 + image_yy.^2).^(1/2)/2;
e_vals(: ,: ,2) = image_xx/2 + image_yy/2 + (image_xx.^2 - 2*image_xx.*image_yy + 4*image_xy.^2 + image_yy.^2).^(1/2)/2;

%Sort eigenvalues based on Frangi 1998 paper
[~, I] = sort(abs(e_vals), 3, 'ascend');
eigs1 = (I(: ,: ,1) == 1).*e_vals(: ,: ,1) + (I(: ,: ,1) == 2).*e_vals(: ,: ,2);
eigs2 = (I(: ,: ,2) == 1).*e_vals(: ,: ,1) + (I(: ,: ,2) == 2).*e_vals(:, : ,2);

%Compute equation 12
S = sqrt(eigs1.^2 + eigs2.^2);

%Obtain result of equation 15
beta = 0.5;
c = max(S(:))/2;
Rb = eigs1./eigs2;
vesselness = exp(-Rb.^2/2/beta^2).*(1 - exp(-S.^2/2/c^2));
vesselness(eigs2 > 0) = 0;

end


    

























