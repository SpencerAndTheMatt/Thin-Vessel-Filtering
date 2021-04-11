%Begin function definition
function [vesselness] = filter3D(image, s, ps)

%Define a meshgrid of x and y points within a range of -6s to 6s
%Going in steps of ps (pixel size)
[x, y, z] = meshgrid(-6*s:ps:6*s);
%Define the Gaussian kernel in D dimensions (Equation 3 Frangi 1998)
%gauss = (1./(sqrt(2.*pi.*s.^2)).^3).*exp(-(x.^2 + y.^2 + z.^2)./2.*(s.^2));

%Set derivatives of all Gauss
%[xx, xy, xz,
% yx, yy, yz
% zx, zy, zz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Previously I used xy = yx, however this led the filter
% to be ineffective. I'm not entirely sure why, however the filter
% seems to work with xy /= yx, hence all 9 derivatives are included.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%xx, xy, xz
gauss_xx = (x.^2.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./(s.^4) - exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2)./(s.^2);
gauss_xy = (x.*y.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4;
gauss_xz = (x.*z.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4;

%yx, yy, yz
gauss_yx = (x.*y.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4;
gauss_yy = (y.^2.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4 - exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2)./s.^2;
gauss_yz = (y.*z.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4;
 
%zx, zy, zz
gauss_zx = (x.*z.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4;
gauss_zy = (y.*z.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4;
gauss_zz = (z.^2.*exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2))./s.^4 - exp(-(x.^2./2 + y.^2./2 + z.^2./2)./s.^2)./s.^2;
 
%Clear x y z values
clear x y z 

%Find Gaussian derivatives of the image
gamma = 1; %Frangi 1998, just under equation (3)

%Images xx, xy, xz
image_xx = convn(s^gamma*image, gauss_xx, 'same');
image_xy = convn(s^gamma*image, gauss_xy, 'same');
image_xz = convn(s^gamma*image, gauss_xz, 'same');

%Images yx, yy, yz
image_yx = convn(s^gamma*image, gauss_yx, 'same');
image_yy = convn(s^gamma*image, gauss_yy, 'same');
image_yz = convn(s^gamma*image, gauss_yz, 'same');

%Images zx, zy, zz
image_zx = convn(s^gamma*image, gauss_zx, 'same');
image_zy = convn(s^gamma*image, gauss_zy, 'same');
image_zz = convn(s^gamma*image, gauss_zz, 'same');

%Conserve memory, clear Gauss derivatives and gamma
clear gauss_xx gauss_xy gauss_xz gauss_yx gauss_yy gauss_yz gauss_zx gauss_zy gauss_zz gamma

%Store eigenvalues within arrays
eVals1 = image_xx./3 + image_yy./3 + image_zz./3 + ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9)./(((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3) + (((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3);
eVals2 = image_xx./3 + image_yy./3 + image_zz./3 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9)./(2.*(((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3)) - (3.^(1./2).*(((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9)./(((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3) - (((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3)).*1i)./2 - (((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3)./2;
eVals3 = image_xx./3 + image_yy./3 + image_zz./3 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9)./(2.*(((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3)) + (3.^(1./2).*(((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9)./(((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3) - (((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3)).*1i)./2 - (((((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 - (image_xx + image_yy + image_zz).^3./27 - (image_xx.*image_yy.*image_zz)./2 + (image_xx.*image_yz.*image_zy)./2 + (image_xy.*image_yx.*image_zz)./2 - (image_xy.*image_zx.*image_yz)./2 - (image_yx.*image_xz.*image_zy)./2 + (image_xz.*image_yy.*image_zx)./2).^2 - ((image_xy.*image_yx)./3 - (image_xx.*image_yy)./3 - (image_xx.*image_zz)./3 + (image_xz.*image_zx)./3 - (image_yy.*image_zz)./3 + (image_yz.*image_zy)./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - ((image_xx + image_yy + image_zz).*(image_xx.*image_yy - image_xy.*image_yx + image_xx.*image_zz - image_xz.*image_zx + image_yy.*image_zz - image_yz.*image_zy))./6 + (image_xx + image_yy + image_zz).^3./27 + (image_xx.*image_yy.*image_zz)./2 - (image_xx.*image_yz.*image_zy)./2 - (image_xy.*image_yx.*image_zz)./2 + (image_xy.*image_zx.*image_yz)./2 + (image_yx.*image_xz.*image_zy)./2 - (image_xz.*image_yy.*image_zx)./2).^(1./3)./2;
 
%Conserve memory, clear image convolutions
clear image_xx image_xy image_xz image_yx image_yy image_yz image_zx image_zy image_zz

%Arrange eigenvalues in ascending order and sort by absolute value
%Create array to hold eigenvalues, ev
ev = zeros(size(eVals1, 1), size(eVals1, 2), size(eVals1, 3), 3);

%Put eigenvalues in ev array
ev(:, :, :, 1) = eVals1;
ev(:, :, :, 2) = eVals2;
ev(:, :, :, 3) = eVals3;

%Sort eigenvalues using an index sort, based on the 4th dimension
[~,I]=sort(abs(ev),4, 'ascend');

%Sort eigenvalues using the index
eVals1=(I(:,:,:,1)==1).*ev(:,:,:,1) + (I(:, :, :,1) == 2).*ev(:,:,:,2) + (I(:,:,:,1)==3).*ev(:,:,:,3);
eVals2=(I(:,:,:,2)==1).*ev(:,:,:,1) + (I(:, :, :,2) == 2).*ev(:,:,:,2) + (I(:,:,:,2)==3).*ev(:,:,:,3);
eVals3=(I(:,:,:,3)==1).*ev(:,:,:,1) + (I(:, :, :,3) == 2).*ev(:,:,:,2) + (I(:,:,:,3)==3).*ev(:,:,:,3);

%Set eigenvalues to real values only
eigs1 = real(eVals1);
eigs2 = real(eVals2);
eigs3 = real(eVals3);

%Clear eVals1/2/3 and ev, I
clear eVals1 eVals2 eVals3 ev I

%Compute equation 12, returning the 'second order structureness'
S = sqrt(eigs1.^2 + eigs2.^2 + eigs3.^2);

%Set constants to use for equation 13
alpha = 0.5;
beta = 0.5;
c = max(S(:))/2;
Rb = abs(eigs1)./(sqrt(abs(eigs2.*eigs3))); 
Ra = abs(eigs2)./abs(eigs3);

%Set values for computing equation 13
Ra_exp = (1 - exp(-((Ra.^2)./(2*alpha.^2))));
Rb_exp = exp(-((Rb.^2)./(2.*beta.^2)));
S_exp = (1 - exp(-((S.^2)/(2.*c.^2))));

%Clear values used to calculate Ra_exp, Rb_exp, S_exp
clear alpha beta c Rb Ra 

%Compute equation 13
vesselness = Ra_exp.*Rb_exp.*S_exp;
vesselness(eigs2 > 0) = 0;
vesselness(eigs3 > 0) = 0;
vesselness(~isfinite(vesselness)) = 0;
vesselness(isnan(vesselness)) = 0;
end

