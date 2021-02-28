%Begin function definition
function [vesselness] = filter3D(image, s, ps)

%Define a meshgrid of x and y points within a range of -6s to 6s
%Going in steps of ps (pixel size)
[x, y, z] = meshgrid(-6*s:ps:6*s);
%Define the Gaussian kernel in D dimensions (Equation 3 Frangi 1998)
%gauss = (1./(sqrt(2.*pi.*s.^2)).^3).*exp(-(x.^2 + y.^2 + z.^2)./2.*(s.^2));

%Set derivatives of all Gauss
%[xx, xy, xz,
% yx, yy, yz,
% zx, zy, zz]

%xx, xy, xz
gauss_xx = (2.^(1./2).*s.^4.*x.^2.*exp(-s.^2.*(x.^2./2 + y.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2)) - (2.^(1./2).*s.^2.*exp(-s.^2.*(x.^2./2 + y.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2));
gauss_xy = (2.^(1./2).*x.*y.*exp(-(x.^2./2 + y.^2./2)./s.^2))./(4.*s.^4.*pi.^(3./2).*(s.^2).^(3./2));
gauss_xz = (2.^(1./2).*s.^4.*x.*z.*exp(-s.^2.*(x.^.2/2 + y.^.2/2 + z.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2));

%yy, yz
gauss_yy = (2.^(1./2).*s.^4.*y.^2.*exp(-s.^2.*(x.^2./2 + y.^2./2 + z.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2)) - (2.^(1./2).*s.^2.*exp(-s.^2.*(x.^2./2 + y.^2./2 + z.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2));
gauss_yz = (2.^(1./2).*s.^4.*y.*z.*exp(-s.^2.*(x.^2./2 + y.^2./2 + z.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2));

%zz
gauss_zz = (2.^(1./2).*s.^4.*z.^2.*exp(-s.^2.*(x.^2./2 + y.^2./2 + z.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2)) - (2.^(1./2).*s.^2.*exp(-s.^2.*(x.^2./2 + y.^2./2 + z.^2./2)))./(4.*pi.^(3./2).*(s.^2).^(3./2));

%Clear x y z values
clear x y z

%Find Gaussian derivatives of the image
gamma = 1; %Frangi 1998, just under equation (3)

%Note: Using dxy == dyx, dyz == dzy, etc
%Images xx, xy, xz
image_xx = convn(s^gamma*image, gauss_xx, 'same');
image_xy = convn(s^gamma*image, gauss_xy, 'same');
image_xz = convn(s^gamma*image, gauss_xz, 'same');

%Images yy & yz
image_yy = convn(s^gamma*image, gauss_yy, 'same');
image_yz = convn(s^gamma*image, gauss_yz, 'same');

%Image zz
image_zz = convn(s^gamma*image, gauss_zz, 'same');

%Conserve memory, clear Gauss derivatives and gamma
clear gauss_xx gauss_xy gauss_xz gauss_yy gauss_yz gauss_zz gamma

%Store eigenvalues within array
% Sort eigenvalues according to the paper
eVals1 = image_xx./3 + image_yy./3 + image_zz./3 + (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9)./(((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3) + (((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3);
eVals2 = image_xx./3 + image_yy./3 + image_zz./3 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9)./(2.*(((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3)) - (3.^(1./2).*((image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9)./(((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3) - (((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3)).*1i)./2 - (((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3)./2;
eVals3 = image_xx./3 + image_yy./3 + image_zz./3 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9)./(2.*(((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3)) + (3.^(1./2).*((image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9)./(((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3) - (((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3)).*1i)./2 - (((((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx.*image_yz.^2)./2 + (image_xz.^2.*image_yy)./2 + (image_xy.^2.*image_zz)./2 - (image_xx + image_yy + image_zz).^3./27 - image_xy.*image_xz.*image_yz - (image_xx.*image_yy.*image_zz)./2).^2 - (image_xy.^2./3 - (image_xx.*image_zz)./3 - (image_yy.*image_zz)./3 - (image_xx.*image_yy)./3 + image_xz.^2./3 + image_yz.^2./3 + (image_xx + image_yy + image_zz).^2./9).^3).^(1./2) - (image_xx.*image_yz.^2)./2 - (image_xz.^2.*image_yy)./2 - (image_xy.^2.*image_zz)./2 - ((image_xx + image_yy + image_zz).*(- image_xy.^2 - image_xz.^2 - image_yz.^2 + image_xx.*image_yy + image_xx.*image_zz + image_yy.*image_zz))./6 + (image_xx + image_yy + image_zz).^3./27 + image_xy.*image_xz.*image_yz + (image_xx.*image_yy.*image_zz)./2).^(1./3)./2;

%Conserve memory, clear image convolutions
clear image_xx image_xy image_xz image_yy image_yz image_zz

%Sort values according to Frangi 1998
%Form logical array
truthArray = eVals1 > eVals2;

%Make zeros arrays to hold new data
newArray1 = zeros(size(eVals1));
newArray2 = zeros(size(eVals2));

%Sort values based on the truthArray
newArray1(truthArray==1) = eVals2(truthArray==1);
newArray2(truthArray==1) = eVals1(truthArray==1);

%Set eigenvalues to real values only
eigs1 = real(newArray1);
eigs2 = real(newArray2);
eigs3 = real(eVals3);

%Clear eVals1/2/3
clear eVals1 eVals2 eVals3 newArray1 newArray2 truthArray

%Compute equation 12, returning the 'second order structureness'
S = sqrt(eigs1.^2 + eigs2.^2 + eigs3.^2);

%Set constants to use for equation 13
alpha = 0.5;
beta = 0.5;
c = max(S(:))/2;
Rb = abs(eigs1)./(sqrt(abs(eigs2.*eigs3))); 
Ra = abs(eigs2)./abs(eigs3);

%Compute equation 13
vesselness = (1 - exp(-((Ra.^2)./2.*alpha.^2)).*(exp(-(Rb.^2)./(2.*beta.^2))).*(1-exp(-(S.^2)./(2.*c.^2))));
vesselness(eigs2 > 0) = 0;
vesselness(eigs3 > 0) = 0;
vesselness(~isfinite(vesselness)) = 0;
end















