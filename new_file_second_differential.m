%Start of file
%This file is designed to convolute a generated 2D array with the second
%differentials of the Gaussian function
%The Gaussians used will be 2D, depending upon x & y
%15/10/2020 @SpencerAndTheMatt

%Define symbols for differentiation
syms x y s D real %real not a variable, makes syms real

%Define Gauss equation
gauss = (1/(sqrt(pi.*s.^2)).^D).*exp(-((x.^2 + y.^2))/2.*s.^2);

%Compute first differentials of gauss
gauss_dx = diff(gauss, x);
gauss_dy = diff(gauss, y);

%Compute second differentials
gauss_dx2 = diff(gauss_dx, x);
gauss_dy2 = diff(gauss_dy, y);
gauss_dxdy = diff(gauss_dx, y);

%Aside
%Plot 2d second differential of gauss, with respect to y
% figure('Name', 'Two dimensional Gauss second derivative')
% points = linspace(-5, 5, 1000);
% plot(points, gauss_2d(points, 1))

%Define meshgrid
X = linspace(-3, 3, 100);
Y = X;
[x y] = meshgrid(X, Y); %2 2D arrays of 1s in the centre column for x, and centre row for y

%Meshgrid defined, define second differentials
s = 1; %All differentials depend on s
D = 2; %D is 2 as the graph is 2 dimensional 
dx2 = (s.^4.*x.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2)) - (s.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2));
dy2 = (s.^4.*y.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2)) - (s.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2));
dxdy = (s.^4.*x.*y.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2));

%Show mesh
%Mesh function is mesh(x 2D array, y 2D array, z 2D array)
% figure('Name', 'dy2 graph');
% mesh(x, y, dy2);
% figure('Name', 'dx2 graph');
% mesh(x, y, dx2);
% figure('Name', 'dxdy graph');
% mesh(x, y, dxdy);

%Show surf
%Surf function is just a surface plot, but takes the same arguments as mesh
% figure('Name', 'dy2 graph, surf')
% surf(x, y, dy2)
% figure('Name', 'dx2 graph, surf')
% surf(x, y, dx2)
% figure('Name', 'dxdy graph, surf')
% surf(x, y, dxdy)

%Make test array to model a blood vessel
%i.e zeros everywhere with a line of 1s down the middle
testArray = zeros(100);
testArray(:,50) = 1;

%Convolute test array with dx2, dy2, dxdy
figure('Name', 'Test array convoluted with dy2')
mesh(conv2(testArray, dy2))
colorbar

figure('Name', 'Test array convoluted with dx2')
mesh(conv2(testArray, dx2))
colorbar

figure('Name', 'Test array convoluted with dxdy')
mesh(conv2(testArray, dxdy))
colorbar



% %Define function for 2D plot
% function y = gauss_2d(x, s)
%     y = (s.^4.*x.^2.*exp(-s.^2.*(x.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2)) - (s.^2.*exp(-s.^2.*(x.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2));
% end



