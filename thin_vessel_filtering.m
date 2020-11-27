% New File
% SpencerAndTheMatt 08/11/2020
% This file is a re-write of new_file_second_differential
% This file will differentiate a gauss function twice, partially, with respect
% to x and y and xy.
% The differentiated Gauss' will then be plotted
% I will then make a 'blood vessel' out of zeros and ones, convolute it, 
% and then attempt to form a Hessian matrix

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

%Define meshgrid
X = linspace(-1, 1, 10);
Y = X;
[x y] = meshgrid(X, Y);

%Define second differentials non-symbolicly to form an array
s = 5; %Width of vessel to be detected
D = 2; %D is 2 as the graph is 2 dimensional 
dx2 = (s.^4.*x.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2)) - (s.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2));
dy2 = (s.^4.*y.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2)) - (s.^2.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2));
dxdy = (s.^4.*x.*y.*exp(-s.^2.*(x.^2/2 + y.^2/2)))/(pi.^(1/2).*(s.^2).^(1/2));

%Show surf
%Surf function is just a surface plot, but takes the same arguments as mesh
figure('Name', 'dy2 graph, surf')
surf(x, y, dy2)
figure('Name', 'dx2 graph, surf')
surf(x, y, dx2)
figure('Name', 'dxdy graph, surf')
surf(x, y, dxdy)

%Make fake blood vessel
%2D array of zeros with 1s in the centre column 
fakeBloodVessel = zeros(50);
for index = 23:27;
    fakeBloodVessel(:,index) = 1;
end;

%Add blob to fake blood vessel, in top left corner
for index = 1:6
    for jndex = 1:6
        fakeBloodVessel(index, jndex) = 1;
    end
end

%Convolute dx2, dy2, dxdy, with fakeBloodVessel
% figure('Name', 'dx2 convoluted with fakeBloodVessel')
convdx2 = conv2(dx2, fakeBloodVessel);
% mesh(convdx2)
% colorbar

% figure('Name', 'dy2 convoluted with fakeBloodVessel')
convdy2 = conv2(dy2, fakeBloodVessel);
% mesh(convdy2)
% colorbar

% figure('Name', 'dxdy convoluted with fakeBloodVessel')
convdxdy = conv2(dxdy, fakeBloodVessel);
% mesh(convdxdy)
% colorbar

%Setting up Hessian Matrix
syms Hxx Hyy Hxy real
H = [Hxx, Hxy; Hxy, Hyy];

%Find eigenvalues in symbolic form
eigenvalues_H = eig(H);

%Substitute in x and y arrays into eigenvalues_H
eigenvalues_1 = convdx2./2 + convdy2./2 - convdx2.^2 - 2.*convdx2.*convdy2 + 4.*convdxdy.^2 + convdy2.^2.^(1./2)./2;
eigenvalues_2 = convdx2./2 + convdy2./2 + convdx2.^2 - 2.*convdx2.*convdy2 + 4.*convdxdy.^2 + convdy2.^2.^(1./2)./2;

%Switch round eigenvalues
%The absolute value of lambda 1 << absolute value of lambda 2
for i = 1:numel(eigenvalues_1);
    if abs(eigenvalues_1(i)) < abs(eigenvalues_2(i));
        temp = eigenvalues_1(i);
        eigenvalues_1(i) = eigenvalues_2(i);
        eigenvalues_2(i) = temp;
    end;
end;

%Plotting eigenvalue matrices
figure('Name', 'Contour, eigenvalues_1')
pcolor(real(eigenvalues_1))
colorbar

figure('Name', 'Contour, eigenvalues_2')
contourf(real(eigenvalues_2))
colorbar

%Plotting a quiver plot, out of curiousity more than anything
figure('Name', 'Quiver plot')
quiver(real(eigenvalues_1), real(eigenvalues_2))


%Apply Rb equation
% number = 0.8; %Measure of vesselness
% for index = 1:numel(eigenvalues_1);
%     if eigenvalues_1(index)/eigenvalues_2(index) > number;
%         eigenvalues_1(index) = 0;
%         eigenvalues_2(index) = 0;
%     end
% end

%Plot eigenvalues with Rb equation applied
% figure('Name', 'Contour, eigenvalues_1 filtered')
% contourf(real(eigenvalues_1))
% colorbar
%
% figure('Name', 'Contour, eigenvalues_2 filtered')
% contourf(real(eigenvalues_2))
% colorbar

%Use equation 15 from Frangi paper
tempArray = zeros(size(eigenvalues_1));
beta = 1; %S, c, beta are 'thresholds which control the sensitivity
c = 1; %of the line filter to the measures Ra, Rb and S' [Frangi paper, 
S = 1; %just under equation (13)].
for index = 1:numel(tempArray);
    if eigenvalues_2(index) > 0;
        tempArray(index) = 0;
    else;
        tempArray(index) = exp((-(eigenvalues_1(index)/eigenvalues_2(index)).^2)/(2.*beta.^2)) .* (1 - exp((-(S).^2)./(2 .* c .^ 2)));
    end;
end;

%Plot equation 15
figure('Name', 'Equation 15 plot')
contourf(tempArray)
colorbar

















