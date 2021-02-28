profile on
%First test
%Make zeros
image = zeros(20, 20, 20);

%Make blood vessel, horizontal, in centre of image
%Width = 3
image(9:11, 9:11, :) = 1;

%Call filter
vessel = filter3D(image, 3, 1);

figure
pcolor(image(:, :, 10));
figure
pcolor(vessel(:, :, 10));

%Second test
%Make zeros
image2 = zeros(50, 50, 50);

%Make two blood vessels, one horizontal, one vertical
image2(:, 24:26, 24:26) = 1;
image2(24:26, :, 24:26) = 1;

%Call filter
vessel2 = filter3D(image2, 3, 1);

%Plot
figure('name', 'blood vessel 2');
stem3(image2(:, :, :));
figure('name', 'filtered blood vessel 2');
stem3(vessel2(:, :, :));

profile viewer

%Compute dice score using binary values
diceImage = de2bi(image);
diceVessel = de2bi(vessel);
similarity = dice(diceImage, diceVessel)








