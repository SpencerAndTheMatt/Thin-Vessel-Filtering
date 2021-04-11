%test3dFilter8
%An extension of test 6, however this will involve multiple vessels and
%noise
image = zeros(20, 20, 20);
image(:, 9:11, 9:11) = 1;
image(9:11, 9:11, :) = 1;
image(9:11, :, 9:11) = 1;
imageCopy = image;
image(4:6, 4:6, 10) = 1;

noise = 0.5;
image = image + (image ~= 1).*noise.*rand(size(image));

%Apply multiscale filter
s=1:0.1:6; % values of s to test
ps=1; % pixel size in mm
V0=zeros(size(image, 1), size(image, 2), size(image, 3) ,length(s)); % allocate space for each output
for kk=1:length(s) % loop over values of s
    V0(:,:,:,kk)=filter3D(image,s(kk),ps);
end

%Compute equation 14
V0f=max(V0,[],4); % Compute eq. 14

%Plot V0f and original image
%Slice of original image
figure,pcolor(image(:, :, 10)),axis image, ,title('Original image'), xlabel('x-coordinates'), ylabel('y-coordinates')
shading interp
colormap spring
colorbar

%Plot a slice of filtered image
figure,pcolor(V0f(:, :, 10)),axis image, ,title('Filtered image, V0f'), xlabel('x-coordinates'), ylabel('y-coordinates')
shading interp
colormap spring
colorbar

%Find snr
signal_to_noise_ratio = psnr(image, imageCopy)

%Find dice
imageLog = imageCopy > 0;
[dice, p] = MaxDiceScore(imageLog, V0f)

%Plot dice
txt = sprintf('Filtered image, V0f > %.2f', p)
figure, pcolor(V0f(:, :, 10) > p), axis image, title(txt)
xlabel('x-coordinates'), ylabel('y-coordinates')
shading interp, colormap spring, colorbar







