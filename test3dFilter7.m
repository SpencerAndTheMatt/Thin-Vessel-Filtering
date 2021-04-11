%test3dFilter7
%Tests overlapping vessels in noisy conditions
%Computes signal to noise ratio and potential dice score

%Form image and imageLog
image = zeros(20, 20, 20);
imageLog = zeros(20, 20, 20);

%Add noise to image
noise = 0.6;
image = image + noise * rand(size(image));

%Add overlapping vessels to image and imageLog
image(9:11, 9:11, :) = 1; %all through z
imageLog(9:11, 9:11, :) = 1;

image(5:8, :, 5:8) = 1;
imageLog(5:8, :, 5:8) = 1;

image(:, 14:15, 5:8) = 1;
imageLog(:, 14:15, 5:8) = 1;

%Activate filter
s=1:0.1:6; % values of s to test
ps=1; % pixel size in mm
V0=zeros(size(image, 1), size(image, 2), size(image, 3) ,length(s)); % allocate space for each output
for kk=1:length(s) % loop over values of s
    V0(:,:,:,kk)=filter3D(image,s(kk),ps);
end

%Compute equation 14
V0f=max(V0,[],4); % Compute eq. 14

%Find signal to noise ratio
signal_to_noise_ratio = psnr(image, imageLog)

%Compute max dice score
[maxDiceScore, P] = MaxDiceScore(image > noise, V0f)

%Plot results
figure,pcolor(image(:, :, 7)),axis image, ,title('Original image')
colormap spring
colorbar

%Plot a slice of filtered image
figure,pcolor(V0f(:, :, 7)),axis image, ,title('Filtered image, V0f')
colormap spring
colorbar

%Plot v0f > p
figure, pcolor(V0f(:, :, 7) > P), axis image, title('Filtered image, V0f > P')
colormap spring
colorbar











