load('projs.mat'); % 144x513 sinogram array called p

                    % origin of each projection in middle (257th element)

                    % each row represents angles from 0 to 178.75 degrees in 1.25 degree increments

N = 512;


% ART Method - Iterative Reconstruction Technique

guess_img = zeros(N+1); % Guess at object matrix (all zeros)

guess_Projection = zeros(1,N+1); % define projection array


rotation_angle = 1.25; % angle step(phi) 1.25 degrees


tic %start timer for ART reconstruction process

for k=1:10

    for i=1:144

        guess_Projection(1,:) = sum(guess_img); % calculate projection of guess at angle phi

        error = (p(i,:) - guess_Projection(1,:))/(N+1); % compare calculated projection of guess with actual projection data from object at same angle

        y = repmat(error,N+1,1); % replicates row to form matrix

        guess_img = y + guess_img; % adjust guess to match actual data projection for that angle

        guess_img = imrotate(guess_img,-rotation_angle,'bicubic','crop'); % rotate "guess" in clockwise direction

    end

      guess_img = imrotate(guess_img,180,'bicubic','crop'); % rotate "guess back to original 0 degrees and repeat process"

end

toc %end timer


tic %start timer for Direct Method reconstruction process

% Direct Method (using Central Slice Theorem)

p = p(:,1:N); % truncate last column of zeros to perform fft

p = fftshift(p,2); % flip along columns so that origin of each row is first element

CentralSlices = fft(p,[],2); % 1D FFT for each row to get Central Slices

CentralSlices = fftshift(CentralSlices,2); % relocate origin to center of object; flip along columns


% Information for griddata() to interpolate to cartesian grid %

Angles = (0:1.25:178.75).'; % make column vector of angles in degrees with step 1.25

cosAngles = cosd(Angles); % cosine of angles

sinAngles = sind(Angles); % sine of angles

wxp = (-N/2:(N/2-1)); % each column in Central Slice array indicates wxp

wx = cosAngles*wxp; % omega_x locations corresponding to "known" data locations

                      % omega_x = omega_x'cos(phi)

wy = sinAngles*wxp; % omega_y locations corresponding to "known" data locations

                      % omega_y = omega_y'sin(phi)

wx_target = zeros(N); % "target" omega_x should be (2^p)x(2^p) which is size of reconstructed image

for i = 1:N

    wx_target(i,:) = wxp;

end

wy_target = wx_target.'; % "target" omega_y is transpose of omega_x


interpolation = griddata(wx,wy,CentralSlices,wx_target,wy_target,'linear');

% origin element is in location (N/2+1,N/2+1) %

interpolation = fftshift(interpolation); % put origin in upper left corner for ifft


% search for the NaN elements and replace them with zeros %

nanloc=find(isnan(interpolation)); %find NaN locations

interpolation(nanloc)=zeros(size(nanloc)); %now zeros where NaN's were


reconstructed = ifft2(interpolation); % ifft to get back reconstructed image

reconstructed = fftshift(reconstructed); % put origin back to middle of object


reconstructed = flipud(reconstructed); % flip image due to Matlab considering +y-axis vertically down

                                     % In original image, +y direction is the top of image

reconstructed(513,513) = 0; % pad zeros to match dimension (N+1)x(N+1)

abs_reconst = abs(reconstructed); % remove imaginary part for error calculations

toc %end timer


Original_img = imread('orig.tiff'); % read Original Image data


figure(1) % Compare Original Image with ART and Direct Reconstructions

subplot(1,3,1); imshow(Original_img,[]); title('Original Image');                                    

subplot(1,3,2); imshow(guess_img,[]); title('ART Reconstruction');                                    

subplot(1,3,3); imshow((abs_reconst),[]); title('Direct Method Reconstruction');


%Reconstruction Comparison through Error Analysis%

abs_Original_img = Original_img - min(min(Original_img)); %remove negative valued data from Original Image

scaled_Original_img = abs_Original_img/max(max(Original_img)); %scale Original image by normalizing with max value


abs_guess_img = guess_img - min(min(guess_img)); %remove negative valued data from ART reconstructed image

scale_guess_img = abs_guess_img/max(max(guess_img)); %scale ART reconstructed image by normalizing with max value

 

scale_reconstructed = abs_reconst/max(max(abs_reconst)); %scale Direct reconstructed image by normalizing with max value


error1 = double(scaled_Original_img) - scale_guess_img; % compare Original vs ART reconstruction

error2 = double(scaled_Original_img) - scale_reconstructed; % compare Orginal vs Direct reconstruction

error3 = scale_guess_img - scale_reconstructed; % compare ART vs Direct reconstructions


% mesh plot error analysis

figure(2)

mesh(abs(error1));title("Original vs ART");xlabel("rows");ylabel("columns");zlabel("ERROR")

figure(3)

mesh(abs(error2));title("Orginal vs Direct");xlabel("rows");ylabel("columns");zlabel("ERROR")

figure(4)

mesh(abs(error3));title("ART vs Direct");xlabel("rows");ylabel("columns");zlabel("ERROR")