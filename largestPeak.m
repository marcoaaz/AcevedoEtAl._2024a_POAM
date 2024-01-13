function [theta_inclination, pixel_cAxisI_discrete_rs] = largestPeak(...
    img_temp2, pol_angle1, pixel_cAxis, mask_fg, intensityConstraint)

%Input:
%img_temp2 = multi-angular polarized image stack that was converted to greyscale and downscaled (0.5)
%pol_angle1 = angular value of interpolation range
%pixel_cAxis = image with discrete anglular value of fitted curve greatest maxima
%mask_fg = target grains

%intensityConstraint: 
%top row: min/max possible intensity in 3D rotation (universal stage)
%bottom row: min/max found intensity in 2D rotation (near vertical/horizontal crystal)

%Output: 
% pixel_cAxisI_discrete = image with discrete empirical value of fitted curve
%maximum 
% theta_inclination = image with continuous quartz c-axis inclination 
% values (relative to the horizontal plane of the section surface). 
% Following Goodchild (1998) MSc thesis

%%
n_layers = size(img_temp2, 3); %4 dim (past version)
n_rows1 = size(pixel_cAxis, 1);
n_cols1 = size(pixel_cAxis, 2);

%Finding: max_f (largest peak)
diff_1 = abs(pixel_cAxis - repmat(reshape(pol_angle1, 1, 1, []), [n_rows1, n_cols1]));
[~, closestIndex] = min(diff_1, [], 3);

vec1 = reshape(img_temp2, [], n_layers);
vec2 = [closestIndex(:)]; %repeat for each channel
vec3 = zeros(size(vec1, 1), 1, 'double');
for m = 1:n_layers
    temp_idx = (vec2 == m);
    vec3(temp_idx) = vec1(temp_idx, m);
end

pixel_cAxisI_discrete = reshape(vec3, n_rows1, n_cols1, 1); 
intensity_img = pixel_cAxisI_discrete(mask_fg);
% multiply intensity by 255 if using HSV-V

%Assumed limits (calibration thin section required)
I_min = intensityConstraint(1, 1); 
I_max = intensityConstraint(1, 2); 
low_T = intensityConstraint(2, 1); 
high_T = intensityConstraint(2, 2); 

undersaturated_px = (pixel_cAxisI_discrete <= low_T) & mask_fg;
oversaturated_px = (pixel_cAxisI_discrete >= high_T) & mask_fg;

mask_RGB = 255*repmat(mask_fg, [1, 1, 3]);
R_temp = mask_RGB(:, :, 1);
G_temp = mask_RGB(:, :, 2);
B_temp = mask_RGB(:, :, 3);
G_temp(oversaturated_px) = 0;
B_temp(oversaturated_px) = 0;
R_temp(undersaturated_px) = 0;
G_temp(undersaturated_px) = 0;
mask_RGB(:, :, 1) = R_temp;
mask_RGB(:, :, 2) = G_temp;
mask_RGB(:, :, 3) = B_temp;

% Crystallographic c-axis inclination 
% Goodchild (1998); Fueten and Goodchild (2001)
pixel_cAxisI_discrete_rs = rescale(pixel_cAxisI_discrete, low_T, high_T, ...
    'InputMin', low_T, 'InputMax', high_T);

% %measured from the vertical
% theta_inclination = (180/pi)*0.5*...
%     acos(1 - 2*(pixel_cAxisI_discrete_rs - I_min)./(I_max - I_min)); 

%measured from the horizontal
theta_inclination = (180/pi)*0.5*...
    acos(2*(pixel_cAxisI_discrete_rs - I_min)./(I_max - I_min) - 1); 

disp('Discrete inclination and largest peak intensity found.')

%% Informative plots

hFig = figure;
hFig.Position = [100, 100, 1400, 650];

subplot(1, 2, 1)
histogram(intensity_img)
xline([low_T, high_T], 'LineWidth', 2, 'Color', 'black')
title('Channel histogram')
subplot(1, 2, 2)
imshow(mask_RGB)
title('Oversaturated pixels')

end