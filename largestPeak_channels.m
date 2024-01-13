function [pixel_maxPeakI_discrete] = largestPeak_channels(...
    img_temp2, pol_angle1, pixel_maxPeak)

%[theta_inclination, pixel_cAxisI_discrete_rs] = 

%Input:
%img_temp2 = multi-angular polarized image stack (e.g., downscaled 0.5)
%pol_angle1 = angular value of interpolation range
%pixel_maxPeak = image with continuous anglular value of fitted curve greatest maxima

%%
n_rows1 = size(pixel_maxPeak, 1);
n_cols1 = size(pixel_maxPeak, 2);
n_channels = size(img_temp2, 3); 
n_layers = size(img_temp2, 4); 

%Finding: max_f (largest peak)
diff_1 = abs(pixel_maxPeak - repmat(reshape(pol_angle1, 1, 1, []), [n_rows1, n_cols1]));
[~, closestIndex] = min(diff_1, [], 3);

vec1 = reshape(img_temp2, [], n_channels, n_layers);
vec2 = [closestIndex(:)]; %repeat for each channel
vec3 = zeros(size(vec1, 1), n_channels, 'double');
for m = 1:n_layers
    temp_idx = (vec2 == m);
    vec3(temp_idx, :) = vec1(temp_idx, :, m);
end

pixel_maxPeakI_discrete = reshape(vec3, n_rows1, n_cols1, n_channels);

end