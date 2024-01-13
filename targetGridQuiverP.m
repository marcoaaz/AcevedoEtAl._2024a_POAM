function [quiver_vectors] = targetGridQuiverP(img_fullRes_crop, mask_bg1_crop, ...
    pixel_max_x, theta_inclination, mask_fg, scale_found, pos_ROI)
%Reads the calculated orientation images (masked) and show a quiver plot of
%the grid for the target mineral mask

%Input:
%img_fullRes_crop = cropped blended image
%pixel_max_x = c-axis orientation
%theta_inclination = angle of c-axis with horizontal plane
%scale_found = downscaling scale to restore
%mask_fg = downscaled foreground
%pos_ROI = bounding box of upscaled background
%mask_bg1_crop = upscaled and cropped background

%Output:
%img_fullRes_crop = fused image <5000 and cropped
%quiver_vectors = full resolution vectors of the target mineral
%scale_found = approx. scale between 

%% vector origin
[n_rows1, n_cols1] = size(mask_fg);
scale_restore = 1/scale_found; %upscale grid
offset_fix = 0.5; %0.5

[X, Y] = meshgrid(1:n_cols1, 1:n_rows1); %downscaled grid (= processing, not resample)
X_new = scale_restore*(X) - (pos_ROI(2) - offset_fix); %match w/ pixel centre
Y_new = scale_restore*(Y) - (pos_ROI(1) - offset_fix);

%unitary vector projection
x_proj = 1*cos((pi/180)*pixel_max_x).*cos((pi/180)*theta_inclination);
y_proj = 1*sin((pi/180)*pixel_max_x).*cos((pi/180)*theta_inclination);

%cleaning up
R_temp = img_fullRes_crop(:, :, 1);
G_temp = img_fullRes_crop(:, :, 2);
B_temp = img_fullRes_crop(:, :, 3);
R_temp(mask_bg1_crop ~= 0) = 0;
G_temp(mask_bg1_crop ~= 0) = 0;
B_temp(mask_bg1_crop ~= 0) = 0;
img_fullRes_crop(:, :, 1) = R_temp;
img_fullRes_crop(:, :, 2) = G_temp;
img_fullRes_crop(:, :, 3) = B_temp;
X_new(mask_fg == 0) = [];
Y_new(mask_fg == 0) = [];
x_proj(mask_fg == 0) = [];
y_proj(mask_fg == 0) = [];

quiver_vectors = cat(1, X_new, Y_new, x_proj, y_proj);

end