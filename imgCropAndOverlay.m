function [pos_ROI, img_overlay] = imgCropAndOverlay(mask_bg1, img_input, img_cAxis_rgb2)

%Input:
%mask_bg1 = upscaled background
%img_input = full resolution optical image
%img_cAxis_rgb2 = orientation map following colour wheel (HSV)

%Output:
%pos_ROI = bounding box for cropping all maps and relocating coordinates
%img_overlay = aesthetic background to orientation maps

%%

%Cropping ROI in full resolution
ROI = regionprops(~mask_bg1, 'BoundingBox'); %top left X, Y, W, H
ROI = ROI.BoundingBox;
tl_row = ceil(ROI(2));
tl_col = ceil(ROI(1));
br_row = tl_row + ROI(4) - 1;
br_col = tl_col + ROI(3) - 1;
pos_ROI = [tl_row, tl_col, br_row, br_col];

%Map aesthetics
% img_fullRes_crop = img_fullRes(pos_ROI(1):pos_ROI(3), pos_ROI(2):pos_ROI(4), :);
img_fullRes_crop = img_input(pos_ROI(1):pos_ROI(3), pos_ROI(2):pos_ROI(4), :);
 
% Combining FoV optical images for maps
%medicine: matching sizes (for plot)
[r_num1, c_num1, ~] = size(img_cAxis_rgb2); %smaller
[r_num2, c_num2, ~] = size(img_fullRes_crop);
r_num = min(r_num1, r_num2);
c_num = min(c_num1, c_num2);

r = img_cAxis_rgb2(1:r_num, 1:c_num, 1);
g = img_cAxis_rgb2(1:r_num, 1:c_num, 2);
b = img_cAxis_rgb2(1:r_num, 1:c_num, 3);
mask_overlay = ~((r == 0) & (g == 0) & (b == 0));
R = img_fullRes_crop(1:r_num, 1:c_num, 1);
G = img_fullRes_crop(1:r_num, 1:c_num, 2);
B = img_fullRes_crop(1:r_num, 1:c_num, 3);
R(mask_overlay) = r(mask_overlay);
G(mask_overlay) = g(mask_overlay);
B(mask_overlay) = b(mask_overlay);

img_overlay = cat(3, R, G, B);

end