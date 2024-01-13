function [img_rgb2, mask_mineral_crop, pos_ROI_ds] = plotTargetAzimuthMap(...
    mask_fg, mask_bg, pixel_cAxis, theta_inclination, ...
    plotOption, destFile_targetGridStereo, saveOption)

%Inteprets the foreground and target mask and plots a 2D orientation
%(azimuth) map in a turbo colormap [0, 180]. It optionally saves all the
%available pixels in the target mask (too many pixels might be difficult to
%plot in Stereonet)

%Input:
%label_map_rs = downscaled label map (partial size 5000x5000)
%pixel_cAxis = low-resolution azimuth map of maximum peak (c-axis)
%mask_bg = low-resolution background (-1 value, 70% of area)

%Output:
%label_map_rs_crop = downscaled and cropped label map
%img_cAxis_rgb = cropped and coloured Azimuth map

%Cropping ROI (downscaled)
ROI = regionprops(~mask_bg, 'BoundingBox'); %top left X, Y, W, H
ROI = ROI.BoundingBox;
tl_row = ceil(ROI(2));
tl_col = ceil(ROI(1));
br_row = tl_row + ROI(4) - 1;
br_col = tl_col + ROI(3) - 1;
pos_ROI_ds = [tl_row, tl_col, br_row, br_col]; %downscaled

%Cropping: Label map and orientation: Cropping
mask_mineral_crop = mask_fg(pos_ROI_ds(1):pos_ROI_ds(3), pos_ROI_ds(2):pos_ROI_ds(4));
pixel_cAxis_crop = pixel_cAxis(pos_ROI_ds(1):pos_ROI_ds(3), pos_ROI_ds(2):pos_ROI_ds(4));
theta_inclination_crop = theta_inclination(pos_ROI_ds(1):pos_ROI_ds(3), pos_ROI_ds(2):pos_ROI_ds(4));

%Colour direction wheel (Wilson et al., 2007)
top_blackness = 0.5; %wheel center
img_hue = rescale(pixel_cAxis_crop, 0, 1, 'InputMin', 0, 'InputMax', 180);
img_saturation = 1*ones([ROI(4), ROI(3)]);
img_value = 1- rescale(theta_inclination_crop, 0, top_blackness, 'InputMin', 0, 'InputMax', 90);
img_hsv = cat(3, img_hue, img_saturation, img_value);
img_rgb = hsv2rgb(img_hsv); %[0-1]

%cleaning up
R_temp = img_rgb(:, :, 1);
G_temp = img_rgb(:, :, 2);
B_temp = img_rgb(:, :, 3);
R_temp(~mask_mineral_crop) = 0; %ignored phases
G_temp(~mask_mineral_crop) = 0;
B_temp(~mask_mineral_crop) = 0;
img_rgb(:, :, 1) = R_temp;
img_rgb(:, :, 2) = G_temp;
img_rgb(:, :, 3) = B_temp;

img_rgb2 = 255*img_rgb;

%% Plot

cmap2 = hsv(255); %hsv for 360, 'turbo' better than 'jet' for 180

if plotOption == 1
    
    %informative plots
    hFig = figure;
    hFig.Position = [100, 100, 1200, 600];

    subplot(1, 2, 1)
    histogram(theta_inclination_crop(mask_mineral_crop))
    title('Discrete inclination (degrees)')
    subplot(1, 2, 2)
    histogram(img_value(mask_mineral_crop))
    title('HSV - Value')

    %wheel coloured map
    min_bar = 0; 
    max_bar = 180;
    interval_val = 20; 
    thick_vals_orig = min_bar:interval_val:max_bar;
    thick_labels = strsplit(num2str(thick_vals_orig));
    thick_vals = rescale(thick_vals_orig, 0, 1); %scaled to triplet
    
    hFig2 = figure;    
    ax = gca;

    imshow(img_rgb, cmap2)
    cb = colorbar(gca, 'Location', 'eastoutside', 'Ticks', thick_vals, ...
        'TickLabels', thick_labels);    
    title({'Pixel orientation by azimuth (hue)', ...
        'and inclination (value)'})

    hFig2.WindowState = 'maximized'; 
end

%% Save

if saveOption == 1   
    
    %medicine (un-segmented fractures and dark pixels are vertical)
    mask_vertical = (theta_inclination_crop == 90);
    mask_mineral_crop2 = mask_mineral_crop & ~mask_vertical;
    %optional: save target grid
    n_px_filtered = sum(mask_mineral_crop2, 'all');
    
    filtered_orientation3 = pixel_cAxis_crop(mask_mineral_crop2);
    filtered_inclination3 = 90 - theta_inclination_crop(mask_mineral_crop2);
       
    %{'phi1', 'Phi', 'phi2'}
    orientation_save = [filtered_orientation3, filtered_inclination3, zeros(n_px_filtered, 1)];
    writematrix(orientation_save, destFile_targetGridStereo)
    
    plotMTEXstereonet(destFile_targetGridStereo)    
end

end