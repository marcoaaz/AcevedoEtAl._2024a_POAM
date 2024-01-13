function [S_sub_filtered, stats_add_filtered, quiver_vectors...
    ] = targetObjectQuiverP(img_fullRes, S_sub, pixel_max_x, theta_inclination, ...
    mask_fg, scale_found, pos_ROI, filter_expression, plotOption, destFile_ori_target, saveOption)

%Reads the calculated orientation images and show a quiver plot for the
%target mineral objects

offset_fix = 1; %pos_ROI was ceil()

%upscaling
pixel_max_x_us = imresize(pixel_max_x, 1/scale_found, "nearest"); 
theta_inclination_us = imresize(theta_inclination, 1/scale_found, "nearest"); 
mask_fg_us = imresize(mask_fg, 1/scale_found, "nearest"); 

%medicine: matching sizes (for plot)
[r_num1, c_num1, ~] = size(pixel_max_x_us); %smaller
[r_num2, c_num2, ~] = size(img_fullRes);
r_num = min(r_num1, r_num2);
c_num = min(c_num1, c_num2);

%Using masks (~3 min) (from qupathPhaseMap_v6.m)
n_poly_sub = length(S_sub);
stats_add = cell(1, n_poly_sub);
% for i = 1:n_poly_sub
parfor i = 1:n_poly_sub
    sprintf('i= %d', i)

    X = S_sub(i).X;
    Y = S_sub(i).Y;
   
    stopX = find(isnan(X)) - 1;
    stopY = find(isnan(Y)) - 1;
    temp = poly2mask(X(1:stopX), Y(1:stopY), r_num, c_num);%not NaN
    
    %medicine 1: polygon vertices X and Y
    S_sub(i).X = X - (pos_ROI(2) - offset_fix); %correcting for crop
    S_sub(i).Y = Y - (pos_ROI(1) - offset_fix);

    %medicine 2: isolated pixels
    CC = bwconncomp(temp);
    numPixels = cellfun(@numel, CC.PixelIdxList);
    [~, idx] = max(numPixels);
    temp = zeros(size(temp), 'logical');
    temp(CC.PixelIdxList{idx}) = 1;        

    %object information
    mask_valid = temp & mask_fg_us(1:r_num, 1:c_num);
    px_O = pixel_max_x_us(mask_valid);
    px_I = theta_inclination_us(mask_valid);
    
    %medicine: filtering out -1 population
    area_temp = sum(temp, 'all');
    fg_temp = (px_O ~= -1) & (px_I ~= -1); %useful data
    area_pos_temp = sum(fg_temp, 'all');
    foreground = (area_pos_temp > 0.5*area_temp); %test    
    
    %Mode of the histograms
    [N_o, edges_o] = histcounts(px_O(fg_temp), 40);
    [~, idx_o] = max(N_o);
%     idx_o = find(idx_o, 1); %only care about the first peak location (debug most cases)
    mode_O = (edges_o(idx_o) + edges_o(idx_o + 1))/2;
        
    [N_i, edges_i] = histcounts(px_I(fg_temp), 40);
    [~, idx_i] = max(N_i);
%     idx_i = find(idx_i, 1); 
    mode_I = (edges_i(idx_i) + edges_i(idx_i + 1))/2;    

    stats = regionprops(temp, 'Area', 'Perimeter', 'Centroid', ...
        'Orientation', 'Circularity', 'Eccentricity', ...
        'Extent', 'Solidity', 'EulerNumber');    
    stats.mode_O = mode_O;
    stats.mode_I = mode_I; 
    stats.foreground = foreground;

%     stats_add = [stats_add; stats];
    stats_add{i} = stats;
end
stats_add = [stats_add{:}]; %as structure

%Filtering
filter_stats = eval(filter_expression);
S_sub_filtered = S_sub(filter_stats);
stats_add_filtered = stats_add(filter_stats);
pct_passed = 100*sum(filter_stats)/length(filter_stats);
sprintf('Percentage of target objects: %s', num2str(pct_passed))

%vector origin
data = reshape([stats_add_filtered.Centroid], 2, []); %match w/ object centroid
X_new = data(1, :) - (pos_ROI(2) - offset_fix); %correcting for crop
Y_new = data(2, :) - (pos_ROI(1) - offset_fix);

%Sexagesimal degrees (rotation exp.) to unitary vector projection
x_proj = 1*cos((pi/180)*[stats_add_filtered.mode_O]).*cos((pi/180)*[stats_add_filtered.mode_I]);
y_proj = 1*sin((pi/180)*[stats_add_filtered.mode_O]).*cos((pi/180)*[stats_add_filtered.mode_I]);

quiver_vectors = cat(1, X_new, Y_new, x_proj, y_proj);

%% Plot

if plotOption == 1 %preliminary to evalue shape threshold
    
    nbins = 40;

    figure %filter histograms
    
    subplot(4, 1, 1)
    histogram([stats_add.Circularity], nbins)
    title('Circularity')
    subplot(4, 1, 2)
    histogram([stats_add.Eccentricity], nbins)
    title('Eccentricity')
    subplot(4, 1, 3)
    histogram([stats_add.Extent], nbins)
    title('Extent')
    subplot(4, 1, 4)
    histogram([stats_add.Solidity], nbins)
    title('Solidity')
end

%% Save

if saveOption == 1
    %medicine (un-segmented fractures and dark pixels are vertical)
    filtered_orientation = [stats_add_filtered.mode_O];
    filtered_inclination = [stats_add_filtered.mode_I];
    
    idx_out = (filtered_inclination > 89); %criteria (default == 90)
    filtered_orientation2 = filtered_orientation(~idx_out);
    filtered_inclination2 = filtered_inclination(~idx_out);
    
    %medicine (MTEX Euler angles)
    filtered_orientation3 = filtered_orientation2';
    filtered_inclination3 = 90 - filtered_inclination2';
    
    %optional: save target grid
    n_obj_filtered = sum(~idx_out);
    sprintf('Total objects >89 degrees = %d (removed)', sum(idx_out))     
    
    %{'phi1', 'Phi', 'phi2'}
    orientation_save = [filtered_orientation3, filtered_inclination3, zeros(n_obj_filtered, 1)];
    writematrix(orientation_save, destFile_ori_target)
    
    plotMTEXstereonet(destFile_ori_target)
end

end