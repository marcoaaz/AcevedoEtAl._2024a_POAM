%Root directory
%Writen by Marco Acevedo Zamora, QUT - 8-Aug-2022
%Updated: 
% M.A., 14 april 2023
% M.A., 8 may 2023

% Reads on the registered image stack (from TrakEM2), it models stage
% rotation TL/RL-PPL/XPL/XPL-lambda waves, produces 2D c-axis orientation,
% and does ray tracing, image blending, and plot object quiver plots.

%dependency for stereonets
mtex_folder = 'E:\Alienware_March 22\current work\00-new code May_22\mtex-5.9.0';
cd(mtex_folder);
startup_mtex

%% User input:
close all
clear 
clc

scriptDir = 'E:\Alienware_March 22\current work\00-new code May_22';
scriptDir2 = 'E:\Alienware_March 22\scripts_Marco\updated MatLab scripts\';
scriptDir3 = fullfile(scriptDir2, 'ROI');

addpath(scriptDir)
addpath(fullfile(scriptDir, 'rayTracing/'))     
addpath(fullfile(scriptDir2, 'ROI'))
addpath(fullfile(scriptDir2, 'plots_miscellaneous'))

%Paper datasets
%E:\paper 2_datasets\nikon LV1000ND\mylonite\final stacks
%E:\paper 2_datasets\nikon LV1000ND\mylonite\zoom-in registeredStacks
%E:\paper 2_datasets\nikon LV1000ND\cpx\stacked\saved_aligned\stack_modality
%E:\paper 2_datasets\nikon LV1000ND\granulite 6KB-67\composite
%E:\paper 2_datasets\nikon LV1000ND\TS-J2KB\stacks\composite

%Dataset definition
workingDir = 'E:\paper 2_datasets\nikon LV1000ND\mylonite\zoom-in registeredStacks';
fileName = 'all_modalities.tif'; %TrakEM2 image stack
blendedImage_suffix = 'test1'; %SuperSIAT input
fileName_spf = 'sample2.spf'; %SuperSIAT training points
fileName_shp = 'sample2_segmentation.shp'; %objects
destDir_suffix = 'review3';%trial saving directory
plotSetup = struct;
plotSetup.lineWidth = 0.6; %vector thickness
plotSetup.autoSFactor = .4; %quiver vectors size
plotSetup.edgeAlpha = 1; %edge transparency

%plotSetup.lineWidth = 0.6; %mylonite and garnet
%plotSetup.autoSFactor = 0.4; %mylonite and garnet
%plotSetup.edgeAlpha = 1; %coarse grained mineral

%Define: Multi-pol imaging spectra
sel_modality = {'RL PPL', 'RL XPL', 'TL PPL', 'TL XPL', 'TL XPL-lambda'}; %acquired series
rlPPL_range = 1:36; 
rlXPL_range = 37:72;
tlPPL_range = 73:108;
tlXPL_range = 109:144;
tlXPLlambda_range = 145:180;

% %Example (paper 2)
% rlPPL_range = 1:36; %2:7
% rlXPL_range = 37:72; %8:19
% tlPPL_range = 73:108; %20:25
% tlXPL_range = 109:144; %26:n_list
% tlXPLlambda_range = 145:180;

%Define image blending
type = {'max'}; %used ray tracing method

%Mylonite
modality = {'RL PPL', 'TL PPL', 'TL XPL', ...
    'rescaled8_pca_25april_pctOut0.5_denoised3x3_registered'}; %'RL XPL', 'TL PPL', 'TL XPL-lambda'
weights = [15, 20, 20, 45];
% %Every sample (not the harzburgites)
% modality = {'RL PPL', 'TL PPL', 'TL XPL'}; %'RL XPL', 'TL PPL', 'TL XPL-lambda'
% weights = [50, 25, 25];

%Define Algorithm 1
sel = 5; %select idx of 'TL-XPL-lambda' modality 
pol_angle1 = [ 0:10:350 ]'; % experiment range (transposed to match sequezeed px)
scale = 0.5; % Decreases computational cost
cluster = 10; %#cores >8 crashed 32GB (dont forget to turn off debugging plots!)
period_modality = [180, 90, 180, 90, 180]; % optical_period; XPL=90, PPL=180

%Saving file name of Algorithms 1 & 2 loop
descripTxt = {'30Dec23'}; %name of modulation image file
% descripTxt = strcat({'R', 'G', 'B'}, '_');

%% Script: Importing single image stacks (patch)
cd(workingDir)

%Note: <4GB tif, or requires Bioformat Exporter plugin
imageName = fullfile(workingDir, fileName);
struct1 = imfinfo(imageName); 
n_rows = struct1(1).Height;
n_cols = struct1(1).Width;
n_channels = struct1(1).BitDepth/8;
n_images = length(struct1);
% list = setdiff(1:n_images, 2); %exclude images
list = setdiff(1:n_images, []); %ALL
n_list = length(list);

%Creating folders
destFolder = fullfile(workingDir, strrep(fileName, '.tif', '_falseColor'));
destFolder_rt = fullfile(workingDir, strrep(fileName, '.tif', '_rayTracing'));
mkdir(destFolder)
mkdir(destFolder_rt)

n_modalities = length(sel_modality);

sel_range = {rlPPL_range, rlXPL_range, tlPPL_range, tlXPL_range, tlXPLlambda_range};
sel_numberLayers = zeros(1, n_modalities);
for i = 1:n_modalities
    sel_numberLayers(i) = length(sel_range{i});
end
%Optional (uncomment for single): Data reduction parameters
n_layers = sel_numberLayers(sel); %

%Pre-check orientation map output size (assuming all images are equal)
img_temp_check = zeros(n_rows, n_cols, 'uint8');
img_temp2_check = imresize(img_temp_check, scale, "bilinear"); %not real
[n_rows1, n_cols1] = size(img_temp2_check);
n_pixels1 = n_rows1*n_cols1;

%Section 1: execute Ray tracing

%Informative structure
info_struct.Height = n_rows; 
info_struct.Width = n_cols;
info_struct.Channels = n_channels;
info_struct.n_tile_layers = n_modalities;
info_struct.sel_range = sel_range;
info_struct.sel_modality = sel_modality;

% {'mean', 'max', 'min', 'range', 'sum', 'std', 'median', 'maxHSV', 'minHSV', 'rangeHSV'}
stats_list = {'max', 'min'};
n_options = length(stats_list);

time_elapsed = zeros(1, n_options);
for k = 1:n_options     
    %parallel computing (default: commented)
%     [time_elapsed] = stats_zProject(imageName, info_struct, stats_list{k}, destFolder_rt);
end

%Section 2: execute Image Blending, fuse ray tracing images (OBIAS input)
n_mode = length(modality);
strMode = [];
for j = 1:n_mode
    strMode = strcat(strMode, '_', modality{j});
end
destFile = fullfile(destFolder_rt, strcat('imgSum', strMode, '_', blendedImage_suffix, '.tif')); %test1_W5000

%Save (run once)
saveOption = 0;
[img_sum_w_rs] = stack_imageFusion(type, modality, weights, destFile, saveOption);

%Defining directories
%Blended image
[img_wkDir, temp] = fileparts(destFile);
img_file = fullfile(img_wkDir, strcat(temp, '.tif')); 
%SuperSIAT files
workingDir_obias = strrep(img_file, '.tif', ''); %SuperSIAT convention
fileName_TIF = strrep(fileName_shp, '.shp', '.tif');
temp_name = fullfile(workingDir_obias, fileName_spf);

%% Optional: Spectra live mode

sel_modes = [1:n_modalities];

img_A = img_sum_w_rs;
[A_cols, A_rows, ~] = size(img_A);

%ROI: 3000x3000x3x180 uint8 is 5GB and slow
roi_width = 2000;
roi_height = roi_width;
% roi_tl_row = 2500; %TS-J2KB 5X garnet_location
% roi_tl_col = 2500;
% roi_tl_row = 1500; %mylonite 5x
% roi_tl_col = 1500;
roi_tl_row = 500; %default
roi_tl_col = 500;
roi_br_row = roi_tl_row + roi_height - 1;
roi_br_col = roi_tl_col + roi_width - 1;

%Importing full stack (3 min)
img_full_temp = zeros(roi_height, roi_width, n_channels, n_list, 'uint8');

k = 0;
for i = 1:n_list
    k = k + 1;    
    temp_img = imread(imageName, i); 
    %ROI is linearised inside ROI_modulation_data.m     

    %RGB - all
    img_full_temp(:, :, :, k) = temp_img(roi_tl_row:roi_br_row, roi_tl_col:roi_br_col, :);
end  
%%
close all 
%Plot ROI
figure

hImage = imshow(img_A(roi_tl_row:roi_br_row, roi_tl_col:roi_br_col, :));
ax = gca; %alternative: fig = gcf; ax = fig.CurrentAxes;

%'Center', [n_cols/2 n_rows/2],...
hCircle = images.roi.Circle(...
    'Center', [roi_height/2, roi_width/2],...
    'Radius', 12, 'Parent', ax, 'Color', 'r');

%pol_angle1: has to be a row
color_space = 1; %1= RGB, 2=CIElab, 3=HSV
sel_modes = [1:5];
addlistener(hCircle, 'MovingROI',...
    @(varargin)ROI_modulation_data(hCircle, img_full_temp, ...
    color_space, sel_range, sel_modes, period_modality, pol_angle1')); %'MovingROI'

%sel_modality = {'RL PPL', 'RL XPL', 'TL PPL', 'TL XPL', 'TL XPL-lambda'};
% sel_modes = [3, 4]; %for paper (Fig aspect was also fixed to host 2 panels)
sel_modes = [5];
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROI_modulation_graph(hCircle, sel_modes, sel_modality)); 

selected_modality = 5;
selected_channel = 4; %R, G, B, Greyscale
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROI_modulationAlgorithm_graph(hCircle, ...
    selected_modality, selected_channel)); 

%% Section 3: Subsetting image and Fitting Sin wave

cluster= 8; %default = 10

%For 1 channel
n_RGB = 1;

% %Alternative: For 3 channels
% sel_channel = [1, 2, 3];
% n_RGB = length(sel_channel);

%channel
for uu = 1:n_RGB %uu = 1:n_RGB
    
    %modalities % for sel = 1:5
    for sel = 5
        optical_period = period_modality(sel);
    
        % Preparing data for PLM modeling
        varRange = sel_range{sel};
        n_layers = length(varRange);   
        
        %select channels
        img_temp = zeros(n_rows, n_cols, n_layers, "double"); %HSV or channel
    
        k = 0;
        for i = varRange
            k = k + 1;    
            temp_img = imread(imageName, i);                       
            
            %Optional: linearize
            temp_img = rgb2lin(temp_img, 'ColorSpace', 'srgb');

            %Greyscale
            prior_temp = double(temp_img);
            img_temp(:, :, k) = mean(prior_temp, 3);

% %             RGB - each
%             img_temp(:, :, k) = double(temp_img(:, :, sel_channel(uu)));
        end            
        img_temp2 = imresize(img_temp, scale, "bilinear"); %downscaling        

        tic;

        %Algorithm 1 and 2 (run in parallel)
        [pixel_calc] = imageFourierS_optim(img_temp2, pol_angle1, ...
            optical_period, cluster);         
        
        %info:
%         pixel_calc = cat(3, pixel_avg, pixel_phase, pixel_range1, ...
%             pixel_range2, pixel_maxPeak, pixel_maxPeakI, pixel_SSE);
        % Save model 
        imgSpectraModel = struct;            
        imgSpectraModel.pixel_avg = pixel_calc(:, :, 1);
        imgSpectraModel.pixel_phase = pixel_calc(:, :, 2);
        imgSpectraModel.pixel_range1 = pixel_calc(:, :, 3);
        imgSpectraModel.pixel_range2 = pixel_calc(:, :, 4);
        imgSpectraModel.pixel_maxPeak = pixel_calc(:, :, 5); %stage rotation for max peak [0-180> 
        imgSpectraModel.pixel_maxPeakI = pixel_calc(:, :, 6);
        imgSpectraModel.pixel_SSE = pixel_calc(:, :, 7); %sum of squared errors (usefulness for prediction)
        
        file_output = strcat(descripTxt{uu}, '_', ...
            'modulationImage_', sel_modality{sel}, '_', num2str(scale), '.mat'); %'wave_TL-PPL_0.15.mat'
        destFile = fullfile(destFolder, file_output);
        save(destFile, "imgSpectraModel", '-mat', '-v7.3')
        
        time = toc;
        t.FittedImage(uu) = time/3600; %hours
        t.Rate(uu) = n_pixels1/time; %px/sec           
    
    end %remove when only Section 2 required
end

%% Optional: Modulation technique: save RGB image (requirement: Section 2) 

% descripTxt = {''}; %see user input
descripTxt = {'30Dec23'};
uu = 1; %channel
sel = 5;

file_output = strcat(descripTxt{uu}, ...
    'wave2_', sel_modality{sel}, '_', num2str(scale), '.mat');
imageFile = strrep(file_output, '.mat', '.tif');
m = load(fullfile(destFolder, file_output));
    
%Build image. Following Axer et al., 2011
R_double = m.imgSpectraModel.pixel_avg;
G_double = m.imgSpectraModel.pixel_phase;
%debug
% G_double(G_double > optical_period) = G_double(G_double > optical_period) - optical_period;
temp_B_double = m.imgSpectraModel.pixel_range1;
B_double = temp_B_double./R_double;

%zeroed mask    
mask_bg = (R_double == 0) & (G_double == 0) & (B_double == 0); 
R_double(mask_bg) = 0;
G_double(mask_bg) = 0;
B_double(mask_bg) = 0;

R = uint8(rescale(R_double, 0, 255)); %transmittance
G = uint8(rescale(G_double, 0, 255)); %direction (inclination proxy)
B = uint8(rescale(B_double, 0, 255)); %retardation
rgb_summary = cat(3, R, G, B);

imshow(rgb_summary)
imwrite(rgb_summary, fullfile(destFolder, imageFile), 'compression', 'none')
    
%% Section 4: Importing model, annotations, and grain objects 

%Import fused image
img_fullRes = imread(img_file); %full (or partial 5000x5000)

%Importing TL-XPL-lambda stack (3 min)
sel = 5;
varRange = sel_range{sel};

n_layers_temp = sel_numberLayers(sel);
img_temp = zeros(n_rows, n_cols, n_layers_temp, 'double');
k = 0;
for i = varRange
    k = k + 1;    
    temp_img = imread(imageName, i);     

    %Greyscale
    prior_temp = double(temp_img);
    img_temp(:, :, k) = mean(prior_temp, 3); 
    
%     %Option: If wanting to process a specific RGB
%     img_temp(:, :, k) = double(temp_img(:, :, sel_channel(uu)))/255;
end    
img_temp2 = imresize(img_temp, scale, "bilinear"); %downscaling        
[n_rows_ds, n_cols_ds, ~] = size(img_temp2);

%% Obtaining SuperSIAT v2.2 map information (follows version convention)

[annotationNames, cmap] = SuperSIATinterpreter(temp_name); %function

%Object structure 
S = shaperead(fullfile(workingDir_obias, fileName_shp));
range_ID = unique([S.ClassID]);
% S2 = struct2table(S);

%Label map (32-bit), objects and classes
label_maps_SS = imread(fullfile(workingDir_obias, fileName_TIF));
label_map = double(label_maps_SS(:, :, 2)) + 1; %to avoid indexing issues 

% %Alternative: Label map (8 min for 54K objects), partial resolution 5000x5000
% plotOption = 0;
% [label_map] = maskFromPoly(S, img_file, cmap, plotOption); %function 

[n_rows_partial, n_cols_partial] = size(label_map);

%scale factor: orientation maps (downscaled)/segmentation (full resolution or partial 5000x5000) 
scale_found = n_rows1/n_rows_partial; 

annotationNames' %check

%Creating Saving Folder
saveFolder = strcat('obiasExport_', num2str(scale), '_', destDir_suffix);
saveDir = fullfile(destFolder_rt, saveFolder);
mkdir(saveDir)

%% Section 5: Target C-axis estimation (discrete estimate)

sel = 5;%modality: TL-XPL-lambda
mineralSel = 2; %selecting target mineral
intensityConstraint = [135; 135];

%quartz [138; 138];. %384D22_5x
%biotite [75; 75];.
%feldspar [100; 100];.
%albite [140; 140];. 

%quartz [135; 135];. 384D22_10x.
%albite [100; 100];.

%biotite [85; 85];. 6KB-67
%plagioclase [100; 100];.
%garnet [3; 3];. 
%amphibole [80; 80];.

%garnet [12; 12];. TS-J2KB
%quartz [170; 170];. 
%biotite [160; 160];.

mineralText = annotationNames(mineralSel);

descripTxt1 = '';
file_output = strcat(descripTxt1, 'wave2_', sel_modality{sel}, '_', num2str(scale), '.mat');
fileName = fullfile(destFolder, file_output); %falseColor dir

m = load(fileName);
R_double = m.imgSpectraModel.pixel_avg;
mask_bg = (R_double == 0); %%Zeroed mask

%Ranges
pixel_range1 = m.imgSpectraModel.pixel_range1;
pixel_range2 = m.imgSpectraModel.pixel_range2;

%Depends on microscope setup (input: imgSpectraModel.pixel_maxPeak)
pixel_maxPeak = m.imgSpectraModel.pixel_maxPeak; %stage rotation value when max peak
[pixel_cAxis] = stageToAxisReorientation(pixel_maxPeak); 

%%


%Foreground mask: downscaled to match maps (assumes FoV are equal)
label_map_ds = imresize(label_map, [n_rows1, n_cols1], "nearest"); %NaN finds automatically
mask_fg = (label_map_ds == mineralSel) & (~mask_bg); %downscaled foreground (mineral inside)

[pixel_maxPeakI_discrete] = largestPeak_channels(...
    img_temp2, pol_angle1, pixel_maxPeak);

[theta_inclination, pixel_cAxisI_discrete_rs] = estimateInclination_range(...
    pixel_range1, mask_fg, intensityConstraint);

%Save csv for Stereonet (only for downscaled output) 
saveCSV = 1; %yes/no
destFile_fullStereo = fullfile(saveDir, strrep(file_output, '.mat', '_fullNet.csv'));
if saveCSV == 1    
    text_save = [pixel_cAxis(:), theta_inclination(:)];
    text_save(text_save(:, 1) == -1, :) = [];
    writematrix(text_save, destFile_fullStereo)   
end

%Optional: noise reduction
r = 1;
kernel1 = [2*r+1, 2*r+1];
pixel_cAxis2 = medfilt2(pixel_cAxis, kernel1);
theta_inclination2 = medfilt2(theta_inclination, kernel1);

%% Section 6: Target mineral plots
%plot and return downscalled and cropped label map
plotOption = 1;
saveOption = 1;
newTxt = strcat('_pxNet_', mineralText, '.csv');
destFile_targetGridStereo = fullfile(saveDir, strrep(file_output, '.mat', newTxt));

[img_cAxis_rgb, mask_mineral_crop, pos_ROI_ds] = plotTargetAzimuthMap(...
    mask_fg, mask_bg, pixel_cAxis2, theta_inclination2, plotOption, ...
    destFile_targetGridStereo, saveOption);

%for quiver plots
img_cAxis_rgb2 = uint8(imresize(255*img_cAxis_rgb, 1/scale_found, 'nearest'));

% %% Build grid quiver plot

%greyscale image: inform about texture (w/ translucent OBIAS map)
img_grey = cat(3, rgb2gray(img_fullRes), ...
    rgb2gray(img_fullRes), ...
    rgb2gray(img_fullRes)); %required to preserve cmap

%Filter; default (all) = ~mask_bg
mask_bg1 = imresize(mask_bg, [n_rows_partial, n_cols_partial]); %upscaled background (full resolution)

%Cropping ROI in full resolution
ROI = regionprops(~mask_bg1, 'BoundingBox'); %top left X, Y, W, H
ROI = ROI.BoundingBox;
tl_row = ceil(ROI(2));
tl_col = ceil(ROI(1));
br_row = tl_row + ROI(4) - 1;
br_col = tl_col + ROI(3) - 1;
pos_ROI = [tl_row, tl_col, br_row, br_col];

img_grey_crop = img_grey(pos_ROI(1):pos_ROI(3), pos_ROI(2):pos_ROI(4), :);
img_fullRes_crop = img_fullRes(pos_ROI(1):pos_ROI(3), pos_ROI(2):pos_ROI(4), :);
mask_bg1_crop = mask_bg1(pos_ROI(1):pos_ROI(3), pos_ROI(2):pos_ROI(4));

%Option: coloured or grey version for plots
img_fullRes_crop = img_grey_crop;

%calculate vectors
[quiver_vectors_grid] = targetGridQuiverP(img_fullRes_crop, mask_bg1_crop, ...
    pixel_cAxis2, theta_inclination2, mask_fg, scale_found, pos_ROI); %function

%image fusion
img_blend = imfuse(img_fullRes_crop, img_cAxis_rgb2, 'blend'); %option 1

%medicine: matching sizes (for plot)
[r_num1, c_num1, ~] = size(img_cAxis_rgb2); %smaller
[r_num2, c_num2, ~] = size(img_fullRes_crop);
r_num = min(r_num1, r_num2);
c_num = min(c_num1, c_num2);

R = img_fullRes_crop(1:r_num, 1:c_num, 1);
G = img_fullRes_crop(1:r_num, 1:c_num, 2);
B = img_fullRes_crop(1:r_num, 1:c_num, 3);
r = img_cAxis_rgb2(1:r_num, 1:c_num, 1); %option 2
g = img_cAxis_rgb2(1:r_num, 1:c_num, 2);
b = img_cAxis_rgb2(1:r_num, 1:c_num, 3);
mask_overlay = ~((r == 0) & (g == 0) & (b == 0));
R(mask_overlay) = r(mask_overlay);
G(mask_overlay) = g(mask_overlay);
B(mask_overlay) = b(mask_overlay);
img_overlay = cat(3, R, G, B);

img_fused = img_overlay; %alternative: img_blend

%Calculate object quiver plot
S_sub = S([S.ClassID] == mineralSel - 1); %indexing issue

%Run object calculations (2 min) Requires setting a manual filter for each
plotOption = 0;
saveOption = 1;
newTxt = strcat('_objNet_', mineralText, '.csv'); %Saving filtered data
destFile_ori_target = fullfile(saveDir, strrep(file_output, '.mat', newTxt));

[S_sub_filtered, stats_add_filtered, quiver_vectors_obj...
    ] = targetObjectQuiverP(img_fullRes, S_sub, pixel_cAxis2, theta_inclination2, ...
    mask_fg, scale_found, pos_ROI, plotOption, destFile_ori_target, saveOption);

%Note: Each mineral has to be informed by the plots (run only first time).
%Change the object filtering criteria inside targetObjectQuiverP following
%the plots. For instance (quartz grain contacts need):
%filter_stats = ([stats_add.foreground]' == 1 & [stats_add.Solidity]' > .35 & [stats_add.Circularity]' > .15);

%Save unfiltered object orientation vectors (used in quiver plots)
destFile_objectQuiver = strrep(destFile_ori_target, ...
    '.csv', '_centroidAndVector2.csv'); %for optional validation with EBSD
writematrix(quiver_vectors_obj, destFile_objectQuiver)

%% Plot: full OBIAS map & Gridded object quiver plot

% %trial other settings
% plotSetup.lineWidth = 0.3; %0.001 vector thickness
% plotSetup.autoSFactor = .4; %0.2 quiver vectors size

plotSetup.cmap = cmap;
plotSetup.range_ID = range_ID; 

plotSetup.transparentVal = 1; %SuperSIAT class map colour strength (default=0.3)
% plotOBIASmap(img_fullRes, S, plotSetup, annotationNames, saveDir);

%img_fused or img_gray_crop
plotSetup.transparentVal = 0;
plotTargetObjectQuiver(img_fused, S_sub_filtered, quiver_vectors_obj,  ...
    plotSetup, mineralText, saveDir)

annotationNames' %check

%% Optional (quartz only): Add-on to validate POAM with EBSD

matched_EBSD = [407, 230, 74, 210, 234, 446, 268, 369];
matched_coordinates = [
    610.623, 984.149;
    974.507, 1170.25;
    628.073, 1307.14;
    803.882, 1300.01;
    653.647, 1214.08;
    976.897, 896.574;
    501.776, 1153.44;
    833.484, 1057.95;    
];

% matched_EBSD = [402, 407, 369, 446, 230, 210, 234, 325];
% matched_coordinates = [
%     2139.14	1687.27; 
%     2227.7	1729.91;
%     2331.38	1761.37;
%     2404.69	1676.02;
%     2399.89	1816.85;
%     2351.04	1852.98;
%     2243.14	1842.83;
%     2168.26	1752.09];

%Finding ebsd
ebsdPolar = readtable('ebsdPolar.csv');
ebsdID = ebsdPolar.ID;

n_grains = length(matched_EBSD);
ebsdPolar2 = [];
for k = 1:n_grains
    idx = (ebsdID == matched_EBSD(k));
    row_temp = ebsdPolar(idx, :);
    ebsdPolar2 = [ebsdPolar2; row_temp];
end
ebsdID2 = ebsdPolar2.ID;

%Make EBSD comparable to Optical modelling
theta_temp = ebsdPolar2.theta;
rho_temp = ebsdPolar2.rho;
lowerHemi = (theta_temp > pi/2);
theta_temp = abs(pi/2 - theta_temp);
rho_temp(lowerHemi) = rho_temp(lowerHemi) + pi;
rho_temp(rho_temp < 0) = rho_temp(rho_temp < 0) + pi;
rho_temp(rho_temp > pi) = rho_temp(rho_temp > pi) - pi;

%Finding optical objects
opticalPolar2 = [];
for k = 1:n_grains
    idx = (abs(matched_coordinates(k, 1) - quiver_vectors_obj(1, :)') < 0.005) & ...
        (abs(matched_coordinates(k, 2) - quiver_vectors_obj(2, :)') < 0.005);

    row_temp = quiver_vectors_obj(:, idx)';
    opticalPolar2 = [opticalPolar2; row_temp];
end

%Converting U and V to inclination (a1) and azimuth (b1)
a = opticalPolar2(:, 3);
b = opticalPolar2(:, 4);
a1 = (180/pi)*acos(sqrt(a.^2 + b.^2)); %optical
b1 = (180/pi)*atan2(b, a);

a2 = (180/pi)*theta_temp; %ebsd
b2 = (180/pi)*rho_temp;

%Fitting lines
axisLimit1 = 90;
axisLimit2 = 180;

ft1 = fittype({'x'}); %forzed zero intercept
ft2 = fittype({'x', '1'});

x= a1'; %inclination
y= a2';
a_fit = linspace(0, axisLimit1, 10);
[p1_a, p1_gof_a] = fit(x', y', ft1);
[p2_a, p2_gof_a] = fit(x', y', ft2);
y1_fitted_a = feval(p1_a, a_fit);
y2_fitted_a = feval(p2_a, a_fit);

%str
coeffval_temp = coeffvalues(p1_a);
a_str_eq1 = sprintf('y = %0.2f*x, r-square = %0.2f', coeffval_temp(1), p1_gof_a.rsquare);

coeffval_temp = coeffvalues(p2_a);
a_str_eq2 = sprintf('y = %0.2f*x + %0.2f, r-square = %0.2f', ...
    coeffval_temp(1), coeffval_temp(2), p2_gof_a.rsquare);

x= b1'; %azimuth
y= b2';
b_fit = linspace(0, axisLimit2, 10);
[p1_b, p1_gof_b] = fit(x', y', ft1);
[p2_b, p2_gof_b] = fit(x', y', ft2);
y1_fitted_b = feval(p1_b, b_fit);
y2_fitted_b = feval(p2_b, b_fit);

%str
coeffval_temp = coeffvalues(p1_b);
b_str_eq1 = sprintf('y = %0.2f*x, r-square = %0.2f', coeffval_temp(1), p1_gof_b.rsquare);

coeffval_temp = coeffvalues(p2_b);
b_str_eq2 = sprintf('y = %0.2f*x + %0.2f, r-square = %0.2f', ...
    coeffval_temp(1), coeffval_temp(2), p2_gof_b.rsquare);

%Plot
nrows = 1;
ncolumns = 2;
d_offset = 0.01;
textAlignment = 'left';
fontSize = 14;
tl_pct = 0.03;
colorSel = {[1, 0, 0, 0.5], [0, 0, 1, 0.5]};
greySel = 0.8;

hFig = figure(12);
hFig.Position = [50, 50, 1500, 500];
clf('reset')

t = tiledlayout(nrows, ncolumns, TileSpacing= "tight");

nexttile(1)
ax = gca;

%fits
plot(a_fit, y1_fitted_a,'-', 'LineWidth', 2, 'Color', colorSel{1})
hold on
plot(a_fit, y2_fitted_a,'b-', 'LineWidth', 2, 'Color', colorSel{2})
%1:1 line
% plot([0, axisLimit1], [0, axisLimit1], '-', 'LineWidth', 1, 'Color', [0, 0, 0, 0.5])
%points
plot(a1, a2, '.black', 'MarkerSize', 15)
%labels
text(a1 + d_offset*axisLimit1, a2 - d_offset*axisLimit1, num2str(ebsdID2), ...
    'HorizontalAlignment', textAlignment, 'Color', 'black', 'FontSize', fontSize)
%fitted eq
text(axisLimit1*tl_pct, axisLimit1*(1-5*tl_pct), a_str_eq1, ...
    'HorizontalAlignment', textAlignment, 'Color', colorSel{1}(1:end-1), ...
    'FontSize', fontSize, 'BackgroundColor', [greySel, greySel, greySel])
text(axisLimit1*tl_pct, axisLimit1*(1-8*tl_pct), a_str_eq2, ...
    'HorizontalAlignment', textAlignment, 'Color', colorSel{2}(1:end-1), ...
    'FontSize', fontSize, 'BackgroundColor', [greySel, greySel, greySel])

hold off
xlim([0, axisLimit1])
ylim([0, axisLimit1])
xlabel('Optical method')
ylabel('SEM-EBSD mapping')
title('C-axis inclination (degrees)', 'Units', 'normalized', 'Position', [0.5, 0.92, 0])
grid on
% grid(ax,'minor')

nexttile(2)
ax = gca;

%fits
plot(b_fit, y1_fitted_b,'-', 'LineWidth', 2, 'Color', colorSel{1})
hold on
plot(b_fit, y2_fitted_b,'b-', 'LineWidth', 2, 'Color', colorSel{2})
%1:1 line
% plot([0, axisLimit2], [0, axisLimit2], '-', 'LineWidth', 1, 'Color', [0, 0, 0, 0.5])
%points
plot(b1, b2, '.black', 'MarkerSize', 15)
%labels
text(b1 + d_offset*axisLimit2, b2 - d_offset*axisLimit2, num2str(ebsdID2), ...
    'HorizontalAlignment', textAlignment, 'Color', 'black', 'FontSize', fontSize)
%fitted eq
text(axisLimit2*tl_pct, axisLimit2*(1-5*tl_pct), b_str_eq1, ...
    'HorizontalAlignment', textAlignment, 'Color', colorSel{1}(1:end-1), ...
    'FontSize', fontSize, 'BackgroundColor', [greySel, greySel, greySel])
text(axisLimit2*tl_pct, axisLimit2*(1-8*tl_pct), b_str_eq2, ...
    'HorizontalAlignment', textAlignment, 'Color', colorSel{2}(1:end-1), ...
    'FontSize', fontSize, 'BackgroundColor', [greySel, greySel, greySel])
hold off
xlim([0, axisLimit2])
ylim([0, axisLimit2])
xlabel('Optical method')
ylabel('SEM-EBSD mapping')
title('C-axis azimuth (degrees)', 'Units', 'normalized', 'Position', [0.5, 0.92, 0])
grid on
% grid(ax, 'minor')

%Hiding axes labels
[row, col] = tilerowcol(t.Children);
if (max(row) > 1)
    xticklabels(t.Children(row < t.GridSize(1)),"") 
    xlabel(t.Children(row < t.GridSize(1)),"")        
end
if (max(col) > 1)
    % yticklabels(t.Children(col > 1), "") %do not activate
    ylabel(t.Children(col > 1), "")
end

