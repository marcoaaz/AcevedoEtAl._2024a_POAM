
%Copyright: Marco Acevedo Zamora, 
%Written, M.A.: 8-Aug-2022
%Updated: M.A.: 14 april 2023; 8 may 2023; 4 jan 2024; 

% stack_spectra_leica_v17.m script:
%
% Reads on the registered image stack (from TrakEM2), it models stage
% rotation TL/RL-PPL/XPL/XPL-lambda waves, produces optic-axis orientation
% maps (with algorithm 3 or 4), and does ray tracing, image blending, and
% plot object-based image analysis (OBIAS) object quiver plots.

close all
clear 
clc

%dependency for stereonets
mtex_folder = 'E:\Alienware_March 22\current work\00-new code May_22\mtex-5.9.0';
cd(mtex_folder);
startup_mtex

%Dependencies
scriptDir = 'E:\Alienware_March 22\current work\00-new code May_22';
scriptDir2 = 'E:\Alienware_March 22\scripts_Marco\updated MatLab scripts\';
scriptDir3 = fullfile(scriptDir2, 'ROI');

addpath(scriptDir)
addpath(fullfile(scriptDir, 'rayTracing/'))     
addpath(fullfile(scriptDir2, 'ROI'))
addpath(fullfile(scriptDir2, 'plots_miscellaneous'))

%Only for Algorithm 4: data directory for Michel-Levy modelling
sourceFolder_Michel = 'E:\paper 2_datasets\nikon LV1000ND\M-L_recolouring';

%Notes:
%Paper datasets (algorithm 3)
%E:\paper 2_datasets\nikon LV1000ND\mylonite\final stacks
%E:\paper 2_datasets\nikon LV1000ND\mylonite\zoom-in registeredStacks
%E:\paper 2_datasets\nikon LV1000ND\cpx\stacked\saved_aligned\stack_modality
%E:\paper 2_datasets\nikon LV1000ND\granulite 6KB-67\composite
%E:\paper 2_datasets\nikon LV1000ND\TS-J2KB\stacks\composite

%setting examples:
%plotSetup.lineWidth = 0.6; %mylonite and garnet
%plotSetup.autoSFactor = 0.4; %mylonite and garnet
%plotSetup.edgeAlpha = 1; %coarse grained mineral

%harzburgite 035
% plotSetup.lineWidth = 0.001; %vector thickness
% plotSetup.autoSFactor = .2; %quiver vectors size
% plotSetup.edgeAlpha = 1; %edge transparency

%Paper datasets (algorithm 4)
%E:\paper 2_datasets\nikon LV1000ND\harzburgites Carl\02\composite
%E:\paper 2_datasets\nikon LV1000ND\harzburgites Carl\03\composite
%E:\paper 2_datasets\nikon LV1000ND\17BSK035\stacks_save\composite_trial4

%% Section 1: User input

%Dataset definition
workingDir = 'E:\paper 2_datasets\nikon LV1000ND\harzburgites Carl\02\composite';
fileName = 'all_modalities.tif'; %TrakEM2 image stack
blendedImage_suffix = 'test2'; %SuperSIAT input
fileName_spf = 'sample7.spf'; %SuperSIAT training points
fileName_shp = 'sample7_segmentation.shp'; %objects
destDir_suffix = 'linear';%saving directory (change to avoid overwriting)
plotSetup = struct;
plotSetup.lineWidth = 0.05;
plotSetup.autoSFactor = .5;
plotSetup.edgeAlpha = 0;

%plotSetup.lineWidth = 0.6; %mylonite and garnet
%plotSetup.autoSFactor = 0.4; %mylonite and garnet
%plotSetup.edgeAlpha = 1; %coarse grained mineral

%Define: Multi-pol imaging spectra
sel_modality = {...
    'RL PPL', 'RL XPL', ...
    'TL PPL', 'TL XPL', 'TL XPL-lambda'}; %acquired series
rlPPL_range = 1:36; 
rlXPL_range = 37:72;
tlPPL_range = 73:108;
tlXPL_range = 109:144;
tlXPLlambda_range = 145:180;
%Note: fitting a Fourier-2 model requires at least 6 data points

%Define ray tracing
% {'mean', 'max', 'min', 'range', 'sum', 'std', 'median', 'maxHSV', 'minHSV', 'rangeHSV'}
stats_list = {'max', 'min'};

%Define image blending
type = {'max'}; %used ray tracing method

% %Mylonite
% modality = {'RL PPL', 'TL PPL', 'TL XPL', ...
%     'rescaled8_pca_25april_pctOut0.5_denoised3x3_registered'}; %'RL XPL', 'TL PPL', 'TL XPL-lambda'
% weights = [15, 20, 20, 45];
%Every sample (not the harzburgites)
modality = {'RL PPL', 'TL PPL', 'TL XPL'}; %'RL XPL', 'TL PPL', 'TL XPL-lambda'
weights = [50, 25, 25];

%Define Algorithm 1
pol_angle1 = [ 0:10:350 ]'; % experiment range (transposed to match sequezeed px)
period_modality = [180, 90, 180, 90, 180]; % optical_period; XPL=90, PPL=180
sel = 5; %select idx of 'TL-XPL-lambda' modality 
scale = 0.5; % Decreases computational cost
cluster = 10; %logical processors >8 crashed 32GB (dont forget to turn off debugging plots!)

importType_loop = 1; %greyscale = 1; rgb = 2
n_RGB = 1;
% sel_channel = [1, 2, 3]; %separate the channels
% n_RGB = length(sel_channel);

%Saving file name of Algorithms 1 & 2 loop
descripTxt = {'6-Jan-24'}; %name of modulation image file
% descripTxt = strcat({'R', 'G', 'B'}, '_'); %when looping each RGB channel

%Define POAM maps
mappingAlgorithm = 4; %3= Algorihtm 3; 4= Algorithm 4
fileName_equation = 'EQst_CieLAB_17-BSK-062-02_linearised.m'; %_linearised
obias_folderName = 'obiasExport';

%% Section 2: Scripted image processing
%Importing single image stacks (patch)

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

%Pre-check number of layers
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

%Section 1: Ray tracing

%Informative structure
info_struct.Height = n_rows; 
info_struct.Width = n_cols;
info_struct.Channels = n_channels;
info_struct.sel_range = sel_range;
info_struct.sel_modality = sel_modality; %str

n_options = length(stats_list);
time_elapsed = zeros(1, n_options);
for k = 1:n_options     
    %parallel computing (default: commented)
%     [time_elapsed] = stats_zProject(imageName, info_struct, stats_list{k}, destFolder_rt);
end

%Section 2: Image Blending of ray tracing images (OBIAS input)
n_mode = length(modality);
strMode = [];
for j = 1:n_mode
    strMode = strcat(strMode, '_', modality{j});
end

destFile_blended = fullfile(destFolder_rt, strcat('imgSum', strMode, '_', blendedImage_suffix, '.tif')); %test1_W5000
saveOption = 0; %no need to save twice
[img_sum_w_rs] = stack_imageFusion(type, modality, weights, destFile_blended, saveOption);

%Defining SuperSIAT directories after OBIAS (object-based image analysis)
workingDir_obias = strrep(destFile_blended, '.tif', ''); %SuperSIAT folder name convention
fileName_TIF = strrep(fileName_shp, '.shp', '.tif'); %labelled map
shapeFile_dir = fullfile(workingDir_obias, fileName_spf); %shape file

%% Section 3 (optional): ROI Tool in live mode

img_A = img_sum_w_rs;

%Edit ROI: 3000x3000x3x180 uint8 is 5GB and slow
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
sel_modes = [1:5]; %
addlistener(hCircle, 'MovingROI',...
    @(varargin)ROI_modulation_data(hCircle, img_full_temp, ...
    color_space, sel_range, sel_modes, period_modality, pol_angle1')); %'MovingROI'

% sel_modes = [3, 4]; %for paper (Fig aspect was [1, 2])
sel_modes = [5]; %1:n_modalities
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROI_modulation_graph(hCircle, sel_modes, sel_modality)); 
%sel_modality = {'RL PPL', 'RL XPL', 'TL PPL', 'TL XPL', 'TL XPL-lambda'};

selected_modality = 5;
selected_channel = 4; %R, G, B, Greyscale
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROI_modulationAlgorithm_graph(hCircle, ...
    selected_modality, selected_channel)); 

%% Section 4: Modulation image calculation (spectral modelling)

%channel
for uu = 1:n_RGB %uu = 1:n_RGB
    
    %modality 
    for sel = 5 % for sel = 1:5
        optical_period = period_modality(sel);
    
        % Preparing image stack (spectra)
        varRange = sel_range{sel};                
        [img_temp2, quality_mask2] = stackImportLoop(imageName, ...
            info_struct, varRange, importType_loop, scale);        

        tic;

        %Algorithm 1 and 2 (parallel loop)                          
        [imgSpectraModel] = imageFourierS_optim_ver2(img_temp2, pol_angle1, ...
            optical_period, quality_mask2, cluster);

        %Save progress
        file_output = strcat(descripTxt{uu}, ...
            'modulationImage_', sel_modality{sel}, '_', num2str(scale), '.mat'); %'wave_TL-PPL_0.15.mat'
        destFile = fullfile(destFolder, file_output);
        save(destFile, "imgSpectraModel", '-mat', '-v7.3')
        
        time = toc;
        t_spectra.FittedImage(uu) = time/3600; %hours
        t_spectra.Rate(uu) = n_pixels1/time; %px/sec           
    
    end %remove when only Section 2 required
end
t_spectra

%% Section 5 (optional): Modulation RGB image (requires Section 2) 
%Build modulation image. Following Axer et al., 2011

uu = 1; %channel
sel = 5;

%edit descriptTxt (when working with channel names)
file_output2 = strcat(descripTxt{uu}, 'modulationImage_', ...
    sel_modality{sel}, '_', num2str(scale), '.mat'); 
%old: 'wave2_', new: '{30Dec}_modulationImage_'

imageFile1 = strrep(file_output2, '.mat', '.tif'); %false coloured modulation image
imageFile2 = strrep(imageFile1, '.tif', '_bgMask.tif');

m = load(fullfile(destFolder, file_output2));    
R_double = m.imgSpectraModel.pixel_avg; %transmittance
G_double = m.imgSpectraModel.pixel_phase; %direction (inclination proxy)
temp_B_double = m.imgSpectraModel.pixel_range1;
B_double = temp_B_double./R_double; %retardation (division by 0 = NaN)

%Mask of background (below 'quality ratio')    
mask_bg = (R_double == 0) & (G_double == 0) & isnan(B_double); 
%Note: it can change cropping coordinates within a ~6 pixel radius and
%therefore the object centroid x-y 'accuracy' (in POAM validation step)

R_double(mask_bg) = 0;
G_double(mask_bg) = 0;
B_double(mask_bg) = 0;
R = uint8(rescale(R_double, 0, 255)); 
G = uint8(rescale(G_double, 0, 255)); 
B = uint8(rescale(B_double, 0, 255)); 
rgb_summary = cat(3, R, G, B);

imshow(rgb_summary)
imwrite(rgb_summary, fullfile(destFolder, imageFile1), 'compression', 'none')
imwrite(mask_bg, fullfile(destFolder, imageFile2), 'compression', 'none')
    
%% Section 6: Importing blended image and mineral grain objects 

%Import fused image
img_fullRes = imread(destFile_blended); %full (or partial 5000x5000)
%greyscale image: can inform about texture w/ translucent OBIAS map
img_grey = cat(3, rgb2gray(img_fullRes), ...
    rgb2gray(img_fullRes), ...
    rgb2gray(img_fullRes)); %repetition required to preserve cmap


%Obtaining SuperSIAT v2.2 map information (follows version convention)
[annotationNames, cmap] = SuperSIATinterpreter(shapeFile_dir); %function
annotationNames' %check

%Object structure 
S = shaperead(fullfile(workingDir_obias, fileName_shp));
range_ID = unique([S.ClassID]);

%Label map (32-bit), objects and classes
label_maps_SS = imread(fullfile(workingDir_obias, fileName_TIF));
label_map = double(label_maps_SS(:, :, 2)) + 1; %to avoid indexing issues 
[n_rows_partial, n_cols_partial] = size(label_map);

%scale factor: 
% orientation maps (downscaled) / segmentation (full or partial resolution 5000x5000) 
scale_found = n_rows1/n_rows_partial; 

%Creating Saving Folder
saveFolder = strcat(obias_folderName, '_', num2str(scale), '_', destDir_suffix);
saveDir = fullfile(destFolder_rt, saveFolder);
mkdir(saveDir)

%% Section 7: Target C-axis estimation (discrete estimate)
close all

sel = 5;%modality: TL-XPL-lambda
mineralSel = 3; %selecting target mineral
Imax_vector = [135; 135]; %max intensity found in [3D universal stage or 2D section]
mineralText = annotationNames(mineralSel);


%Configuring 

% %mylonite quartz
% filter_expression = "([stats_add.foreground]' == 1 & [stats_add.Solidity]' > .35 & [stats_add.Circularity]' > .15)";

%default (trial and error approach)
filter_expression = "([stats_add.foreground]' == 1)";

%Note: The object filtering criteria follows the quality control in the maps at the
%end of the script. This expression enters 'targetObjectQuiverP.m'.

%Background mask (within all_modalities_falseColor dir)
mask_bg = imread(fullfile(destFolder, imageFile2)); 
mask_bg1 = imresize(mask_bg, ...
    [n_rows_partial, n_cols_partial]); %upscaled bg (full resolution)
%Note: Careful. The upscaled background defines the cropping ROI (see below).

%Foreground mask: downscaling to match maps (assumes FoV are equal)
label_map_ds = imresize(label_map, [n_rows1, n_cols1], "nearest"); %NaN finds automatically
mask_fg = (label_map_ds == mineralSel) & (~mask_bg); %downscaled foreground (mineral inside)

%Importing image stack (3 min): optional in Algorithm 3 & mandatory in Algorithm 4
sel = 5; %TL-XPL-lambda 
varRange = sel_range{sel};

switch mappingAlgorithm
    
    case 3 % Algorithm 3

        %Ranges
        pixel_range1 = m.imgSpectraModel.pixel_range1;        
        
        %Peak 1
        pixel_maxPeak = m.imgSpectraModel.pixel_maxPeak; %continuous max 'peak 1' stage rotation value
        [pixel_cAxis] = stageToAxisReorientation(pixel_maxPeak); 
        
        % %Optional: closest discrete value        
        % importType = 1; %greyscale
        % [img_temp2, ~] = stackImportLoop(imageName, info_struct, varRange, importType, scale);        
        % [pixel_maxPeakI_discrete] = largestPeak_channels(...
        %     img_temp2, pol_angle1, pixel_maxPeak); 
        
        [theta_inclination, pixel_cAxisI_discrete_rs] = estimateInclination_range(...
            pixel_range1, mask_fg, Imax_vector);
        
        %Save csv for Stereonet (only for downscaled output) 
        saveCSV = 0; %yes/no
        destFile_fullStereo = fullfile(saveDir, strrep(file_output2, '.mat', '_fullNet.csv'));
        if saveCSV == 1    
            text_save = [pixel_cAxis(:), theta_inclination(:)];
            text_save(text_save(:, 1) == -1, :) = [];
            writematrix(text_save, destFile_fullStereo)   
        end
        
        %Optional: noise reduction
        r = 1;
        kernel1 = [2*r+1, 2*r+1];
        theta_inclination2 = medfilt2(theta_inclination, kernel1);
        pixel_cAxis2 = medfilt2(pixel_cAxis, kernel1); 
        %pixel_cAxis2 (algorithm 3) is equivalent to cAxis_peak1_med(algorithm 4)
        
        
        %POAM maps
        %plot and return downscaled and cropped label map
        plotOption = 1;
        saveOption = 1;
        newTxt = strcat('_pxNet_', mineralText, '.csv');
        destFile_targetGridStereo = fullfile(saveDir, strrep(file_output2, '.mat', newTxt));
        
        [img_cAxis_rgb, mask_mineral_crop, pos_ROI_ds] = plotTargetAzimuthMap(...
            mask_fg, mask_bg, pixel_cAxis2, theta_inclination2, plotOption, ...
            destFile_targetGridStereo, saveOption); %image as double
        %Notes: requires MTEX
        
        %upscaling, cropping, and making overlay of orientation map (for quiver plots)
        img_cAxis_rgb2 = uint8(imresize(img_cAxis_rgb, 1/scale_found, 'bilinear'));
        [pos_ROI, img_fused] = imgCropAndOverlay(mask_bg1, img_grey, img_cAxis_rgb2); %option 1
        % img_fused = imfuse(img_fullRes_crop, img_cAxis_rgb2, 'blend'); %option 2 (overwrite img_fused)
        
        % %Optional: calculate orientation map vectors ('grid quiver plot')
        % mask_bg1_crop = mask_bg1(pos_ROI(1):pos_ROI(3), pos_ROI(2):pos_ROI(4));
        % [quiver_vectors_grid] = targetGridQuiverP(img_fullRes_crop, mask_bg1_crop, ...
        %     pixel_cAxis2, theta_inclination2, mask_fg, scale_found, pos_ROI);
        
        %Target object calculations (2 min) 
        %filter expression = "([stats_add.foreground]' == 1)"; %default
        S_sub = S([S.ClassID] == mineralSel - 1); %-1 fixes indexing issue
        plotOption = 0;
        saveOption = 1;
        newTxt = strcat('_objNet_', mineralText, '.csv'); %Saving filtered data
        destFile_ori_target = fullfile(saveDir, strrep(file_output2, '.mat', newTxt));
        
        [S_sub_filtered, stats_add_filtered, quiver_vectors_obj...
            ] = targetObjectQuiverP(...
            img_fullRes, S_sub, pixel_cAxis2, theta_inclination2, ...
            mask_fg, scale_found, pos_ROI, filter_expression, ...
            plotOption, destFile_ori_target, saveOption);         
        
        %Save unfiltered object orientation vectors (used in quiver plot below)
        destFile_objectQuiver = strrep(destFile_ori_target, ...
            '.csv', '_centroidAndVector.csv'); %for validation with EBSD
        writematrix(quiver_vectors_obj, destFile_objectQuiver)
    
    case 4 % Algorithm 4

        %Ranges
        pixel_range1 = m.imgSpectraModel.pixel_range1;
        pixel_range2 = m.imgSpectraModel.pixel_range2; %used in algorithm 4
        
        importType = 2;
        [img_temp2_channels, ~] = stackImportLoop(imageName, ...
            info_struct, varRange, importType, scale);   
        
        %stage rotation angle when max peak
        pixel_maxPeak1 = m.imgSpectraModel.pixel_maxPeak;% peak 1   
        pixel_maxPeak2 = pixel_maxPeak1 + 90; %peak 2
        pixel_maxPeak2(pixel_maxPeak2 > 180) = pixel_maxPeak2(pixel_maxPeak2 > 180) - 180;

        [pixel_cAxis1] = stageToAxisReorientation(pixel_maxPeak1); 
        [pixel_cAxis2] = stageToAxisReorientation(pixel_maxPeak2);         
        
        %Colour adaptation estimation with polynomial fitting (saving for adaptingColors.m)
        file_peak1 = fullfile(workingDir, 'intensity_peak1.tif');
        file_peak2 = fullfile(workingDir, 'intensity_peak2.tif');
        if isfile(file_peak1) && isfile(file_peak2)

            pixel_maxPeakI_discrete1 = imread(file_peak1);
            pixel_maxPeakI_discrete2 = imread(file_peak2);        
        else
            %closest discrete value
            [pixel_maxPeakI_discrete1] = largestPeak_channels(...
                img_temp2_channels, pol_angle1, pixel_maxPeak1);        
            [pixel_maxPeakI_discrete2] = largestPeak_channels(...
                img_temp2_channels, pol_angle1, pixel_maxPeak2);
        
            imwrite(uint8(pixel_maxPeakI_discrete1), file_peak1, 'compression', 'none')
            imwrite(uint8(pixel_maxPeakI_discrete2), file_peak2, 'compression', 'none')
        end

        %Colour matching        
        file_closest1 = fullfile(workingDir, 'img_closest1.tif');
        file_closest2 = fullfile(workingDir, 'img_closest2.tif');
        if isfile(file_closest1) && isfile(file_closest2)

            img_closest1 = imread(file_closest1);
            img_closest2 = imread(file_closest2);        
        else
            tic;        
            %run only once
            chartFile = fullfile(sourceFolder_Michel, 'M-L chart.tif'); %sRGB
            modelFile = fullfile(sourceFolder_Michel, fileName_equation); %equation
            %requisite: running 'adaptingColors.m' script to obtain modelFile
        
            plotOption = 0; %1 when processing one pixel during QC & debugging
            [img_closest1, img_closest2] = highestOrderPeak_ver2(...
                pixel_maxPeakI_discrete1, pixel_maxPeakI_discrete2, ...
                modelFile, chartFile, workingDir, plotOption);
            
            [n_rows_ds, n_cols_ds] = size(img_closest1);
            time_chart = toc;
            t_finding = time_chart/60; %min
            t_finding_px = (n_rows_ds*n_cols_ds)/time_chart; %56K px/sec 
            sprintf('Colour matching step. Total time = %.1f min, time per pixel = %.1f', ...
                t_finding, t_finding_px)
        
        end

        %Rearrange peak images
        [cAxis_peak1, cAxis_peak2, ...
            pixel_peak1, pixel_peak2, ...
            range_peak1, range_peak2] = rearrangePeakImages(...
            img_closest1, img_closest2, ...
            pixel_cAxis1, pixel_cAxis2, ...
            pixel_maxPeakI_discrete1, pixel_maxPeakI_discrete2, ...
            pixel_range1, pixel_range2);
        
        range_peak1_grey = mean(range_peak1, 3); 
        [theta_inclination, pixel_cAxisI_discrete_rs] = estimateInclination_range(...
            range_peak1_grey, mask_fg, Imax_vector);
        %Note: the slow-axis does not have inclination (this only contributes to
        %improve the maps)
        %range_peak1 is sRGB (not linear) until recalculating the orientation maps.

        %Save csv for Stereonet (only for downscaled output) 
        saveCSV = 1; %yes/no
        destFile_fullStereo = fullfile(saveDir, strrep(file_output2, '.mat', '_fullNet.csv'));
        if saveCSV == 1    
            text_save = [cAxis_peak1(:), theta_inclination(:)];
            text_save(text_save(:, 1) == -1, :) = [];
            writematrix(text_save, destFile_fullStereo)   
        end
        
        %Optional: noise reduction
        r = 1;
        kernel1 = [2*r+1, 2*r+1];
        cAxis_peak1_med = medfilt2(cAxis_peak1, kernel1);
        theta_inclination2 = medfilt2(theta_inclination, kernel1);        
        
        %POAM maps
        %plot and return downscaled and cropped label map
        plotOption = 1;
        saveOption = 1; %Notes: requires MTEX
        newTxt = strcat('_pxNet_', mineralText, '.csv');
        destFile_targetGridStereo = fullfile(saveDir, strrep(file_output2, '.mat', newTxt));
        
        [img_cAxis_rgb, mask_mineral_crop, pos_ROI_ds] = plotTargetAzimuthMap(...
            mask_fg, mask_bg, cAxis_peak1_med, theta_inclination2, plotOption, ...
            destFile_targetGridStereo, saveOption); %image as double        
        
        %upscaling, cropping, and making overlay of orientation map (for quiver plots)
        img_cAxis_rgb2 = uint8(imresize(img_cAxis_rgb, 1/scale_found, 'bilinear'));
        [pos_ROI, img_fused] = imgCropAndOverlay(mask_bg1, img_grey, img_cAxis_rgb2); %option 1
        % img_fused = imfuse(img_fullRes_crop, img_cAxis_rgb2, 'blend'); %option 2 (overwrite img_fused)
        
        % %Optional: calculate orientation map vectors ('grid quiver plot')
        % mask_bg1_crop = mask_bg1(pos_ROI(1):pos_ROI(3), pos_ROI(2):pos_ROI(4));
        % [quiver_vectors_grid] = targetGridQuiverP(img_fullRes_crop, mask_bg1_crop, ...
        %     cAxis_peak1_med, theta_inclination2, mask_fg, scale_found, pos_ROI); 
        
        %Target object calculations (2 min) 
        %filter expression = "([stats_add.foreground]' == 1)"; %default
        S_sub = S([S.ClassID] == mineralSel - 1); %-1 fixes indexing issue
        plotOption = 1;
        saveOption = 1;
        newTxt = strcat('_objNet_', mineralText, '.csv'); %Saving filtered data
        destFile_ori_target = fullfile(saveDir, strrep(file_output2, '.mat', newTxt));
        
        [S_sub_filtered, stats_add_filtered, quiver_vectors_obj...
            ] = targetObjectQuiverP(...
            img_fullRes, S_sub, cAxis_peak1_med, theta_inclination2, ...
            mask_fg, scale_found, pos_ROI, filter_expression, ...
            plotOption, destFile_ori_target, saveOption);        
        
        %Save unfiltered object orientation vectors (used in quiver plot below)
        destFile_objectQuiver = strrep(destFile_ori_target, ...
            '.csv', '_centroidAndVector.csv'); %for validation with EBSD
        writematrix(quiver_vectors_obj, destFile_objectQuiver)
end
%Note: pixel_cAxis2 (algorithm 3) is equivalent to cAxis_peak1_med(algorithm 4)

%% Section 8: Plot OBIAS map & object optic-axis orientation quiver map

% %Optional: trial other settings
% plotSetup.lineWidth = 0.3; %0.001 vector thickness
% plotSetup.autoSFactor = .4; %0.2 quiver vectors size

plotSetup.cmap = cmap;
plotSetup.range_ID = range_ID; 

plotSetup.transparentVal = 1; %SuperSIAT class map colour strength (default=0.3)
plotOBIASmap(img_fullRes, S, plotSetup, annotationNames, saveDir);

%img_fused or img_gray_crop
plotSetup.transparentVal = 0;
plotTargetObjectQuiver(img_fused, S_sub_filtered, quiver_vectors_obj,  ...
    plotSetup, mineralText, saveDir)

annotationNames' %check
