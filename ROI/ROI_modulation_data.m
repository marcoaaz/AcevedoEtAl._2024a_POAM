
function ROI_modulation_data(roiHandle, img, color_space, range_modes, ...
    sel_modes, period_modality, pol_angle)

%Input:

%roiHandle = region of interest
%img = field of view cropped as 2000x2000x3x180 uint8
%color_space = selected RGB, HSV, CieLAB
%range_modes = data acquisition series indexes
%sel_modes = selected modality indexes
%period_modality = XPL=90, XPL-lambda & PPL=180
%pol_angle = nominal stage angle series (manual experiment). input as row

%%
%extracts logical mask
logicalMask = createMask(roiHandle);
sz = size(logicalMask);
ind = find(logicalMask); %column-major order
[row, col] = ind2sub(sz, ind);

n_pixels = length(row);
n_dataPoints = length(pol_angle);

%data
px_rgb_test = double(img(row, col, :, :));
px_rgb = squeeze(mean(px_rgb_test, [1, 2])); %integrating ROI pixels (3x180)

%Optional: change colour space
[px, channelNames] = switchColourSpace(px_rgb, color_space); %output [0-255]

%Optional: linearize
px_lin = rgb2lin(px'/255, 'ColorSpace', 'srgb'); %input [0-1]
px = 255*px_lin'; %output [0-255]

%% Fitting curves

x = pol_angle'; %as column   
n_modes = length(sel_modes);
n_channels = 3; 

fitCell = cell(1, n_modes);
xCell = cell(1, n_modes);
yCell = cell(1, n_modes);
periodCell = cell(1, n_modes);
legendCell = cell(1, n_modes);
k= 0;
for sel = sel_modes %[2, 3], ppl and xpl
    k = k + 1;

    %Axis            
    data = px(:, range_modes{sel}); %sRGB    
    
    %Fit model
    algorithm1_fit = fittype('a_0 + a_1*cos(x*w) + b_1*sin(x*w) + a_2*cos(2*x*w) + b_2*sin(2*x*w)');
    optical_period = period_modality(sel);        
    algorithm1_w = 2*pi/(optical_period); %radians

    xData = cell(1, n_channels + 1);
    yData = cell(1, n_channels + 1);
    fitresult = cell(1, n_channels + 1);
    gof = cell(1, n_channels + 1);
    period = cell(1, n_channels + 1);    
    legend_str = cell(1, n_channels + 1);
    
    %% Colour    
    for c = 1:n_channels %1:n_channels
        y = data(c, :);        
                
        %Fourier-2
        % [fitresult{c}, gof{c}, preparedData] = pixelFourierS(x, y, optical_period);
        % [fitresult{c}, gof{c}, preparedData] = pixelFourierS_ver2(x, y, optical_period);
        [fitresult{c}, gof{c}] = pixelFourierS_ver2(x, y, algorithm1_w, algorithm1_fit);
        
        [xData{c}, yData{c}] = prepareCurveData(x', y);
        
        w = fitresult{c}.w;
        period{c} = 2*pi/w; %period
        legend_str{c} = {sprintf('%s: %0.1f', channelNames{c}, period{c}); 
            sprintf('%s_fit', channelNames{c})};    
    end    

    %% Greyscale
    y = mean(data, 1); %Note: when moving ROI, 'y' is left incomplete [1, 2].
    c = n_channels + 1; %fourth channel   
          
    %Fourier-2
    % [fitresult{c}, gof{c}, preparedData] = pixelFourierS(x, y, optical_period);
    % [fitresult{c}, gof{c}, preparedData] = pixelFourierS_ver2(x, y, optical_period);
    [fitresult{c}, gof{c}] = pixelFourierS_ver2(x, y, algorithm1_w, algorithm1_fit);
    
    [xData{c}, yData{c}] = prepareCurveData(x', y);          
     
    w = fitresult{c}.w;
    period{c} = 2*pi/w; %period
    legend_str{c} = {sprintf('%s: %0.1f', 'greyscale' , period{c}); 
        sprintf('%s_fit', 'greyscale')};  
    
    % fitresult{4}

    %% pass to plotting loop
    
    fitCell{k} = fitresult;
    xCell{k} = xData;
    yCell{k} = yData;
    periodCell{k} = optical_period;
    legendCell{k} = cat(1, legend_str{:});

    %% Debug: save greyscale spectra
   
    %2 x range_modes (36 rotation steps)
    % writematrix(preparedData, 'ROI_modulation_data_TL-XPL-lambda_site1.csv')

end

roiHandle.UserData.n_dataPoints = n_dataPoints;
roiHandle.UserData.n_pixels = n_pixels;
roiHandle.UserData.fitRGB = fitCell;
roiHandle.UserData.xCell = xCell;
roiHandle.UserData.yCell = yCell;
roiHandle.UserData.periodCell = periodCell;
roiHandle.UserData.legendCell = legendCell;

end