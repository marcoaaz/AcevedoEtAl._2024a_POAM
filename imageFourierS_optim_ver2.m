
function [imgSpectraModel] = imageFourierS_optim_ver2(img_temp2, pol_angle1, optical_period, quality_mask, cluster)

%Input:
%img_temp2 (n_rows1, n_cols1, n_layers) = downscaled spectral cube (1 channel)
%pol_angle (36x1) = nominal stage angle series (manual experiment). Value in sexagesimal degrees. 
%optical_period = 90 for PPL; 180 for XPL and XPL-lambda plate
%quality_mask = fraction of data completeness [0-1]
%cluster = number of logical processors to use (find sweet spot)

%Output:
%imgSpectraModel = modelled and described spectra

%%

imgDim = size(img_temp2);
n_rows1 = imgDim(1);
n_cols1 = imgDim(2);
n_layers = imgDim(3);
quality_ratio = 0.8; %fraction of non-zero values in pixels

% pol_angle1 = pol_angle1'; %(1x36)

%Algorithm 1 (Fourier-2)
% algorithm1_fit = fittype('a_0 + a_1*cos(x*w) + b_1*sin(x*w) + a_2*cos(2*x*w) + b_2*sin(2*x*w)');
algorithm1_w = 2*pi/(optical_period); %radians

%Coefficients
angle_temp = reshape(pol_angle1*algorithm1_w, 1, 1, []); %36x1
yData_temp1 = img_temp2*(2/n_layers);%(n_rows1, n_cols1, 36)

a_0 = mean(img_temp2, 3); %(n_rows1, n_cols1, 1)
a_1 = sum(cos(angle_temp).*yData_temp1, 3);
b_1 = sum(sin(angle_temp).*yData_temp1, 3);
a_2 = sum(cos(2*angle_temp).*yData_temp1, 3);
b_2 = sum(sin(2*angle_temp).*yData_temp1, 3);
model_temp = cat(3, a_0, a_1, b_1, a_2, b_2); %~100MB

%Algorithm 2
search_interval = 30; 
%Note: min. expected angular distance between features (TL-XPL-lambda)
x0_range = 0:search_interval:optical_period-search_interval;

pixel_avg = zeros(n_rows1, n_cols1, 'double');
pixel_phase = zeros(n_rows1, n_cols1, 'double');
pixel_range1 = zeros(n_rows1, n_cols1, 'double');
pixel_range2 = zeros(n_rows1, n_cols1, 'double');
pixel_maxPeak = zeros(n_rows1, n_cols1, 'double');
pixel_maxPeakI = zeros(n_rows1, n_cols1, 'double');
% pixel_SSE = zeros(n_rows1, n_cols1, 'double');

tic;
time1 = 0;
for i = 1:n_cols1 %real

% for i = 760:960 %934 %
% for i = 760% debug
% parfor (i = 1:n_cols1, cluster)    
    
    %Parallel loop
    model_col = model_temp(:, i, :);
    mask_col = quality_mask(:, i);

    col_avg = zeros(n_rows1, 1, 'double');
    col_phase = zeros(n_rows1, 1, 'double');
    col_range1 = zeros(n_rows1, 1, 'double');
    col_range2 = zeros(n_rows1, 1, 'double');    
    col_maxPeak = zeros(n_rows1, 1, 'double'); 
    col_maxPeakI = zeros(n_rows1, 1, 'double');    

    parfor (j = 1:n_rows1, cluster) %real

    % for j = 1341%1:n_rows1 %debug
    % for j = 1:n_rows1

        % %debug
        % sprintf('row= %d', j)

        model_px = squeeze(model_col(j, 1, :)); %5x1        
        mask_px = mask_col(j, 1);

        if mask_px > quality_ratio
            
            % Algorithm 1            
            x = sym('x');

            f = (model_px(1) + ...
                model_px(2)*cos(x*algorithm1_w) + model_px(3)*sin(x*algorithm1_w) + ...
                model_px(4)*cos(2*x*algorithm1_w) + model_px(5)*sin(2*x*algorithm1_w));
            f_d1 = diff(f, x);
            f_d2 = diff(f_d1, x);
            %Note: In you input a constant function, f_d1 becomes zero and the code
            %wont run properly. Thus, fitted background (0 cfit) cannot be an input.            
            
            %Function handles for x25 speed boost in Algorithm 2
            F = matlabFunction(f); 
            F_d1 = matlabFunction(f_d1);
            F_d2 = matlabFunction(f_d2);
                        
            % Algorithm 2 (computes 10 times longer than Algorithm 1)               
            s = sinDescriptor_ver3(x0_range, F, F_d1, F_d2);

            % %debug section    
            % xData = pol_angle1;
            % yData = squeeze(img_temp2(j, i, :));
            % sinDescriptorPlot_ver2(f, s, xData, yData, optical_period)           
            
            B2 = s.MinMax;        
            col_avg(j) = s.avg_f;
            col_phase(j) = s.avg_x;
            col_range1(j) = s.range1;        
            col_range2(j) = s.range2;                      
            col_maxPeak(j) = B2(1, 1); %max_x; Max. of the two max peaks found in greyscale  
            col_maxPeakI(j) = B2(1, 2); %max_f (continuous estimate)                  
           
        else
            col_avg(j) = 0;
            col_phase(j) = 0;
            col_range1(j) = 0;        
            col_range2(j) = 0;                      
            col_maxPeak(j) = 0; 
            col_maxPeakI(j) = 0;
        end
    end
    pixel_avg(:, i) = col_avg; %R
    pixel_phase(:, i) = col_phase; %G
    pixel_range1(:, i) = col_range1; %B
    pixel_range2(:, i) = col_range2;    
    pixel_maxPeak(:, i) = col_maxPeak; %stage angle for c-axis estimate
    pixel_maxPeakI(:, i) = col_maxPeakI; 
    
    %Stopwatch steps
    time2 = toc;
    time_interval = time2 - time1;
    time1 = time2;
    sprintf('column= %d, elapsed time= %.1f, interval time= %.1f', ...
        i, time2, time_interval)
end
% Save model 
imgSpectraModel = struct;            
imgSpectraModel.pixel_avg = pixel_avg;
imgSpectraModel.pixel_phase = pixel_phase;
imgSpectraModel.pixel_range1 = pixel_range1;
imgSpectraModel.pixel_range2 = pixel_range2;
imgSpectraModel.pixel_maxPeak = pixel_maxPeak; %stage rotation for max peak [0-180> 
imgSpectraModel.pixel_maxPeakI = pixel_maxPeakI;

end