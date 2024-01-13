
function [pixel_calc] = imageFourierS_optim(img_temp2, pol_angle1, optical_period, cluster)

%pol_angle (36x1) = nominal stage angle series (manual experiment). 
% Value in sexagesimal degrees. 

%optical_period = 90 for PPL; 180 for XPL and XPL-lambda plate

imgDim = size(img_temp2);
n_rows1 = imgDim(1);
n_cols1 = imgDim(2);

pol_angle1 = pol_angle1'; %test

algorithm1_fit = fittype('a_0 + a_1*cos(x*w) + b_1*sin(x*w) + a_2*cos(2*x*w) + b_2*sin(2*x*w)');
algorithm1_w = 2*pi/(optical_period); %radians
  
pixel_avg = zeros(n_rows1, n_cols1, 'double');
pixel_phase = zeros(n_rows1, n_cols1, 'double');
pixel_range1 = zeros(n_rows1, n_cols1, 'double');
pixel_range2 = zeros(n_rows1, n_cols1, 'double');
pixel_maxPeak = zeros(n_rows1, n_cols1, 'double');
pixel_maxPeakI = zeros(n_rows1, n_cols1, 'double');
pixel_SSE = zeros(n_rows1, n_cols1, 'double');

for i = 1:n_cols1 %934 %
% for i = 400:430 %934 %
% parfor (i = 1:n_cols1, cluster)
    sprintf('column: %d', i)

    px_col = img_temp2(:, i, :);

    col_avg = zeros(n_rows1, 1, 'double');
    col_phase = zeros(n_rows1, 1, 'double');
    col_range1 = zeros(n_rows1, 1, 'double');
    col_range2 = zeros(n_rows1, 1, 'double');    
    col_maxPeak = zeros(n_rows1, 1, 'double'); 
    col_maxPeakI = zeros(n_rows1, 1, 'double');
    col_SSE = zeros(n_rows1, 1, 'double');
    parfor (j = 1:n_rows1, cluster)            
%     for j = 1177%1:n_rows1

        px = squeeze(px_col(j, 1, :)); %36x1
        
        % [fitresult, gof] = pixelFourierS(pol_angle1, px, optical_period); %Algorithm 1
        [fitresult, gof] = pixelFourierS_ver2(pol_angle1, px, algorithm1_w, algorithm1_fit);

        s = sinDescriptor(fitresult); %Algorithm 2     
%         s.approx %debug        

        col_avg(j) = s.avg_f;
        col_phase(j) = s.avg_x;
        col_range1(j) = s.range1;        
        col_range2(j) = s.range2;
        col_SSE(j) = gof.sse;
        
        %Saving only one peak (the max. of the two peaks found in greyscale)
        B2 = s.MinMax;
        max_x = B2(1, 1);
        max_f = B2(1, 2);
        col_maxPeak(j) = max_x; 
        col_maxPeakI(j) = max_f; %continuous estimate             
    end
    pixel_avg(:, i) = col_avg; %R
    pixel_phase(:, i) = col_phase; %G
    pixel_range1(:, i) = col_range1; %B
    pixel_range2(:, i) = col_range2;    
    pixel_maxPeak(:, i) = col_maxPeak; %stage angle for c-axis estimate
    pixel_maxPeakI(:, i) = col_maxPeakI; 
    pixel_SSE(:, i) = col_SSE; %fitting, sum of squared errors    
end

pixel_calc = cat(3, pixel_avg, pixel_phase, pixel_range1, ...
    pixel_range2, pixel_maxPeak, pixel_maxPeakI, pixel_SSE);

end