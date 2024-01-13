function [pixel_avg, pixel_range, pixel_phase, pixel_max_x, pixel_SSE] = imageFourierS(img_temp2, pol_angle1, optical_period, cluster)

imgDim = size(img_temp2);
n_rows1 = imgDim(1);
n_cols1 = imgDim(2);
  
% pixelModels = cell(n_rows, n_cols);
pixel_avg = zeros(n_rows1, n_cols1, 'double');
pixel_range = zeros(n_rows1, n_cols1, 'double');
pixel_phase = zeros(n_rows1, n_cols1, 'double');
pixel_max_x = zeros(n_rows1, n_cols1, 'double');
pixel_SSE = zeros(n_rows1, n_cols1, 'double');

parfor (i = 1:n_cols1, cluster)
%         for i = 1:n_cols1
%     for i = 690
    col_avg = zeros(n_rows1, 1, 'double');
    col_range = zeros(n_rows1, 1, 'double');
    col_phase = zeros(n_rows1, 1, 'double');
    col_max_x = zeros(n_rows1, 1, 'double'); %stage angle for c-axis estimate
    col_SSE = zeros(n_rows1, 1, 'double');
    for j = 1:n_rows1
%         for j = 395
        sprintf('%d, %d', j, i)

        px = squeeze(img_temp2(j, i, :, :));
        [fitresult, gof] = pixelFourierS(pol_angle1, px, optical_period);  
                        
        %prevent background zeroing
        temp_coeff = coeffvalues(fitresult);
        temp_cond = (sum(temp_coeff(1:5) == 0) == 5);
        if ~temp_cond            
            s = sinDescriptor(fitresult, optical_period);    
            
    %         pixelModels{j, i} = fitresult;
            col_avg(j) = s.avg_f;
            col_range(j) = s.range_f;
            col_phase(j) = s.avg_x;
            col_max_x(j) = s.max_x;
            col_SSE(j) = gof.sse;

%                 s.approx %debug
        else            
            col_avg(j) = -1;
            col_range(j) = -1;
            col_phase(j) = -1;
            col_max_x(j) = -1;
            col_SSE(j) = -1;
        end
    end
    pixel_avg(:, i) = col_avg;
    pixel_range(:, i) = col_range;
    pixel_phase(:, i) = col_phase;
    pixel_max_x(:, i) = col_max_x;
    pixel_SSE(:, i) = col_SSE;

    %% Debug section
%         figure
%         plot(pol_angle1, px, 'Marker', '+', 'LineStyle','none')
%         hold on
%         plot(fitresult)
%         hold off
%         grid on
% 
%         sinDescriptorPlot(fitresult, s, optical_period)
end

end