function [img_closest1, img_closest2] = highestOrderPeak(...
    pixel_maxPeakI_discrete1, pixel_maxPeakI_discrete2, ...
    EQ_struct, strip_michelL, plotOption)

% waveplateRetardation = 500;
% maxRetardation = 2500;
% strip_threshold = floor(waveplateRetardation*n_levels/maxRetardation);

[n_rows, n_cols, n_channels] = size(pixel_maxPeakI_discrete1);
n_levels = size(strip_michelL, 1);

%CieLAB
pixel_maxPeakI_discrete1 = rgb2lab(pixel_maxPeakI_discrete1./255);
pixel_maxPeakI_discrete2 = rgb2lab(pixel_maxPeakI_discrete2./255);

%plot settings
exponentValue = 2;
aspectRatio = 3;
chartSize = 400;
lineWidth = 2;

%interference order approximation
img_closest1 = zeros(n_rows, n_cols, 'uint16');
img_closest2 = zeros(n_rows, n_cols, 'uint16');
for i = 1:n_rows %n_rows
% for i =  2328%n_cols
    sprintf('row: %d', i)    
    row_closest1 = zeros(1, n_cols, 'uint16'); 
    row_closest2 = zeros(1, n_cols, 'uint16'); 

    parfor j = 1:n_cols %n_cols
    % for j = 2502 %n_cols %for Plot debugging
        
        %sRGB
        px1 = squeeze(pixel_maxPeakI_discrete1(i, j, :));
        px2 = squeeze(pixel_maxPeakI_discrete2(i, j, :));
                
        px1_adapted = [
            feval(EQ_struct(1).function, px1(1)), ...
            feval(EQ_struct(2).function, px1(2)), ...
            feval(EQ_struct(3).function, px1(3))];
        px2_adapted = [
            feval(EQ_struct(1).function, px2(1)), ...
            feval(EQ_struct(2).function, px2(2)), ...
            feval(EQ_struct(3).function, px2(3))];               
        
        %Finding closest colour
        D1 = sqrt(sum( abs(strip_michelL - repmat(px1_adapted, [n_levels, 1])).^exponentValue, 2 ));
        D2 = sqrt(sum( abs(strip_michelL - repmat(px2_adapted, [n_levels, 1])).^exponentValue, 2 ));
        [D1_min, closestIndex1] = min(D1, [], 'all');   
        [D2_min, closestIndex2] = min(D2, [], 'all');   
        
        row_closest1(j) = closestIndex1;
        row_closest2(j) = closestIndex2;
    end
    img_closest1(i, :) = row_closest1;
    img_closest2(i, :) = row_closest2;

     %% Plot debugging (1 pixel)
        if plotOption == 1
                        
            %Info
            px1
            px2
            px1_adapted
            px2_adapted
            
            %Quality control
            square_size = 200;
            img_ones = ones(square_size, square_size, 3);
            img_ones(:, :, 1) = img_ones(:, :, 1).*(px1(1));
            img_ones(:, :, 2) = img_ones(:, :, 2).*(px1(2));
            img_ones(:, :, 3) = img_ones(:, :, 3).*(px1(3));
            img_ones1 = img_ones;
            img_ones = ones(square_size, square_size, 3);
            img_ones(:, :, 1) = img_ones(:, :, 1).*(px2(1));
            img_ones(:, :, 2) = img_ones(:, :, 2).*(px2(2));
            img_ones(:, :, 3) = img_ones(:, :, 3).*(px2(3));
            img_ones2 = img_ones;
            img_ones = ones(square_size, square_size, 3);
            img_ones(:, :, 1) = img_ones(:, :, 1).*(px1_adapted(1));
            img_ones(:, :, 2) = img_ones(:, :, 2).*(px1_adapted(2));
            img_ones(:, :, 3) = img_ones(:, :, 3).*(px1_adapted(3));
            img_ones3 = img_ones;
            img_ones = ones(square_size, square_size, 3);
            img_ones(:, :, 1) = img_ones(:, :, 1).*(px2_adapted(1));
            img_ones(:, :, 2) = img_ones(:, :, 2).*(px2_adapted(2));
            img_ones(:, :, 3) = img_ones(:, :, 3).*(px2_adapted(3));
            img_ones4 = img_ones;
            
            img_ones_montage = lab2rgb([img_ones1, img_ones2; img_ones3, img_ones4]);
    
            x_str= [square_size/2, 3*square_size/2, square_size/2, 3*square_size/2];
            y_str= [square_size/2, square_size/2, 3*square_size/2, 3*square_size/2];
            
            %Intensity difference plot
            figure,
            imshow(img_ones_montage)
            hold on
            text(x_str, y_str, {'px1', 'px2', 'px1 adapted', 'px2 adapted'}, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle')    
            
            hFig = figure;
            hFig.Position = [50, 50, chartSize*aspectRatio, chartSize];
            ax = gca;
                        
            plot(1:n_levels, D1, 'Color', [1, 0, 1], 'LineWidth', lineWidth, 'DisplayName', 'peak1')
            hold on
            plot(1:n_levels, D2, 'Color', [0, 1, 0], 'LineWidth', lineWidth, 'DisplayName', 'peak2')
            
            %min lines
            plot(xlim(ax), [D1_min, D1_min], 'LineWidth', lineWidth*0.8, ...
                'Color', [1, 0, 1, 0.5], 'LineStyle','--');
            plot(xlim(ax), [D2_min, D2_min], 'LineWidth', lineWidth*0.8, ...
                'Color', [0, 1, 0, 0.5], 'LineStyle','--');
            plot([closestIndex1, closestIndex1], ylim(ax), 'LineWidth', lineWidth*0.8, ...
                'Color', [1, 0, 1, 0.5], 'LineStyle','--');
            plot([closestIndex2, closestIndex2], ylim(ax), 'LineWidth', lineWidth*0.8, ...
                'Color', [0, 1, 0, 0.5], 'LineStyle','--');
            hold off

            %customization
            grid(ax, 'minor')
            xlim([0, 1000])
            % ylim([0, 10000])
            legend('Location','eastoutside')            
        end

end


end