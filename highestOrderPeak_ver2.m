function [img_closest1, img_closest2] = highestOrderPeak_ver2(...
    pixel_maxPeakI_discrete1, pixel_maxPeakI_discrete2, ...
    modelFile, chartFile, workingDir, plotOption)

img_michelL = imread(chartFile);
strip_michelL = double(squeeze(img_michelL(1, :, :)))./255; %[0-1]

%sRGB to CieLAB
strip_michelL_lab = rgb2lab(strip_michelL, 'ColorSpace', 'srgb'); %sRGB
pixel_maxPeakI_discrete1_lab = rgb2lab(pixel_maxPeakI_discrete1./255, ...
    "ColorSpace", "linear-rgb"); %linearised
pixel_maxPeakI_discrete2_lab = rgb2lab(pixel_maxPeakI_discrete2./255, ...
    "ColorSpace", "linear-rgb");

%obtaining function handles
k = load(modelFile, "-mat");
EQ_struct = k.EQ_struct;

R_cfit = EQ_struct(1).function;
G_cfit = EQ_struct(2).function;
B_cfit = EQ_struct(3).function;
R_formula_str = strrep(formula(R_cfit), newline, ''); %cfit to sym
G_formula_str = strrep(formula(G_cfit), newline, '');
B_formula_str = strrep(formula(B_cfit), newline, '');
R_f = subs(str2sym(R_formula_str), coeffnames(R_cfit), num2cell(coeffvalues(R_cfit).'));
G_f = subs(str2sym(G_formula_str), coeffnames(G_cfit), num2cell(coeffvalues(G_cfit).'));
B_f = subs(str2sym(B_formula_str), coeffnames(B_cfit), num2cell(coeffvalues(B_cfit).'));
R_F = matlabFunction(R_f); %sym to function handle
G_F = matlabFunction(G_f);
B_F = matlabFunction(B_f);

%colour matching
exponentValue = 2; 
[n_rows, n_cols, n_channels] = size(pixel_maxPeakI_discrete1);
n_levels = size(strip_michelL_lab, 1);

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
        
        px1 = squeeze(pixel_maxPeakI_discrete1_lab(i, j, :));
        px2 = squeeze(pixel_maxPeakI_discrete2_lab(i, j, :));
        
        %colour adaptation in each channel
        px1_adapted = [R_F(px1(1)), G_F(px1(2)), B_F(px1(3))];
        px2_adapted = [R_F(px2(1)), G_F(px2(2)), B_F(px2(3))];           
        
        %Finding closest colour (L2 or Euclidean norm)
        D1 = sqrt(sum( abs(strip_michelL_lab - repmat(px1_adapted, [n_levels, 1])).^exponentValue, 2 ));
        D2 = sqrt(sum( abs(strip_michelL_lab - repmat(px2_adapted, [n_levels, 1])).^exponentValue, 2 ));
        [D1_min, closestIndex1] = min(D1, [], 'all');   
        [D2_min, closestIndex2] = min(D2, [], 'all');   
        
        row_closest1(j) = closestIndex1;
        row_closest2(j) = closestIndex2;

        % %Optional: Info
        % px1
        % px2
        % px1_adapted
        % px2_adapted
    end
    img_closest1(i, :) = row_closest1;
    img_closest2(i, :) = row_closest2;
    %Note: this value depends on the spatial resolution of the Michel-Levy
    %chart simulation ('n_levels') but it can be converted to birefringence
    %It would be an approximation depending on the polynomial fit to the
    %real experiment.

     %% Plot for debugging (use only 1 pixel)

        if plotOption == 1         

            %% Intensity difference plot for quality control

            square_size = 200;
            img_ones1 = ones(square_size, square_size, 3).*reshape(px1, 1, 1, []);
            img_ones2 = ones(square_size, square_size, 3).*reshape(px2, 1, 1, []);
            img_ones3 = ones(square_size, square_size, 3).*reshape(px1_adapted, 1, 1, []);
            img_ones4 = ones(square_size, square_size, 3).*reshape(px2_adapted, 1, 1, []);

            img_ones_montage = lab2rgb([img_ones1, img_ones2; img_ones3, img_ones4]);
    
            x_str= [square_size/2, 3*square_size/2, square_size/2, 3*square_size/2];
            y_str= [square_size/2, square_size/2, 3*square_size/2, 3*square_size/2];            
            
            figure,
            imshow(img_ones_montage)
            hold on
            text(x_str, y_str, {'px1', 'px2', 'px1 adapted', 'px2 adapted'}, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle')    
            
            %% Plot for debugging (requires deactivating the parfor loop)

            %plot settings
            aspectRatio = 3;
            chartSize = 400;
            lineWidth = 2;

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

imwrite(img_closest1, fullfile(workingDir, 'img_closest1.tif'))
imwrite(img_closest2, fullfile(workingDir, 'img_closest2.tif'))

end
