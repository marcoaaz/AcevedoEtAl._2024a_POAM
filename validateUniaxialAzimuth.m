%% Script to validate optic-axis orientation (Azimuth/Dip): uniaxial c-axis with EBSD data
%Sample: granite mylonite 384D22
%Target: quartz
%Date: 3-01-2024
%Last update: Marco Acevedo Z.

%refined EBSD map
ebsd_ori = readtable('E:\paper 2_datasets\nikon LV1000ND\mylonite\zoom-in registeredStacks\ebsdPolar.csv');

%POAM orientation quiver map centroids (after OBIAS)
%10X
% poamFile = 'E:\paper 2_datasets\nikon LV1000ND\mylonite\zoom-in registeredStacks\all_modalities_rayTracing\obiasExport_linear_0.5_review3\30Dec23_modulationImage_TL XPL-lambda_0.5_objNet_quartz_centroidAndVector2.csv';
%5X
poamFile = 'E:\paper 2_datasets\nikon LV1000ND\mylonite\final stacks\all_modalities_rayTracing\obiasExport_0.5_linear\6-Jan-24modulationImage_TL XPL-lambda_0.5_objNet_quartz_centroidAndVector.csv';

% %10X
% matched_EBSD = [407, 230, 74, 210, 234, 446, 268, 369];
% 
% matched_coordinates = [
%     610.623, 984.149;
%     974.507, 1170.25;
%     628.073, 1307.14;
%     803.882, 1300.01;
%     653.647, 1214.08;
%     976.897, 896.574;
%     501.776, 1153.44;
%     833.484, 1057.95;    
% ];

% 5X
matched_EBSD = [402, 407, 369, 446, 230, 210, 234, 325];
matched_coordinates = [
    2139.14	1687.27; 
    2227.7	1729.91;
    2331.38	1761.37;
    2404.69	1676.02;
    2399.89	1816.85;
    2351.04	1852.98;
    2243.14	1842.83;
    2168.26	1752.09];

n_grains = length(matched_EBSD);

%Finding ebsd
ebsdID = ebsd_ori.ID;
ebsd_ori2 = [];
for k = 1:n_grains
    idx = (ebsdID == matched_EBSD(k));
    row_temp = ebsd_ori(idx, :);
    ebsd_ori2 = [ebsd_ori2; row_temp];
end
ebsdID2 = ebsd_ori2.ID;

%Make EBSD comparable to Optical modelling
theta_temp = ebsd_ori2.theta;
rho_temp = ebsd_ori2.rho;
lowerHemi = (theta_temp > pi/2);
theta_temp = abs(pi/2 - theta_temp);
rho_temp(lowerHemi) = rho_temp(lowerHemi) + pi;
rho_temp(rho_temp < 0) = rho_temp(rho_temp < 0) + pi;
rho_temp(rho_temp > pi) = rho_temp(rho_temp > pi) - pi;

%Finding optical objects
quiver_vectors_obj = readmatrix(poamFile);
accuracy_val = 6; %default = 0.005

optical_ori2 = [];
for k = 1:n_grains
    idx = (...
        (abs(matched_coordinates(k, 1) - quiver_vectors_obj(1, :)') < accuracy_val) & ...
        (abs(matched_coordinates(k, 2) - quiver_vectors_obj(2, :)') < accuracy_val)...
        );

    row_temp = quiver_vectors_obj(:, idx)';
    optical_ori2 = [optical_ori2; row_temp];
end
%Note: recalculated modulation images (mask_bg) might have 6 px offset.

%Converting U and V to inclination (a1) and azimuth (b1)
a = optical_ori2(:, 3);
b = optical_ori2(:, 4);
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