%% Script to validate slow-axis Azimuth: biaxial slow-axis azimuth with EBSD data
%Sample: harzburgite basement mylonite 17-BSK-035
%Target: orthopyroxene (enstatite)
%Date: 3-01-2024
%Last update: Marco Acevedo Z.

%refined EBSD map
ebsd_ori = readtable('E:\Alienware_March 22\02-Research geology\05-QUT\Paper 2\EBSD_harzburgite 17-BSK-035_opx\ebsd_objects.csv');

%POAM orientation quiver map centroids (after OBIAS)
poamFile = 'E:\paper 2_datasets\nikon LV1000ND\17BSK035\stacks_save\composite_trial4\all_modalities_rayTracing\obiasExport_0.5_review2\wave2_TL XPL-lambda_0.5_objNet_orthopyroxene_centroidAndVector2.csv';

%5X
matched_EBSD = [7505, 7428, 6128, 5730, 5320, 4671, 4120, 10655, ...
    10029, 8414, 7798, 7025, 4276, 4838];

matched_coordinates = [ %POAM
    2952.09	2878.57;
    3165.73	2854.92;
    2887.08	3039.74;
    3378.81	3123.31;
    3170.98	3163.51;
    3375.71	3321.02;
    3189.32	3354.66;
    3600.34	2265.34;
    3581.26	2365.28;     
    3410.88	2720.2;
    3743.77	2776.56;
    3608.63	2965.8;
    3612.55	3364.36;
    3778.56	3287.55
    ];
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

%Finding optical objects
quiver_vectors_obj = readmatrix(poamFile);
accuracy_val = 0.005;

optical_ori2 = [];
for k = 1:n_grains
    idx = (abs(matched_coordinates(k, 1) - quiver_vectors_obj(1, :)') < accuracy_val) & ...
        (abs(matched_coordinates(k, 2) - quiver_vectors_obj(2, :)') < accuracy_val);

    row_temp = quiver_vectors_obj(:, idx)';
    optical_ori2 = [optical_ori2; row_temp];
end

%Calculating azimuth
optical_u = optical_ori2(:, 3);
optical_v = optical_ori2(:, 4);
optical_in = (180/pi)*acos(sqrt(optical_u.^2 + optical_v.^2));
optical_az = (180/pi)*atan2(optical_v, optical_u);

ebsd_u = ebsd_ori2{:, 'x'}; %table
ebsd_v = ebsd_ori2{:, 'y'};
ebsd_az = 180 + (180/pi)*atan2(ebsd_v, ebsd_u);

%%

%Plot setup
axisLimit_bottom = 1;
axisLimit_top = 180;
nrows = 1;
ncolumns = 2;
textAlignment = 'left';
fontSize = 12;
tl_pct = 0.03;
colorSel = {[1, 0, 0, 0.5], [0, 0, 1, 0.5]};
greySel = 0.8;
d_offset = 0.01; %text offset (*axisLimit_top)
tl1_x = axisLimit_top*tl_pct; %text top left
tl2_x = tl1_x;
tl1_y = axisLimit_top*(1-5*tl_pct);
tl2_y = axisLimit_top*(1-8*tl_pct);

ft1 = fittype({'x'}); %forzed zero intercept
ft2 = fittype({'x', '1'});

%azimuth
x= optical_az'; 
y= ebsd_az';
b_fit = linspace(axisLimit_bottom, axisLimit_top, 150); 
%default = 10; more steps required if log-scale axes

[p1_b, p1_gof_b] = fit(x', y', ft1);
[p2_b, p2_gof_b] = fit(x', y', ft2);
y1_fitted_b = feval(p1_b, b_fit);
y2_fitted_b = feval(p2_b, b_fit);

%str
coeffval_temp = coeffvalues(p1_b);
b_str_eq1 = sprintf('y = %0.2f*x ; rsquare = %0.2f', ...
    coeffval_temp(1), p1_gof_b.rsquare);

coeffval_temp = coeffvalues(p2_b);
b_str_eq2 = sprintf('y = %0.2f*x + %0.2f ; rsquare = %0.2f', ...
    coeffval_temp(1), coeffval_temp(2), p2_gof_b.rsquare);


hFig = figure;
hFig.Position = [100, 100, 1200, 600];
clf('reset')

ax = gca;

%fits
plot(b_fit, y1_fitted_b, '-', 'LineWidth', 2, 'Color', colorSel{1})
hold on
plot(b_fit, y2_fitted_b, '-', 'LineWidth', 2, 'Color', colorSel{2})
%1:1 line
% plot([0, axisLimit2], [0, axisLimit2], '-', 'LineWidth', 1, 'Color', [0, 0, 0, 0.5])
%points
plot(x, y, '.black', 'MarkerSize', 15)
%labels
text(x + d_offset, y - d_offset, num2str(ebsdID2), ...
    'HorizontalAlignment', textAlignment, 'Color', 'black', 'FontSize', fontSize)
%fitted eq
text(tl1_x, tl1_y, b_str_eq1, ...
    'HorizontalAlignment', textAlignment, 'Color', colorSel{1}(1:end-1), ...
    'FontSize', fontSize, 'BackgroundColor', [greySel, greySel, greySel])
text(tl2_x, tl2_y, b_str_eq2, ...
    'HorizontalAlignment', textAlignment, 'Color', colorSel{2}(1:end-1), ...
    'FontSize', fontSize, 'BackgroundColor', [greySel, greySel, greySel])
hold off

xlim([1, axisLimit_top])
ylim([1, axisLimit_top])
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
%appereance
xlabel('Optical method', 'FontSize', fontSize)
ylabel('SEM-EBSD mapping', 'FontSize', fontSize)
ax.XAxis.FontSize = fontSize;
ax.YAxis.FontSize = fontSize;
title(['Slow-axis azimuth between [0-180', char(176), ']'], 'Units', 'normalized', ...
    'Position', [0.5, 0.92, 0], 'FontSize', 1.5*fontSize)
grid on
% grid(ax, 'minor')
