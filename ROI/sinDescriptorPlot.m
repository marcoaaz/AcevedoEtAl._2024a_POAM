function sinDescriptorPlot(sf, s, xData, yData)
%sf = cfit of Fourier-2
%s = structure describing sf

%period = modality sexagesimal period (ppl = 180, xpl=90)
w = sf.w;
period = 2*pi/w; %period
x_plot = 1:1:2*period; %depends on data acquisition

%lines
avg_f = s.avg_f;
avg_x = s.avg_x;
A = s.Inflexion;
B2 = s.MinMax;
% range1 = s.range1;
% range2 = s.range2;

%points
zero_list = A(:, 1); %inflexion
f_list = A(:, 2);
zero_list2 = B2(:, 1); %critical
f_list2 = B2(:, 2);

%lines
min_x = B2(end, 1);
min_f = B2(end, 2);
max_x = B2(1, 1);
max_f = B2(1, 2);

%% (Plot) cfit to sym

%Settings
lineWidth = 2;
pointSize = 13;
fontSize = 13;
gridAlpha = 0.7;

hFig = figure(10);
hFig.Position = [100, 100, 1000, 600];
% hFig.Position = [100, 100, 1200, 800];

clf('reset') %clearing the automatic plot

%empirical data
h1 = plot(xData, yData, '.black', 'Marker', 'o', 'Color', 'magenta', ...
    'MarkerSize', 0.5*pointSize, 'LineStyle', 'none', ...
    'MarkerFaceColor', 'black', 'MarkerEdgeColor','none', ...
    'DisplayName', 'Pixel values');
ax= gca;
hold on

%functions
h2 = plot(feval(sf, x_plot), 'LineWidth', 2*lineWidth, ...
    'Color', [0, 0, 0, 0.5], 'LineStyle', '-', 'DisplayName', 'Fourier-2');
% plot(subs(d1_sf, x_plot)) %evaluate sym expression
% plot(subs(d2_sf, x_plot))

%lines
h5 = plot([min_x, min_x], ylim(ax), 'LineWidth', lineWidth, 'Color', 'b', 'LineStyle','-', ...
    'DisplayName', sprintf('min x = %0.1f', min_x));
h6 = plot(xlim(ax), [min_f, min_f], 'LineWidth', lineWidth, 'Color', 'b', 'LineStyle','-', ...
    'DisplayName', sprintf('min y = %0.1f', min_f));
h7 = plot([avg_x, avg_x], ylim(ax), 'LineWidth', lineWidth, 'Color', 'k', 'LineStyle','-', ...
    'DisplayName', sprintf('mean x = %0.1f', avg_x));
h8 = plot(xlim(ax), [avg_f, avg_f], 'LineWidth', lineWidth, 'Color', 'k', 'LineStyle','-', ...
    'DisplayName', sprintf('mean y = %0.1f', avg_f));
h9 = plot([max_x, max_x], ylim(ax), 'LineWidth', lineWidth, 'Color', 'r', 'LineStyle','-', ...
    'DisplayName', sprintf('max x = %0.1f', max_x));
h10 = plot(xlim(ax), [max_f, max_f], 'LineWidth', lineWidth, 'Color', 'r', 'LineStyle','-', ...
    'DisplayName', sprintf('max y = %0.1f', max_f));
h11 = plot([period, period], ylim(ax), 'LineWidth', lineWidth, 'Color', 'black', 'LineStyle','--', ...
    'DisplayName', sprintf('period = %0.1f', period)); %period

%points
%/Inflection
h3 = plot(zero_list, f_list, 'Marker', "pentagram", 'Color', 'cyan', ...
    'MarkerSize', pointSize, 'LineStyle', 'none', ...
    'MarkerFaceColor', 'cyan', 'MarkerEdgeColor','black', ...
    'DisplayName', 'Inflection points');
%/Max (two peaks in TL-XPL-lambda)
h4 = plot(zero_list2, f_list2, 'Marker', 'o', 'Color', 'magenta', ...
    'MarkerSize', pointSize, 'LineStyle', 'none', ...
    'MarkerFaceColor','magenta', 'MarkerEdgeColor', 'black', ...
    'DisplayName', 'Critical points');

xlim([0, 2*period])
hold off
grid on

grid(ax, 'minor')
% ax.GridColor = [0, 0, 0]; %optional
% ax.GridLineWidth = 0.3*lineWidth;
% ax.GridAlpha = gridAlpha;
% ax.MinorGridAlpha = gridAlpha;

ax.XAxis.FontSize = 1*fontSize;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 1*fontSize;
ax.YAxis.FontWeight = 'bold';

xlabel('Angle°', 'Interpreter', 'none', 'FontSize', fontSize, 'FontWeight', 'bold');
ylabel('Channel 8-bit', 'Interpreter', 'none', 'FontSize', fontSize, 'FontWeight', 'bold'); 
% title('Pixel spectral description', 'FontSize', 1.2*fontSize) %optional

legend([h3, h4, h5, h6, h7, h8, h9, h10, h11], ...
    'Location', 'eastoutside', 'FontSize', fontSize)

end