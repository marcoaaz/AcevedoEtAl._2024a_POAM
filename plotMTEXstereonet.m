function plotMTEXstereonet(destFile_ori_target)

[saveDir, csvName] = fileparts(destFile_ori_target);

% axes convention
setMTEXpref('xAxisDirection', 'south'); %setMTEXpref('yAxisDirection','north');
setMTEXpref('zAxisDirection', 'intoPlane');

cs = crystalSymmetry.load('quartz.cif');

h = Miller({0,0,1}, cs); %interrogated planes {1,0,0}, {0,1,0}, {0,0,1}
ori = orientation.load(destFile_ori_target, cs, 'ColumnNames', {'phi1', 'Phi', 'phi2'}); %odf
odf = calcDensity(ori, 'halfwidth', 2*degree); %orientation density function; 'zero_range', 
% pf = calcPoleFigure(odf, h, equispacedS2Grid('antipodal')); %simulate diffraction counts

% Optical data plots

hFig = figure;

%odf or ori works; 'contour'; 'antipodal'; 'logarithmic'
plotPDF(odf, h, ...
    'points', 'all', 'resolution', 1.5*degree, 'linewidth', 1, ...
    'grid', 'grid_res', 15*degree, 'DisplayName', 'C-axis') 
annotate(zvector, 'label', {'Z'}, 'BackgroundColor', 'w')
mtexColorbar
CLim(gcm, 'equal');
% mtexColorMap
plot(ori, 'MarkerSize', 1, 'MarkerColor', 'magenta', 'add2all')

hFig.WindowState = 'maximized'; 
ax = gca;
imageFile = strcat(csvName, '.pdf');
fullDest = fullfile(saveDir, imageFile); 
exportgraphics(ax, fullDest, 'ContentType', 'vector') %, 'Resolution', 600

end