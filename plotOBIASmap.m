function plotOBIASmap(img_gray, S, plotSetup, annotationNames, saveDir)
%OBIAS map with all phases

%Plot settings
contentType = 'image'; %vector
range_ID = plotSetup.range_ID;
cmap = plotSetup.cmap;
transparentVal = plotSetup.transparentVal;
edgeAlpha = plotSetup.edgeAlpha;
lineWidth = plotSetup.lineWidth; %vector thickness


colorRange = makesymbolspec("Polygon", ...
    {'ClassID', [min(range_ID), max(range_ID)], 'FaceColor', cmap, ...
    'FaceAlpha', transparentVal, 'EdgeColor', [0 0 0], ...
    'EdgeAlpha', edgeAlpha, 'LineWidth', lineWidth});

hFig = figure;

imshow(img_gray)
hold on
mapshow(S, 'SymbolSpec', colorRange) %30 sec plot
hold off

%colorbar
colormap(cmap)
n_masks1 = size(cmap, 1);
posTicks = (1:(n_masks1 + 1)) - 0.5;
Ticks = posTicks/n_masks1; %scaled from 0 to 1

c = colorbar('eastoutside', 'Ticks', Ticks, 'TickLabels', annotationNames, ...
    'TickDirection', 'out', 'TickLength', 0.005); %include empty pixels
set(c, 'YDir', 'reverse');
c.AxisLocation = 'out';
c.FontSize = 8;
c.Label.String = 'Mineral';
c.Label.FontSize = 10;

%accommodate colorbar
decrease_by = 0.05; 
axpos = get(gca, 'position');
axpos(3) = axpos(3) - decrease_by;
set(gca, 'position', axpos);

%save figure window
imageFile = strcat('obiasOverlay.pdf');
fullDest = fullfile(saveDir, imageFile); 
hFig.WindowState = 'maximized'; 
ax = gca;
exportgraphics(ax, fullDest, 'ContentType', contentType, 'Resolution', 600)

% figure_frame = getframe(gcf);
% fulldest = fullfile(destFolder_rt, imageFile); 
% imwrite(figure_frame.cdata, fulldest, 'Compression', 'none');
end