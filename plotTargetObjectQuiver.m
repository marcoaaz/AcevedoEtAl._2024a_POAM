function plotTargetObjectQuiver(...
    img_fused, S_sub_filtered, quiver_vectors_obj,  plotSetup, mineralText, saveDir)
% Object quiver plot overlaid on gridded map

%Plot settings
contentType = 'image'; %vector
range_ID = plotSetup.range_ID;
cmap = plotSetup.cmap;
transparentVal = plotSetup.transparentVal;
edgeAlpha = plotSetup.edgeAlpha;
lineWidth = plotSetup.lineWidth; %vector thickness
autoSFactor = plotSetup.autoSFactor;

colorRange = makesymbolspec("Polygon", ...
    {'ClassID', [min(range_ID), max(range_ID)], 'FaceColor', cmap, ...
    'FaceAlpha', transparentVal, 'EdgeColor', [1 1 1], ...
    'EdgeAlpha', edgeAlpha, 'LineWidth', lineWidth});

%Orientation color bar
min_bar = 0; 
max_bar = 180;
interval_val = 20; 
thick_vals_orig = min_bar:interval_val:max_bar;
thick_labels = strsplit(num2str(thick_vals_orig));
thick_vals = rescale(thick_vals_orig, 0, 1); %scaled to triplet

hFig = figure; %quiver plot
ax = gca;

imshow(img_fused) %option 2
hold on
mapshow(S_sub_filtered, 'SymbolSpec', colorRange) %30 sec plot
quiver(quiver_vectors_obj(1, :), quiver_vectors_obj(2, :), ...
    quiver_vectors_obj(3, :), quiver_vectors_obj(4, :), ...
    'AutoScaleFactor', autoSFactor, 'Color', 'white', 'LineWidth', 2*lineWidth)
axis equal
hold off
colormap(hsv(255))
cb = colorbar(ax, 'Location', 'eastoutside', 'Ticks', thick_vals, ...
    'TickLabels', thick_labels);
title(ax, mineralText)

%save figure window
imageFile = strcat('quiverOverlay_', mineralText, '.pdf');
fullDest = fullfile(saveDir, imageFile); 
hFig.WindowState = 'maximized'; 
exportgraphics(ax, fullDest, 'ContentType', contentType, 'Resolution', 600)

% figure_frame = getframe(gcf);
% fulldest = fullfile(destFolder_rt, imageFile); 
% imwrite(figure_frame.cdata, fulldest, 'Compression', 'none');
end