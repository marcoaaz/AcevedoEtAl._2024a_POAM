function [...
    cAxis_peak1, cAxis_peak2, ...
    pixel_peak1, pixel_peak2, ...
    range_peak1, range_peak2] = rearrangePeakImages(...
    img_closest1, img_closest2, ...
    pixel_cAxis1, pixel_cAxis2, ...
    pixel_maxPeakI_discrete1, pixel_maxPeakI_discrete2, ...
    pixel_range1, pixel_range2)

[n_rows_ds, n_cols_ds, n_channels] = size(pixel_maxPeakI_discrete1);

%Finding true images
idx1 = (img_closest1 > img_closest2);
idx2 = ~idx1;
mask1 = cat(3, idx1, idx1, idx1);
mask2 = cat(3, idx2, idx2, idx2);

%preallocate
cAxis_peak1 = zeros(n_rows_ds, n_cols_ds, 1);
cAxis_peak2 = zeros(n_rows_ds, n_cols_ds, 1);
pixel_peak1 = zeros(n_rows_ds, n_cols_ds, n_channels);
pixel_peak2 = zeros(n_rows_ds, n_cols_ds, n_channels);
range_peak1 = zeros(n_rows_ds, n_cols_ds, 1);
range_peak2 = zeros(n_rows_ds, n_cols_ds, 1);

%cAxis orientation
cAxis_peak1(idx1) = pixel_cAxis1(idx1);
cAxis_peak1(idx2) = pixel_cAxis2(idx2);
cAxis_peak2(~idx1) = pixel_cAxis1(~idx1);
cAxis_peak2(~idx2) = pixel_cAxis2(~idx2);

%Peak intensity (8-bit)
pixel_peak1(mask1) = pixel_maxPeakI_discrete1(mask1);
pixel_peak1(mask2) = pixel_maxPeakI_discrete2(mask2);
pixel_peak2(~mask1) = pixel_maxPeakI_discrete1(~mask1);
pixel_peak2(~mask2) = pixel_maxPeakI_discrete2(~mask2);

%Pixel range (continuous value)
range_peak1(idx1) = pixel_range1(idx1);
range_peak1(idx2) = pixel_range2(idx2);
range_peak2(~idx1) = pixel_range1(~idx1);
range_peak2(~idx2) = pixel_range2(~idx2);

%% Optional: Plot ranges
% 
% hFig = figure;
% hFig.Position = [50, 50, 700, 700];
% 
% his1 = histogram(range_peak2, 'FaceColor', [0, 1, 0], 'EdgeColor', 'none');
% hold on
% his2 = histogram(range_peak1, 'FaceColor', [1, 0, 0], 'EdgeColor', 'none');
% hold off
% xlim([0, 255])
% ylim([0, 50000])
% grid on
% 
% legend([his1, his2], {'Peak2', 'Peak1'})
% title('Greyscale pixel ranges (highest order= peak1)')
% xlabel('Channel 8-bit')
% ylabel('Counts')


end