function [img_sum_w_rs] = stack_imageFusion(type, mode, weights, destFile, saveOption)

workingDir = fileparts(destFile);

n_type = length(type);
n_mode = length(mode);
weights = 100*weights/sum(weights);

fileName_info = fullfile(workingDir, strcat(type{1}, '_', mode{1}, '.tif'));
info_st = imfinfo(fileName_info);
n_rows = info_st.Height;
n_cols = info_st.Width;
n_channels = info_st.SamplesPerPixel;

RT_stack = zeros(n_rows, n_cols, n_channels, n_mode, 'double');
for i = 1:n_type
    for j = 1:n_mode
        fileName = fullfile(workingDir, strcat(type{i}, '_', mode{j}, '.tif'));        
        RT_stack(:, :, :, j) = imread(fileName);
    end
end

img_sum_w = zeros(n_rows, n_cols, n_channels, 'double');
for m = 1:n_mode
    img_sum_w = img_sum_w + RT_stack(:, :, :, m)*weights(m);
end
img_sum_w_rs = uint8(rescale(img_sum_w, 0, 255));

if saveOption == 1
    imwrite(img_sum_w_rs, destFile);
end

%% Further experiments

% %test
% img_sum = sum(RT_stack, 4);
% img_sum_rs = uint8(rescale(img_sum, 0, 255));

% figure
% imshow(img_sum_rs)
% 
% figure
% imshow(img_sum_w_rs)

end


