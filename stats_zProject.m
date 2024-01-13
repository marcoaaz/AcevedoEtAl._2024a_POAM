function [time_elapsed] = stats_zProject(fileDir, info_struct, stats_mode, destFolder)
%Ray tracing on individual tiles (cached in RAM)
%Marco A., 20-Oct-22
tic;

pixel_type = 'uint8';
im_height = info_struct.Height; 
im_width = info_struct.Width;
im_channels = info_struct.Channels;
sel_range = info_struct.sel_range;
sel_modality = info_struct.sel_modality;

n_modalities = length(sel_modality);

M = 4; %number of logical processors
%Note: monitor that there is not RAM memory overload

parfor (i = 1:n_modalities, M)
    temp_mtx = [];
    temp_mtx_grey = [];
    img_stats = []; %clearing (runtime error)
    img_stats_index = [];     
    
    tileName = strcat(sel_modality{i}, '.tif');

    %preallocating      
    n_layers = length(sel_range{i});
    k = 0;
    switch stats_mode
        case {'mean', 'max', 'min', 'range', 'sum', 'std', 'median', 'pca'}
            temp_mtx = zeros(im_height, im_width, im_channels, n_layers, pixel_type);        
            for j = sel_range{i}        
                k = k + 1;
                img_temp = imread(fileDir, j); %from BioFormat exporter (>4 GB)
                temp_mtx(:, :, :, k) = img_temp;        
            end  

        case {'maxHSV', 'minHSV', 'rangeHSV'}
            temp_mtx_grey = zeros(im_height, im_width, n_layers, pixel_type); %only for indexing purposes        
            temp_mtx = zeros(im_height, im_width, im_channels, n_layers, pixel_type);        
            for j = sel_range{i}       
                k = k + 1;
                img_temp = imread(fileDir, j); %from BioFormat exporter (>4 GB)
                temp_mtx(:, :, :, k) = img_temp;     
                
                img_V = rgb2hsv(double(img_temp)/255);
                temp_mtx_grey(:, :, k) = img_V(:, :, 3); %[0-1]                   
            end  
    end
    
    %Maths following Fiji>Stack>Z-project
    switch stats_mode
        case 'mean' %double
            img_stats = mean(temp_mtx, 4);%requires uint()    
            
        case 'max' %uint8
            [max_temp, max_index] = max(temp_mtx, [], 4); 
            img_stats = max_temp;
            img_stats_index = max_index;
           
        case 'min' %uint8
            [min_temp, min_index] = min(temp_mtx, [], 4);
            img_stats = min_temp;
            img_stats_index = min_index;

        case 'range' %uint8
            [max_temp, max_index] = max(temp_mtx, [], 4);
            [min_temp, min_index] = min(temp_mtx, [], 4);
            img_stats = max_temp - min_temp;
            img_stats(img_stats < 0) = 0;

        case 'sum' %double
            img_stats = sum(temp_mtx, 4); 

        case 'std' %double
            img_stats = std(single(temp_mtx), 0, 4); %RAM issue:to single, double

        case 'median' %uint8
            img_stats = median(temp_mtx, 4);        
        
        case 'maxHSV' %uint8
            [~, max_index] = max(temp_mtx_grey, [], 3); 
            vec1 = reshape(temp_mtx, [], n_layers);
            vec2 = [max_index(:); max_index(:); max_index(:)];
            vec3 = zeros(size(vec1, 1), 1, 'uint8');
            for m = 1:n_layers
                temp_idx = (vec2 == m);
                vec3(temp_idx) = vec1(temp_idx, m);
            end
            img_stats = reshape(vec3, im_height, im_width, im_channels);
            img_stats_index = max_index;

        case 'minHSV' %uint8
            [~, min_index] = min(temp_mtx_grey, [], 3); 
            vec1 = reshape(temp_mtx, [], n_layers);
            vec2 = [min_index(:); min_index(:); min_index(:)]; %n_channels
            vec3 = zeros(size(vec1, 1), 1, 'uint8');
            for m = 1:n_layers
                temp_idx = (vec2 == m);
                vec3(temp_idx) = vec1(temp_idx, m);
            end
            img_stats = reshape(vec3, im_height, im_width, im_channels)
            img_stats_index = min_index;

        case 'rangeHSV' %uint8
            [~, max_index] = max(temp_mtx_grey, [], 3);
            [~, min_index] = min(temp_mtx_grey, [], 3); 
            vec1 = reshape(temp_mtx, [], n_layers);
            vec2 = [max_index(:); max_index(:); max_index(:)]; 
            vec3 = [min_index(:); min_index(:); min_index(:)]; 
            vec4 = zeros(size(vec1, 1), 1, 'uint8');
            vec5 = zeros(size(vec1, 1), 1, 'uint8')
            for m = 1:n_layers
                temp_idx_max = (vec2 == m); %max
                temp_idx_min = (vec3 == m); %min
                vec4(temp_idx_max) = vec1(temp_idx_max, m);
                vec5(temp_idx_min) = vec1(temp_idx_min, m);
            end
            img_max = reshape(vec4, im_height, im_width, im_channels);
            img_min = reshape(vec5, im_height, im_width, im_channels);
            img_stats = img_max - img_min; 
            %next work: diff. should be in the HSV space
            img_stats(img_stats < 0) = 0;
    end
    check = size(img_stats_index, 1) > 0; %if matrix

    %Saving  
    destTileName = strcat(stats_mode, '_', tileName);
    fullFileName = fullfile(destFolder, destTileName);
    %for max/min linear index
    destTileName_index = strcat(stats_mode, 'Index_', tileName);
    fullFileName_index = fullfile(destFolder, destTileName_index);

    if sum(ismember({'sum', 'std', 'pca'}, stats_mode))

        img_stats2 = single(img_stats); %changing format
        
        %Configure file saving
        t = Tiff(fullFileName, 'w');
        tagstruct = [];
        tagstruct.Photometric = Tiff.Photometric.RGB;
        tagstruct.BitsPerSample = 32;
        tagstruct.SamplesPerPixel = 3;
        tagstruct.SampleFormat = 3;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Software = 'MATLAB';        
        tagstruct.ImageLength = im_height; %constant dim
        tagstruct.ImageWidth = im_width; 
            
        setTag(t, tagstruct)
        write(t, img_stats2);
        close(t);    
    else
        img_stats = uint8(img_stats);        
        imwrite(img_stats, fullFileName, 'compression', 'none')

        %exclusively for index in max/min
        if check %in transparent workspace
            %For tile
%             img_stats_index = uint8(img_stats_index); %<256 steps
            %For montages
            img_stats_index = uint8(rescale(img_stats_index, ...
                0, 255, 'InputMin', 1, 'InputMax', n_layers)); %<256 steps
            imwrite(img_stats_index, fullFileName_index, 'compression', 'none')
        end
       
    end
    
%     clear img_stats img_stats2
    fprintf('tile #: %d, stat: %s \n', i, stats_mode)
end
fprintf('%s completed.\n', stats_mode)

time_elapsed = toc;
end