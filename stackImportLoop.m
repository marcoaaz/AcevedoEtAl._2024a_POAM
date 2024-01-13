function [stack_output, quality_mask2] = stackImportLoop(imageName, info_struct, varRange, importType, scale)

n_rows = info_struct.Height; 
n_cols = info_struct.Width;
n_channels = info_struct.Channels;
n_layers = length(varRange);      

quality_mask = zeros(n_rows, n_cols, n_layers, "double");

switch importType

    case 1 %greyscale

        stack_temp = zeros(n_rows, n_cols, n_layers, "double"); %3GB if 1651x1651                
        k = 0;
        for i = varRange
            k = k + 1;
            temp_img = imread(imageName, i);
            
            %linearize uint8
            temp_img1 = rgb2lin(double(temp_img)/255, 'ColorSpace', 'srgb');
        
            %Greyscale            
            channel_temp = 255*mean(temp_img1, 3); %greyscale
            % channel_temp = 255*temp_img(:, :, sel_channel(uu)); %RGB channel
            
            stack_temp(:, :, k) = channel_temp;            
        
            %Data quality mask
            quality_mask(:, :, k) = (channel_temp ~= 0);
        end
        
        
    case 2 %rgb

        stack_temp = zeros(n_rows, n_cols, n_channels, n_layers, 'double');        
        k = 0;
        for i = varRange
            k = k + 1;
            temp_img = imread(imageName, i);
            
            %linearize uint8
            temp_img1 = rgb2lin(double(temp_img)/255, 'ColorSpace', 'srgb');
        
            %Greyscale            
            channel_temp = 255*mean(temp_img1, 3); %greyscale
            % channel_temp = 255*temp_img(:, :, sel_channel(uu)); %RGB channel
                        
            stack_temp(:, :, :, k) = 255*temp_img1;
        
            %Data quality mask
            quality_mask(:, :, k) = (channel_temp ~= 0);
        end      

end
quality_mask1 = sum(quality_mask, 3)/k; %fraction

%Downscaling (for performance gain)
stack_output = imresize(stack_temp, scale, "bilinear");
quality_mask2 = imresize(quality_mask1, scale, "bilinear");

end