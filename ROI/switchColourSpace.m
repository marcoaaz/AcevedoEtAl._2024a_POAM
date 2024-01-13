function [px, channelNames] = switchColourSpace(px_rgb, color_space)
%color_space= 1 %1=RGB, 2=CieLAB, 3=HSV

switch color_space
    case 1

        %Optional: Linearize RGB
        % px_rgb = rgb2lin(px_rgb, ...
        %     'OutputType', 'double', 'ColorSpace','srgb'); %linearize (no longer 8-bit)
        %'double' | 'single' | 'uint8' | 'uint16'
        
        px = px_rgb;
        channelNames = {'R', 'G', 'B'};
    case 2
        px_lab = rgb2lab(px_rgb', 'ColorSpace', 'adobe-rgb-1998'); 
        px = px_lab';
        channelNames = {'L', 'a', 'b'};
    case 3
        px_hsv = rgb2hsv(px_rgb'/255);
        px = px_hsv';
        channelNames = {'H', 'S', 'V'};
end

% px = double(px);

end