%https://au.mathworks.com/matlabcentral/answers/448477-how-to-draw-a-l-a-b-colour-space-or-hsv-colour-space-2d

plotradius = 1000;
plotsaturation = 1;
[x, y] = meshgrid(-plotradius:plotradius);  %get x,y coordinates of pixels in image

%convert to angle, radius. Angle is the hue, radius is the saturation in the HSV cylinder
[hue, value] = cart2pol(x, y);  %clock-wise from X-axis

%rescale from -pi:pi to 0:1 since matlab use 0:1 for the range of the hue
hue2 = rescale(hue, 0, 1, 'InputMin', 0, 'InputMax', pi);

%rescale saturation to range 0:1. Not that anything above 1 is outside the cylinder
value = value/plotradius; 

%set value constant for all points in the disk
saturation = ones(size(hue)) * plotsaturation;  

%now set points outside the disk (saturation >1) to white. That'd be a value of 1 and saturation of 0. hue doesn't matter
outsidedisk = (value > 1) | (hue < 0);
value(outsidedisk) = 1;
saturation(outsidedisk) = 0;

%finally, convert hsv to rgb and plot
rgb = hsv2rgb(cat(3, hue2, saturation, value));

from_row = floor(plotradius) + 1;

figure
imshow(rgb(from_row:end, :, :))