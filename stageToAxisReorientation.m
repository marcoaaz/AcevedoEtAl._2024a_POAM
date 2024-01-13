function [pixel_cAxis] = stageToAxisReorientation(pixel_maxPeak)

%Function to convert the stage orientation value at the maximum intensity
%peak and convert it to slow-axis azimuth based on the largest peak found
%for a quartz of known orientation.

%The output depends on lambda-plate transmission axes orientation. We
%followed the convention in: Nikon Eclipse LV1000ND POL has TL-XPL-lambda
%peaks of bluish (135) and yellow (45) colors a definite rotation angles.

%medicine for polarizer orientation; for >=180 values (no longer required)
optical_period = 180; %only applies to TL-XPL-lambda
cond1 = (pixel_maxPeak >= optical_period);
pixel_maxPeak(cond1) = pixel_maxPeak(cond1) - optical_period; %prevent clipping

%Finding c-axis 2D orientation (assuming Pol-EW and Ana-NS)
%convetion: counter clock-wise rotation
cond2 = (pixel_maxPeak >= 0) & (pixel_maxPeak <= 135);
cond3 = (pixel_maxPeak > 135) & (pixel_maxPeak <= 180); %last equal not required
pixel_maxPeak(cond2) = pixel_maxPeak(cond2) + 45;
pixel_maxPeak(cond3) = pixel_maxPeak(cond3) - 135;

pixel_cAxis = pixel_maxPeak;

end