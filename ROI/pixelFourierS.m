function [fitresult, gof, preparedData] = pixelFourierS(pol_angle, px, optical_period)
%Fits a Fourier model to a pixel
%Zeros are included in the frames due to rotation. They affect Fourier-2
%fitting when rotational intervals are irregularly skipped. Lesser points
%with intervals produce random fits that can give negative values. negative
%values. Additionally, the resize interpolation incorporates pixels that
%were zero background

%Input:
%pol_angle = nominal stage angle series (manual experiment)
%px_rgb = pixel in double values
%color_space = 1 %1=RGB, 2=CieLAB, 3=HSV
%optical_period = 90 for PPL; 180 for XPL and XPL-lambda plate

gof = struct; 
gof.sse = 0; %Unrecognized field name "SSE".

%%
quality_ratio = 0.8; %reject misfitting 

n_point = length(px);
n_point_min = 6;
w_fixed = 2*pi/(optical_period); %default period

%skip fringe pixels [0.1 garnet TL-XPL, 3 most minerals] (px == 0 background)
idx_in = ~(px < 0.1); 
n_in = sum(idx_in);

pol_angle = pol_angle(idx_in); %check for transposing
px = px(idx_in);    
[xData_temp, yData_temp] = prepareCurveData(pol_angle, px);     
preparedData = [xData_temp, yData_temp];

if (n_in >= n_point*quality_ratio) && (n_point > n_point_min)            
    %Fourier-2
    [fitresult, gof] = createFit_two(xData_temp, yData_temp, w_fixed);   
else
    f = fittype('a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w)');
    fitresult = cfit(f, 0, 0, 0, 0, 0, w_fixed);       
end

end