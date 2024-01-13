function [fitresult, gof] = pixelFourierS_ver2(pol_angle, px, w_fixed, f)
%Fits a Fourier model to a pixel
%Zeros are included in the frames due to rotation. They affect Fourier-2
%fitting when rotational intervals are irregularly skipped. Lesser points
%with intervals produce random fits that can give negative values. 
% The resize interpolation might incorporate pixels that
%were zero background around the fringes of the FoV

%Input:
%pol_angle = nominal stage angle series (manual experiment). input as column
%px = pixel value as double

quality_ratio = 0.8; %reject misfitting 
n_point_min = 6;

gof = struct; 
gof.sse = 0; %Unrecognized field name "SSE".

%skip fringe pixels [0.1 garnet TL-XPL, 3 most minerals] (px == 0 background)
n_point = length(px);
idx_in = ~(px < 0.1); 
n_in = sum(idx_in);

if (n_in >= n_point*quality_ratio) && (n_point > n_point_min)            
    
    xData_temp = pol_angle(idx_in); %check for transposing
    yData_temp = px(idx_in);
    yData_temp1 = yData_temp*(2/n_in); %for coefficient calculation    
    
    % [xData_temp, yData_temp] = prepareCurveData(pol_angle, px);     
    % preparedData = [xData_temp, yData_temp];
    % %Fourier-2
    % [fitresult, gof] = createFit_two(xData_temp, yData_temp, w_fixed);     
    
    angle_temp = xData_temp*w_fixed;

    a_0 = mean(yData_temp);
    a_1 = yData_temp1*cos(angle_temp); %x is transposed as column
    b_1 = yData_temp1*sin(angle_temp);
    a_2 = yData_temp1*cos(2*angle_temp); 
    b_2 = yData_temp1*sin(2*angle_temp);
    
    %ordered alphabetically: a0, a1, a2, b1, b2, w
    fitresult = cfit(f, a_0, a_1, a_2, b_1, b_2, w_fixed);           
else

    fitresult = cfit(f, 0, 0, 0, 0, 0, w_fixed);       
end

end