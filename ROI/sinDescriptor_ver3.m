function [s] = sinDescriptor_ver3(x0_range, F, F_d1, F_d2)

%Input:
%sf = structure containing algorithm 1 model (cfit)
%F = Fourier-2 function handle
%F_d1 = first derivative of F
%F_d2 = second derivative of F

%Output:
%s =  structure containing algorithm 2 description

%%
search_interval = x0_range(2) - x0_range(1);
s = struct; %pre-allocating

options = optimset('Display', 'off'); %'iter'

info_A = [];
info_B = [];
for bottom = x0_range
    
    top = bottom + search_interval;

    %derivative intercepts in X axis (Y must differ in sign)
    try 
        [zero_critical_temp, ~, ~, ~] = fzero(F_d1, ...
            [bottom, top], options);

        %F at first derivative intercept                
        f_critical_temp = F(zero_critical_temp);
        % Note: in TL-XPL-lambda, the min peak is the second max peak.        
        
        info_B = [info_B; zero_critical_temp, f_critical_temp];
    end

    try
        [zero_inf_temp, ~, ~, ~] = fzero(F_d2, ...
            [bottom, top], options);
        
        %F at second derivative intercept        
        f_inf_temp = F(zero_inf_temp);
        d1_inf_temp = F_d1(zero_inf_temp);
        
        info_A = [info_A; zero_inf_temp, f_inf_temp, d1_inf_temp];
    end           
            
end 
A2 = sortrows(info_A, 3, 'descend'); %inflections (+ to - slope)
B2 = sortrows(info_B, 2, 'descend'); %max to min   
% (assumes zero-search redundancy at intervals is very unlikely)    

%Spectra range (preventing missing data issues)
rows_temp = size(B2, 1);
if rows_temp == 4
    range1 = B2(1, 2) - B2(end, 2);
    range2 = B2(2, 2) - B2(end-1, 2);    
elseif rows_temp == 3    
    range1 = B2(1, 2) - B2(end, 2);
    range2 = B2(1, 2) - B2(end-1, 2);    
else % == 2    
    range1 = B2(1, 2) - B2(end, 2);
    range2 = 0;    
end

avg_f = mean(A2(:, 2)); %mean of inflexions
positive_inf_x = A2(1, :); %inflexion of largest positive slope (f zero intercept)

F2 = @(x)F(x) - avg_f;
avg_x = fzero(F2, positive_inf_x(1), options);
%Note: avg_x can fail if the avg_f is much lower than f due to an artefact

%Optinal medicine (or fix outside):
avg_x(avg_x < 0) = avg_x(avg_x < 0) + 180;
avg_x(avg_x >= 180) = avg_x(avg_x >= 180) - 180;

%Attributing
s.avg_f = avg_f;
s.avg_x = avg_x;
s.range1 = range1;
s.range2 = range2;
s.Inflexion = A2; %inflexion points
s.MinMax = B2; %critical points

%  %debug:
% info_A
% info_B
% A2
% B2

end