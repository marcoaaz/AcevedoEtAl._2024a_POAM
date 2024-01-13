function [s] = sinDescriptor(sf)

%Input:
%sf = cfit

%Output:
%s = structure containing 

%optical period in sexagesimal degrees
w = sf.w;
period = 2*pi/w;
search_interval = 30; %minimum expected distance between features
x0_range = 0:search_interval:period;
n_x0 = length(x0_range);

s = struct; %pre-allocating

%Function: cfit to sym
formula_str = strrep(formula(sf), newline, '');
f = subs(str2sym(formula_str), coeffnames(sf), num2cell(coeffvalues(sf).'));

%Derivates
vars = sym(indepnames(sf));
d1_sf = diff(f, vars(1));
d2_sf = diff(d1_sf, vars(1));

if d1_sf == 0
    A = [0, 0, 0]; %inflexion points
    B2 = [0, 0;
        0, 0]; %critical points
    range1 = 0;
    range2 = 0;
    avg_f = 0;
    avg_x = 0;

    %prevents FZERO cannot continue because user-supplied function_handle
    %==> @()0.0 failed with the error below. Too many input arguments.
else

    %Fitting Sin: input is a function handle
    options = optimset('Display', 'off'); %'iter'
    
    zero_critical_list = zeros(1, n_x0, 'double');
    zero_inf_list = zeros(1, n_x0, 'double');
    f_inf_list = zeros(1, n_x0, 'double');
    d1_inf_list = zeros(1, n_x0, 'double');    
    f_critical_list = zeros(1, n_x0, 'double');
    k = 0;
    for x0_guess = x0_range
        k = k + 1;
        
        %derivatives
        [zero_critical_temp, ~, ~, ~] = fzero(matlabFunction(d1_sf), x0_guess, options);
        [zero_inf_temp, ~, ~, ~] = fzero(matlabFunction(d2_sf), x0_guess, options);
        %[x_zero3, fval3, exitflag, output]
            
        %first
        % Note: in TL-XPL-lambda, the min peak is the second max peak.
        f_critical_temp = feval(sf, zero_critical_temp);
        
        %second
        f_inf_temp = feval(sf, zero_inf_temp);
        d1_inf_temp = double(subs(d1_sf, zero_inf_temp));    
            
        %attributing
        zero_critical_list(k) = zero_critical_temp;
        zero_inf_list(k) = zero_inf_temp;
        f_critical_list(k) = f_critical_temp;
        f_inf_list(k) = f_inf_temp;
        d1_inf_list(k) = d1_inf_temp;            
    end
    
    %medicine
    zero_inf_list(zero_inf_list < 0) = zero_inf_list(zero_inf_list < 0) + period; %not negative
    zero_inf_list(zero_inf_list >= 180) = zero_inf_list(zero_inf_list >= 180) - period; %not above 180
    zero_critical_list(zero_critical_list < 0) = zero_critical_list(zero_critical_list < 0) + period; %not negative
    zero_critical_list(zero_critical_list >= 180) = zero_critical_list(zero_critical_list >= 180) - period; %not above 180
    
    info_function = [x0_range', zero_inf_list', f_inf_list', d1_inf_list', ...
        zero_critical_list', f_critical_list'];
    
    tolerance1 = 0.005;
    [~, idx_A] = uniquetol(info_function(:, [2, 3]), tolerance1, 'ByRows', true'); %stable
    [~, idx_B] = uniquetol(info_function(:, [5, 6]), tolerance1, 'ByRows', true');
    
    A = info_function(idx_A, [2, 3, 4]);
    B = info_function(idx_B, [5, 6]);
    A2 = sortrows(A, 3, 'descend'); %+ to -
    B2 = sortrows(B, 2, 'descend'); %max to min   
   
    %Range
    rows_temp = size(B2, 1);
    if rows_temp == 4
        range1 = B2(1, 2) - B2(end, 2);
        range2 = B2(2, 2) - B2(end-1, 2);    
    elseif rows_temp == 3
        range1 = B2(1, 2) - B2(end, 2);
        range2 = B2(1, 2) - B2(end-1, 2);    
    else
        range1 = B2(1, 2) - B2(end, 2);
        range2 = 0;    
    end

    positive_inf = A2(1, :);% inflexion of positive slope (f zero intercept)
    avg_f = mean(A2(:, 2)); %mean of inflexions
    avg_x = fzero(matlabFunction(f - avg_f), positive_inf(1), options);
    avg_x(avg_x < 0) = avg_x(avg_x < 0) + period; %medicine

end

%Attributing
s.avg_f = avg_f;
s.avg_x = avg_x;
s.range1 = range1;
s.range2 = range2;
s.Inflexion = A; %inflexion points
s.MinMax = B2; %critical points

end