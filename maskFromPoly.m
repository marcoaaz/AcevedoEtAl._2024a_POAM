function [label_map] = maskFromPoly(S, img_file, cmap, plotOption)
%reads the structure from shaperead() and produces semantic segmentation
%masks

n_rows = imfinfo(img_file).Height;
n_cols = imfinfo(img_file).Width;
n_poly = length(S);

%Create masks (~3 min) (from qupathPhaseMap_v6.m)
label_map = zeros(n_rows, n_cols); %can be any spatial resolution
for i = 1:n_poly
    X = S(i).X;
    Y = S(i).Y;
    stopX = find(isnan(X)) - 1;
    stopY = find(isnan(Y)) - 1;
    temp = poly2mask(X(1:stopX), Y(1:stopY), n_rows, n_cols);%not NaN  
    
    %assign label
    temp_ID = S(i).ClassID + 1; %to avoid indexing issues    
    label_map(temp) = temp_ID;    
end

if plotOption == 1
    %semantic segmentation
    figure    
    imshow(label2rgb(label_map, cmap))
end