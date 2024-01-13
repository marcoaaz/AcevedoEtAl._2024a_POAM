function [annotationNames, cmap] = SuperSIATinterpreter(temp_name)
%reads SuperSIAT annotation metadata
%outputs the class names of the annotations, color map, and optionally a
%table with the training points

fid = fopen(temp_name); 
dd_cell = textscan(fid,'%s');
dd = dd_cell{1}; %1x1 cell
fclose(fid);
n_class = str2num(dd{1});

fileID = fopen(temp_name, 'r');
text_char = fscanf(fileID, '%c');

% expression1 = '(?<annotation>[a-zA-Z_-]+)';
expression1 = '(?<annotation>[a-zA-Z][a-zA-Z_-\s]+[a-zA-Z][0-9]?)'; %labels and special symbols '_'
expression4 = '\n+'; %data with new line char(10)

n_classes = str2num(text_char(1:2)); %2-digits
annotationNames = string(regexp(text_char, expression1, 'tokens'));
annotationNames' %display

[idxStart, idxEnd] = regexp(text_char, expression1);

annotationArray = [];
training_dataArray = [];
colour_val = zeros(1, n_classes, 'double');
population_val = zeros(1, n_classes, 'double');
for i = 1:n_classes
    
    %data
    dataStart_temp = idxEnd(i) + 2;
    if i < n_classes
        dataEnd_temp = idxStart(i+1) - 2;
    elseif i == n_classes
        dataEnd_temp = length(text_char) - 2;
    end
    text_temp = text_char(dataStart_temp:dataEnd_temp);    

    from1 = regexp(text_temp, expression4, 'start');    
    
    text_temp %debug: help line 

    %description
    a = from1(1) + 1;
    b = from1(2) - 1; 
    %Error: Index exceeds the number of array elements. Index must not exceed 1.
    %Solution: every class must have at least 1 annotation point
    
    text1 = text_temp(a:b);
    text1 = str2double( strsplit(string(text1)) );
    colour_val(i) = text1(1);
    n_infoRows = text1(2);    
    population_val(i) = n_infoRows;

    %data
    training_data = zeros(n_infoRows, 3, 'double');
    for j = 1:n_infoRows
        c = from1(j + 1) + 1;
        if j < n_infoRows
            d = from1(j + 2) + 1;
        elseif j == n_infoRows
            d = length(text_temp) - 1;
        end
        text2 = text_temp(c:d);
        text2 = str2double( strsplit(string( text2 )) );
        training_data(j, 1) = text2(1); %x
        training_data(j, 2) = text2(2); %y
        training_data(j, 3) =  text2(3); %scale 
    end
    training_dataArray = [training_dataArray; training_data];  

    %labels
    annotations_col = repmat(annotationNames(i), n_infoRows, 1);
    annotationArray = [annotationArray; annotations_col];

end
fclose(fileID);

%Colouring (from qupathPhaseMap_v6.m)
dataType = 'int32'; %'int32'
a = int32(colour_val'); %int32()
B = bitand(a, 255, dataType);
G = bitand(bitsrl(a, 8), 255, dataType);
R = bitand(bitsrl(a, 16), 255, dataType);
alpha = bitand(bitsrl(a, 24), 255, dataType);

rgb_append = double([B, G, R, alpha]); %inverse order (SuperSIAT convention)
rgb_append = array2table(rgb_append);
rgb_append.Properties.VariableNames = {'R', 'G', 'B', 'alpha'};
cmap = rgb_append{:, 1:3}/255;


end