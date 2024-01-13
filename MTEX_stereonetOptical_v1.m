%% Root folder

workingDir2 = 'E:\paper 2_datasets\nikon LV1000ND\mylonite\final stacks\all_modalities_rayTracing'; %optical
cd(workingDir2)
destFolder_rt = workingDir2;

%Import files
%explanation: https://mtex-toolbox.github.io/ODFTheory.html
%https://mtex-toolbox.github.io/PoleFigurePlot.html
%theory discretization: https://mtex-toolbox.github.io/PoleFigure2ODF.html

%Optical data: generating crystals
csvName = 'wave2_TL XPL-lambda_0.5_objNet_quartz.csv';
% csvName = 'wave2_TL XPL-lambda_0.5_pxNet_quartz.csv';
fileName = [workingDir2 filesep csvName]; 
fileName2 = strrep(fileName, '.csv', '_modified.csv');
[~, csvName2] = fileparts(fileName2);

%medicine
temp = readtable(fileName); %1=obj, 2=azimuth, 3=dip
n_objects = size(temp, 1);
temp2 = addvars(temp(:, [2, 3]), zeros(n_objects, 1));
temp2{:, 'Var3'} = 90 - temp2{:, 2}; %dip inversion
writetable(temp2, fileName2);

% plotting convention
setMTEXpref('xAxisDirection','south'); %setMTEXpref('yAxisDirection','north');
setMTEXpref('zAxisDirection','intoPlane');

cs = crystalSymmetry.load('quartz.cif');
h = Miller({0,0,1}, cs); %interrogated planes {1,0,0}, {0,1,0}, {0,0,1}
ori = orientation.load(fileName2, cs, 'ColumnNames', {'phi1', 'Phi', 'phi2'}); %odf
odf = calcDensity(ori, 'halfwidth', 2*degree); %orientation density function; 'zero_range', 
pf = calcPoleFigure(odf, h, equispacedS2Grid('antipodal')); %simulate diffraction counts

% Optical data plots

hFig = figure;

%odf or ori works; 'contour'; 'antipodal'; 'logarithmic'
plotPDF(odf, h, ...
    'points', 'all', 'resolution', 1.5*degree,'linewidth', 1, ...
    'grid', 'grid_res',15*degree, 'DisplayName', 'C-axis') 
annotate(zvector,'label',{'Z'},'BackgroundColor','w')
mtexColorbar
CLim(gcm,'equal');
% mtexColorMap
plot(ori, 'MarkerSize', 1, 'MarkerColor', 'magenta', 'add2all')

%Optional:
%great circles, longitude (Note: difficult to calibrate)
% f = fibre(Miller({1, 1, -2, 1}, cs), vector3d.Y); 
% plot(f,'color','red','linewidth',2,'add2all','DisplayName', 'fibre')
%Annotating segmented object names (Note: takes long)
% for i = 1:n_objects
%     label_temp = num2str(num2str(temp.Var1(i)));
%     double_temp = temp2{i, :};
%     mod_temp = orientation.byEuler(double_temp(1)*degree, ...
%         double_temp(2)*degree, ...
%         double_temp(3)*degree,cs);
% 
%     annotate(mod_temp,...
%         'marker','s','MarkerSize', 0.1, 'MarkerFaceColor', 'r',...
%         'label', label_temp, 'color', 'black', 'FontSize', 0.3)
% end
%Pole plot (Note: the data is not real)
% plotIPDF(odf,xvector,'noLabel');
% hold all % keep plot
% plot(Miller(0,0,0,1,cs),'symmetrised','labeled','backgroundColor','w')
% plot(Miller(1,1,-2,0,cs),'symmetrised','labeled','backgroundColor','w')
% plot(Miller(0,1,-1,0,cs),'symmetrised','labeled','backgroundColor','w')
% plot(Miller(0,1,-1,1,cs),'symmetrised','labeled','backgroundColor','w')
% hold off % next plot command deletes all plots

ax = gca;
hFig.WindowState = 'maximized'; 
imageFile2 = strcat(csvName2, '.pdf');
fulldest2 = fullfile(destFolder_rt, imageFile2); 
exportgraphics(ax, fulldest2, 'ContentType', 'vector')


