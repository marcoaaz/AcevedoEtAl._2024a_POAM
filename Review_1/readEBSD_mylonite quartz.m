%% Install MTEX

cd('E:\Alienware_March 22\current work\00-new code May_22\mtex-5.9.0');
startup_mtex

%%
% clear; 
% close all;

workingDir = 'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\EBSD_mylonite';
cd(workingDir)

% plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('yAxisDirection','north'); 
setMTEXpref('zAxisDirection','intoPlane');

path1 = 'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\EBSD_mylonite'; %ebsd
path2 = 'E:\paper 2_datasets\nikon LV1000ND\granulite 6KB-67\composite\all_modalities_rayTracing'; %optical

destFolder_rt = path2;

%Saving EBSD grains table
opticalDir = 'E:\paper 2_datasets\nikon LV1000ND\mylonite\final stacks';
fileDest = fullfile(opticalDir, 'ebsdPolar.csv');

%% Importing Blended image and XFM

%XFM data
DX = imread('DX_registered.tiff');
DY = imread('DY_registered.tiff');
elasticPeak = imread('elasticPeak_registered.tiff');

magnitude = sqrt(DX.^2 + DY.^2);
angle = atan2(DY, DX);
angle(angle < 0) = 2*pi + angle(angle < 0);
angle2 = angle*(180/pi);

img_h = rescale(angle, 0, 1);
img_s = rescale(elasticPeak, 0, 1);
img_v = rescale(magnitude, 0, 1);
img_hsv = cat(3, img_h, img_s, img_v);
img_rgb = uint8(255*hsv2rgb(img_hsv));
imwrite(img_rgb, 'img_rgb.tif')

%optical data
optical_reg = 'blended_registeredHSV.tif';
img_optical = imread(optical_reg);
% img_optical = imrotate(img_optical, 270); %counter-clockwise

%% Import EBSD data (Christoph data)

minSize = 80; %100
n_smooth_iter = 50;
degreeBottom = 1;
degreeTop = 2;

clear crystalShape

fname = [path1 '\ROI_Mylonite.crc'];

% crystal symmetry (trigonal and hexagonal)
%https://mtex-toolbox.github.io/CrystalReferenceSystem.html

% cs_real = {'notIndexed', ...
%     crystalSymmetry('321', [4.9 4.9 5.5], 'X||a*', 'Y||b', 'Z||c*', ...
%   'mineral', 'Quartz-new', 'color', [0.53 0.81 0.98])}; %*.cif file

cs_real = {'notIndexed', crystalSymmetry.load('quartz')};

ebsd = EBSD.load(fname, cs_real,'interface','crc',...
  'convertEuler2SpatialReferenceFrame');

ebsd('quartz').CS = ebsd('quartz').CS.Laue; %Switches to trigonal symmetry
odf_real = ebsd('indexed').orientations;

%Specifications and calculations
N = Miller({1,0,0},{0,1,0},{0,0,1}, ebsd('quartz').CS); %face normals
N_plot = Miller({0,0,1}, ebsd('quartz').CS);
% dist = [1, 1, 1]; %face distances

% crystal_shape = crystalShape( N ./ dist);
crystal_shape = crystalShape.quartz; %define a crystal shape wtih ideal model (run only once)
% crystal_shape = crystalShape.hex(ebsd.CS);

[grains, ebsd.grainId] = calcGrains(ebsd('indexed'), ...
    'threshold', [degreeBottom*degree, degreeTop*degree]); %[1*degree, 15*degree]

g30 = grains(grains.grainSize > minSize); %Ignore grains with size <= 30 pixels
g30 = smooth(g30, n_smooth_iter);
g80 = g30(g30.grainSize > 1.5*minSize);

stereonet_data = g80('quartz').meanOrientation;
p = stereonet_data*N_plot;
plot(p, 'grid')
% [theta, rho] = polar(p); %data
% theta= Z towards stereonet periphery; 
% rho= X clockwise

temp_ID = [g80.id];
temp_grainSize = [g80.grainSize];
temp_phi1 = [stereonet_data.phi1];
temp_Phi1 = [stereonet_data.Phi];
temp_phi2 = [stereonet_data.phi2];
temp_GOS = [g80.GOS];
temp_theta = [p.theta]; %of Miller direction
temp_rho = [p.rho];
temp_x = [p.x];
temp_y = [p.y];
temp_z = [p.z];

%Save ebsd data
vars_temp = [temp_ID, temp_grainSize, ...
    temp_phi1, temp_Phi1, temp_phi2, temp_GOS, ...
    temp_x, temp_y, temp_z, temp_theta, temp_rho];
varNames_temp = {'ID', 'grainSize', ...
    'phi1', 'Phi1', 'phi2', 'GOS', ...
    'x', 'y', 'z', 'theta', 'rho'};
temp_ebsdPolar = array2table(vars_temp, 'VariableNames', varNames_temp);

writetable(temp_ebsdPolar, fileDest)

%% EBSD Grain boundaries with XFM 2D Quiver

%crystal shape for plot
cSGrains = g80('quartz').meanOrientation *...
    crystal_shape * 0.7 * sqrt(g80('quartz').area); 

INC = 4;
R = 1:INC:length(DX);
[x, y]= meshgrid(1:250,1:250);

%settings
alpha = g80.innerBoundary.misorientation.angle / (5*degree);
autoScale = 1;

%Grain boundaries with ideal crystal
close all

hFig = figure; %transposed/permuted image (it is not real)

plot(g80.boundary, 'linewidth', 2, 'micronbar', 'off', 'linecolor', [1 0 1]);
hold on
% plot(g80.innerBoundary, 'linewidth',1.5, 'edgeAlpha',alpha, 'linecolor',[0. 0.2 .7]);

imshow(img_optical)
% quiver(x(R,R), y(R,R), DX(R,R), -(DY(R,R)), autoScale, 'LineWidth', 1, 'Color', 'yellow'); 

plot(ebsd('indexed'), ebsd('indexed').orientations, 'faceAlpha', 0.5, 'figSize','large')
plot(g80('quartz'), 0.4*crystal_shape, 'linewidth', 2, 'colored', 'FaceAlpha', 0); %0.6 (higher is opaque)
hold off
text(g80, g80.id)

title('Grain boundaries with EBSD model and XFM 2D Quiver')
legend off;
hFig.WindowState = 'maximized'; 

%Stereonet
%https://mtex-toolbox.github.io/CrystalShapes.html

hFig2 = figure;

%for N_plot, 'antipodal'
plotPDF(stereonet_data, N, 'contour', 'resolution', 1.5*degree, 'linewidth', 1, ...
    'grid', 'grid_res', 15*degree, 'DisplayName', 'C-axis') %,'linewidth',2
plot(stereonet_data, 0.006*cSGrains, 'add2all', 'linewidth', 1, 'FaceAlpha', 0.3); %shapes

mtexColorMap('turbo')
mtexColorbar % remove colorbars
mtexColorbar % remove colorbars
CLim(gcm,'equal');
mtexColorbar % add a single colorbar

hFig2.WindowState = 'maximized'; 

%% Old optical stereonets (now plotMTEXstereonet.m)
%explanation: https://mtex-toolbox.github.io/ODFTheory.html
%https://mtex-toolbox.github.io/PoleFigurePlot.html
%theory discretization: https://mtex-toolbox.github.io/PoleFigure2ODF.html

% Import files 
fileName = [path2 filesep 'wave2_TL XPL-lambda_0.5_targetStereo_gt.csv']; 
fileName2 = strrep(fileName, '.csv', '_modified.csv');
[~, csvName] = fileparts(fileName2);

%medicine
temp = readtable(fileName); %1=obj, 2=azimuth, 3=dip
n_objects = size(temp, 1);
temp2 = addvars(temp(:, [2, 3]), zeros(n_objects, 1));
temp2{:, 'Var3'} = 90-temp2{:, 2}; %dip inversion
writetable(temp2, fileName2);

cs = crystalSymmetry.load('quartz.cif');
h = Miller({0,0,1}, cs); %interrogated planes {1,0,0}, {0,1,0}, {0,0,1}
ori = orientation.load(fileName2,cs,'ColumnNames', {'phi1','Phi','phi2'}); %odf
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
imageFile2 = strcat(csvName, '.pdf');
fulldest2 = fullfile(destFolder_rt, imageFile2); 
exportgraphics(ax, fulldest2, 'ContentType', 'vector')

%% Old plots
% Synchrotron plot

% figure, 
% binscatter(DX(:), DY(:))
% xline(0, 'LineWidth', 1, 'Color', 'black')
% yline(0, 'LineWidth', 1, 'Color', 'black')
% grid on

% EBSD data plots 
alpha = g30.innerBoundary.misorientation.angle / (5*degree);

% %Grain boundaries with Good pixels
% figure;
% plot(ebsd('indexed'), ebsd('indexed').orientations, 'faceAlpha',1, 'figSize','large')
% hold on;
% plot(g30.boundary, 'linewidth',2);
% plot(g30.innerBoundary, 'linewidth', 1.5, 'edgeAlpha', alpha, 'linecolor',[0. 0.2 .7]);
% hold off
% title('Grain boundaries with Good pixels')

% %Pole figure
% ipfKey = ipfColorKey(ebsd('Quartz-new'));
% colors = ipfKey.orientation2color(g30('Quartz-new').meanOrientation);

% figure;
% plotIPDF(g30('Quartz-new').meanOrientation, colors, vector3d.Z,...
% 'MarkerSize', 0.1*g30('Quartz-new').area, 'markerEdgeColor', 'k');
% title('Pole figure')

