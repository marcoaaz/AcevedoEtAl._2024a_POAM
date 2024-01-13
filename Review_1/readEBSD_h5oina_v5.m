%% Install MTEX

cd('E:\Alienware_March 22\current work\00-new code May_22\mtex-5.9.0');
startup_mtex

%MTEX plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('yAxisDirection','north'); 
setMTEXpref('zAxisDirection','intoPlane');

%Relevant documentation
%https://mtex-toolbox.github.io/GrainTutorial.html
%https://mtex-toolbox.github.io/EBSD.plot.html
%https://mtex-toolbox.github.io/EBSD.quiver.html
%https://mtex-toolbox.github.io/BirefringenceDemo.html

%%
clear; 
close all;

%Quartz
% workingDir = 'E:\Alienware_March 22\02-Research geology\05-QUT\Paper 2\EBSD_mylonite';
% fileName1 = [workingDir '\ROI_Mylonite.crc'];

%Enstatite
workingDir = 'E:\Alienware_March 22\02-Research geology\05-QUT\Paper 2\EBSD_harzburgite 17-BSK-035_opx';
fileName1 = '17BSK_035 Aggressive-035Subset.h5oina';
cd(workingDir)


ebsd = EBSD.load(fileName1, 'convertEuler2SpatialReferenceFrame');
ebsd.mineralList'

%% Reconstruct grains

min_size = 15; %pixels
n_smooth_iter = 3;
degreeBottom = 3;
degreeTop = 5;

[grains, ebsd.grainId] = calcGrains(ebsd('indexed'), ...
    'threshold', [degreeBottom*degree, degreeTop*degree]);

grains = smooth(grains, n_smooth_iter); %data denoising
grains_min1 = grains(grains.grainSize > min_size);
grains_min1_smoothed = smooth(grains_min1, n_smooth_iter);
grains_min2 = grains_min1_smoothed(grains_min1_smoothed.grainSize > 1.5*min_size);

%Phase map
lineWidth = 1;

% figure
% plot(ebsd)
% hold on 
% plot(grains_min2.boundary, 'lineWidth', lineWidth)
% hold off

%Selected mineral
sel_mineral = 'Enstatite'; %'Quartz-new'; 'Enstatite'
ori_px = ebsd(sel_mineral).orientations;%pixel crystallographic orientations
target = grains_min2(sel_mineral);
ori = target.meanOrientation;%grain crystallographic orientations
[centroid_x, centroid_y] = centroid(target);

%Data for table of EBSD objects for validation
temp_grain_id = target.id;
temp_grainSize = [target.grainSize];
temp_phi1 = [ori.phi1];
temp_Phi1 = [ori.Phi];
temp_phi2 = [ori.phi2];
temp_GOS = [target.GOS];

%Experiment geometry
vprop = vector3d.Z; % the propagation direction
omega = 35; % the polarizer direction
polarizer = rotate(vector3d.X, omega * degree); 
thickness = 100000; %nm

%Target mineral
X_en = 0.92; %Mol% (atomic %)
X_ofe = 0.08; %including 0.01 wollastonite

%Troger 1980 (book, pg.72, 76, 78)
%Enstatite
n_alpha_en = 1.657; %n_x
n_beta_en = 1.659; %n_y
n_gamma_en = 1.665; %n_z

%Ortho-ferrosilite 
n_alpha_ofe = 1.765;
n_beta_ofe = 1.770; 
n_gamma_ofe = 1.788;

%Crystallographic model (a, b, c)
cs = ebsd(sel_mineral).CS; %real crystalSymmetry (objects)
cs_ideal = crystalShape.enstatite; %ideal crystalShape

rI_en = refractiveIndexTensor(diag([n_beta_en  n_alpha_en  n_gamma_en]), cs); 
rI_ofe = refractiveIndexTensor(diag([n_beta_ofe  n_alpha_ofe  n_gamma_ofe]), cs);
rI = X_en*rI_en + (1-X_en) * rI_ofe;
rISpecimen = ori*rI; % transform tensor into specimen coordinates
[dn, pMin, pMax] = rISpecimen.birefringence(vprop); %birefringence vectors
%theta, rho = polar_angle, azimuth_angle
%pMax = fastest polarization direction
%pMin = slowest polarization direction

%Fix (mirroring hemisphere)
vertical_dir = pMin.y;
pMin.y = -1*vertical_dir;

%optic-axis orientation
vOptical = rISpecimen.opticalAxis; %propagation vector
x_oa = vOptical.x; %u unitary vector
y_oa = vOptical.y; %v
z_oa = vOptical.z; 
x_oa1 = x_oa(:, 1);
y_oa1 = y_oa(:, 1);
z_oa1 = z_oa(:, 1);
x_oa2 = x_oa(:, 2);
y_oa2 = y_oa(:, 2);
z_oa2 = z_oa(:, 2);

%Grain slow-axis orientation
x = pMin.x; %u unitary vector
y = pMin.y; %v

%Save ebsd data
vars_temp = [temp_grain_id, temp_grainSize, ...
    temp_phi1, temp_Phi1, temp_phi2, temp_GOS, ...
    centroid_x, centroid_y, x, y];
varNames_temp = {'ID', 'grainSize', ...
    'phi1', 'Phi1', 'phi2', 'GOS', ...
    'centroid_x', 'centroid_y', 'x', 'y'};
temp_ebsdPolar = array2table(vars_temp, 'VariableNames', varNames_temp);

fileDest = fullfile(workingDir, 'ebsd_objects.csv');
writetable(temp_ebsdPolar, fileDest)

%% Plot

fontSize = 10;

colorKey  = spectralTransmissionColorKey(rI,thickness);
colorKey.propagationDirection = vprop; %vertical
colorKey.polarizer = []; %[] for circular; polarizer for linear (E-W)
colorKey.phi = 90 * degree; %XPL=90 to PPL=0 (rotation from the analyzer)
%not used: tau, angle between polarization direction and the slowest
%polarization direction pMin

rgb = colorKey.orientation2color(ori_px);

hFig = figure;

%assemblage
plot(grains( ...
    {'Forsterite, ferroan', ...
    'Pyrope',...
    'Diopside', ...
    'Chromite'}), ...
    'FaceAlpha', .3, 'EdgeAlpha', 0.2)
hold on

%target
plot(ebsd(sel_mineral), rgb) %pixels according to spectral transmission
plot(grains_min2(sel_mineral).boundary, 'lineWidth', lineWidth) %grain boundaries
% text(grains_min2(sel_mineral), grains_min2(sel_mineral).id, ...
%     'FontSize', fontSize, 'Color', 'black', ...
%     'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
%     'BackgroundColor', 'w')

%slow-axis
plot(centroid_x, centroid_y, '.', 'MarkerSize', 10, ...
    'Color', 'red', 'DisplayName', 'Centroid')
quiver(centroid_x, centroid_y, x, y, ...
    'LineWidth', 2, 'AutoScaleFactor', .2, ...
    'Color', 'red', 'DisplayName', 'Slow-axis')

%optic-axes (projection into 2D plane)
length_val = 200;
% quiver(centroid_x, centroid_y, length_val*x_oa1, length_val*y_oa1, ...
%     'LineWidth', 2, 'AutoScale', 'off', 'AutoScaleFactor', 1, ...
%     'Color', 'red', 'DisplayName', 'Optic-axis 1')
% quiver(centroid_x, centroid_y, length_val*x_oa2, length_val*y_oa2, ...
%     'LineWidth', 2, 'AutoScale', 'off', 'AutoScaleFactor', 1, ...
%     'Color', 'green', 'DisplayName', 'Optic-axis 2')

%ideal crystal
% plot(grains_min2(sel_mineral), cs_ideal, 'colored', 'linewidth', 1, 'FaceAlpha', 0.5) 

hold off
legend('Location', 'eastoutside')

% hFig.WindowState = 'maximized';

