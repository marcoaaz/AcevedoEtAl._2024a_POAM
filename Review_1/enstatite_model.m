cd('E:\Alienware_March 22\current work\00-new code May_22\mtex-5.9.0');
startup_mtex
%%
cS = crystalShape.enstatite; %ideal crystalShape
cS_symmetry = cS.CS;

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

rI_en = refractiveIndexTensor(diag([n_beta_en  n_alpha_en  n_gamma_en]), cS_symmetry); 
rI_ofe = refractiveIndexTensor(diag([n_beta_ofe  n_alpha_ofe  n_gamma_ofe]), cS_symmetry);
rI = X_en*rI_en + (1-X_en) * rI_ofe;

% compute the optical axes
vOptical = rI.opticalAxis;
% m = vector3d(Miller({1,0,0}, {0,1,0}, {0,0,1}, cS_symmetry));
a = cS_symmetry.aAxis;
b = cS_symmetry.bAxis;
c = cS_symmetry.cAxis;

% and plot it
close all

arrowWidth = 0.001;

plot(cS, 'colored', 'FaceAlpha', 0.7)
hold on
arrow3d(vOptical,'antipodal','facecolor','red', 'DisplayName', 'optic-axes')
arrow3d(a, 'facecolor', 'blue', 'DisplayName', 'a', 'arrowWidth', arrowWidth)
arrow3d(b, 'facecolor', 'green', 'DisplayName', 'b','arrowWidth', arrowWidth)
arrow3d(c, 'facecolor', 'purple', 'DisplayName', 'c','arrowWidth', arrowWidth)
hold off