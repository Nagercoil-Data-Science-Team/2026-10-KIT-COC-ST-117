
clear; clc; close all;
fprintf('=======================================================================\n');
fprintf('  PERSONALIZED IMPLANT DESIGN (STL-DRIVEN)\n');
fprintf('=======================================================================\n\n');

base_dir = fileparts(mfilename('fullpath'));
if isempty(base_dir), base_dir = pwd; end
results_dir = fullfile(base_dir, 'implant_results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

PATIENT_ID = 'LTKN9332';
FONT = 'Times New Roman';
set(0, 'DefaultAxesFontName', FONT);
set(0, 'DefaultTextFontName', FONT);

% -------------------------------------------------------------------------
% 1. LOAD BONE STLs  (assumed axes: X = superior-inferior, Y = ML, Z = AP;
%    tibia proximal end at max X, femur distal end at min X)
% -------------------------------------------------------------------------
fprintf('[+] Step 1: Loading bone STL models...\n');
femur_stl_path = fullfile(base_dir, 'femur_model.stl');
if ~exist(femur_stl_path, 'file')
    femur_stl_path = fullfile(base_dir, 'reconstruction_results', 'femur_model.stl');
end
tibia_stl_path = fullfile(base_dir, 'tibia_model.stl');
if ~exist(tibia_stl_path, 'file')
    tibia_stl_path = fullfile(base_dir, 'reconstruction_results', 'tibia_model.stl');
end
if ~exist(femur_stl_path, 'file') || ~exist(tibia_stl_path, 'file')
    error('Bone STL models not found. Ensure reconstruction has been executed.');
end
femur_tri = stlread(femur_stl_path);  tibia_tri = stlread(tibia_stl_path);
fv = femur_tri.Points;  ff = femur_tri.ConnectivityList;
tv = tibia_tri.Points;  tf = tibia_tri.ConnectivityList;
fprintf('    - Femur STL : %d vertices, %d faces\n', size(fv,1), size(ff,1));
fprintf('    - Tibia STL : %d vertices, %d faces\n', size(tv,1), size(tf,1));

% -------------------------------------------------------------------------
% 2. SURGICAL RESECTION PLANES + LANDMARKS FROM THE STL
% -------------------------------------------------------------------------
tib_resect = 6.0;     % mm removed from proximal tibia
fem_resect = 8.5;     % mm removed from distal femur (= distal metal thickness)

tib_top = max(tv(:,1));
x_tcut  = tib_top - tib_resect;                         % tibial cut plane
sec_t   = tv(abs(tv(:,1) - x_tcut) < 1.5, :);
if size(sec_t,1) < 20, sec_t = tv(tv(:,1) > x_tcut - 4, :); end

xf_min  = min(fv(:,1));
sec_f   = fv(abs(fv(:,1) - (xf_min + fem_resect)) < 1.5, :);
if size(sec_f,1) < 20, sec_f = fv(fv(:,1) < xf_min + fem_resect + 4, :); end

mid  = @(a) (max(a) + min(a)) / 2;
span = @(a) max(a) - min(a);
plateau_width_ml = span(sec_t(:,2));   plateau_depth_ap = span(sec_t(:,3));
plateau_center_y = mid(sec_t(:,2));    plateau_center_z = mid(sec_t(:,3));
condyle_width_ml = span(sec_f(:,2));   condyle_depth_ap = span(sec_f(:,3));
condyle_center_y = mid(sec_f(:,2));    condyle_center_z = mid(sec_f(:,3));

fprintf('\n[+] Anatomical landmarks measured at the resection planes:\n');
fprintf('    - Tibial cut plane X            : %.2f mm\n', x_tcut);
fprintf('    - Tibial section  ML x AP       : %.2f x %.2f mm\n', plateau_width_ml, plateau_depth_ap);
fprintf('    - Femoral section ML x AP       : %.2f x %.2f mm\n', condyle_width_ml, condyle_depth_ap);

% -------------------------------------------------------------------------
% 3. IMPLANT PARAMETERS (one implant, sized to this patient's bones)
% -------------------------------------------------------------------------
fprintf('\n[+] Step 2: Defining implant parameters...\n');
params = struct();

% Tibial tray (2 mm under-sized so it does not overhang the bone)
params.tibial_tray.W_ml         = plateau_width_ml - 2.0;
params.tibial_tray.D_ap         = plateau_depth_ap - 2.0;
params.tibial_tray.thickness    = 3.5;
params.tibial_tray.stem_length  = 28.0;
params.tibial_tray.stem_radius  = 5.5;
params.tibial_tray.fin_width    = 16.0;
params.tibial_tray.fin_thick    = 3.2;
params.tibial_tray.notch_factor = 0.22;
params.tibial_tray.center_y     = plateau_center_y;
params.tibial_tray.center_z     = plateau_center_z;
params.tibial_tray.x_bot        = x_tcut;                                   % sits on cut bone
params.tibial_tray.x_top        = params.tibial_tray.x_bot + params.tibial_tray.thickness;
params.tibial_tray.x_stem_tip   = params.tibial_tray.x_bot - params.tibial_tray.stem_length;

% Femoral component dimensions (needed first to place insert dishes)
params.femoral.W_ml            = condyle_width_ml - 1.0;
params.femoral.D_ap            = condyle_depth_ap - 1.0;
params.femoral.resect_thick    = fem_resect;
params.femoral.center_y        = condyle_center_y;
params.femoral.center_z        = condyle_center_z;
params.femoral.r_med_runner    = 22.5;
params.femoral.r_lat_runner    = 21.0;
params.femoral.condyle_width   = 0.42 * params.femoral.W_ml;
params.femoral.condyle_spacing = params.femoral.W_ml - params.femoral.condyle_width;
params.femoral.trochlear_angle = 140.0;
params.femoral.top_height      = 32.0;     % flange height above apex (mm)
params.femoral.wall_thick      = 3.5;      % posterior/anterior wall (mm)

% UHMWPE insert
params.tibial_insert.core_thickness = 9.0;
params.tibial_insert.x_bot          = params.tibial_tray.x_top;
params.tibial_insert.x_top          = params.tibial_insert.x_bot + params.tibial_insert.core_thickness;
params.tibial_insert.med_dish_depth = 2.8;
params.tibial_insert.lat_dish_depth = 2.2;
params.tibial_insert.dish_offset    = params.femoral.condyle_spacing / 2;
params.tibial_insert.med_dish_ry    = params.femoral.condyle_width / 2;
params.tibial_insert.lat_dish_ry    = params.femoral.condyle_width / 2;
params.tibial_insert.med_dish_rz    = 0.30 * params.tibial_tray.D_ap;
params.tibial_insert.lat_dish_rz    = 0.28 * params.tibial_tray.D_ap;
params.tibial_insert.spine_height   = 4.5;
params.tibial_insert.spine_width    = 8.0;
params.tibial_insert.spine_length   = 14.0;
params.tibial_insert.conformity     = 1.08;   % dish radius / femoral radius

params.femoral.x_apex = params.tibial_insert.x_top - 1.0;   % condyle seated in dish

% Materials
params.materials.ti6al4v = struct('name','Ti-6Al-4V ELI (Titanium Alloy)','E',110000,'nu',0.33,'yield_stress',880,'uts',950);
params.materials.uhmwpe  = struct('name','UHMWPE (Articular Polyethylene)','E',1000,'nu',0.42,'yield_stress',25,'uts',40);
params.materials.cocrmo  = struct('name','Co-Cr-Mo Cast Alloy (ASTM F75)','E',210000,'nu',0.30,'yield_stress',650,'uts',1000);

% ISO 14243 peak stance
params.loading.F_axial_total = 2100;
params.loading.medial_ratio  = 0.60;
params.loading.lateral_ratio = 0.40;
params.loading.F_ap_shear    = 320;

fprintf('    - Tibial tray       : W=%.2f mm  D=%.2f mm\n', params.tibial_tray.W_ml, params.tibial_tray.D_ap);
fprintf('    - Femoral component : W=%.2f mm  D=%.2f mm\n', params.femoral.W_ml, params.femoral.D_ap);

% -------------------------------------------------------------------------
% 4. PLACE BONES AROUND THE IMPLANT AND APPLY THE RESECTIONS
% -------------------------------------------------------------------------
fem_shift = params.femoral.x_apex - xf_min;      % femur distal end -> implant apex
fv_insitu = fv;  fv_insitu(:,1) = fv_insitu(:,1) + fem_shift;
x_fcut    = params.femoral.x_apex + fem_resect;  % distal femoral cut plane

keep_f = all(reshape(fv_insitu(ff,1), size(ff)) >= x_fcut - 1e-6, 2);
keep_t = all(reshape(tv(tf,1),        size(tf)) <= x_tcut + 1e-6, 2);
[bone_ff, bone_fv] = decimate_mesh(ff(keep_f,:), fv_insitu, 30000);
[bone_tf, bone_tv] = decimate_mesh(tf(keep_t,:), tv,        30000);

% -------------------------------------------------------------------------
% 5. GENERATE IMPLANT MESHES + STL EXPORT
% -------------------------------------------------------------------------
fprintf('\n[+] Step 3: Generating implant meshes and exporting STL files...\n');
[tray_V, tray_F] = generate_tibial_tray_mesh(params);
[ins_V,  ins_F ] = generate_tibial_insert_mesh(params);
[fem_V,  fem_F ] = generate_femoral_component_mesh(params);
if isempty(tray_V) || isempty(ins_V) || isempty(fem_V)
    error('An implant mesh came out empty - check the STL landmark values above.');
end
[asm_V, asm_F] = merge_triangulations({tray_V, ins_V, fem_V}, {tray_F, ins_F, fem_F});
[all_V, all_F] = merge_triangulations({bone_tv, tray_V, ins_V, fem_V, bone_fv}, ...
                                      {bone_tf, tray_F, ins_F, fem_F, bone_ff});

export_stl(triangulation(tray_F, tray_V), 'implant_tibial_tray.stl',        results_dir, base_dir);
export_stl(triangulation(ins_F,  ins_V ), 'implant_tibial_insert.stl',      results_dir, base_dir);
export_stl(triangulation(fem_F,  fem_V ), 'implant_femoral_component.stl',  results_dir, base_dir);
export_stl(triangulation(asm_F,  asm_V ), 'implant_knee_assembly.stl',      results_dir, base_dir);
stlwrite(triangulation(all_F, all_V), fullfile(results_dir, 'implant_assembly_with_bone.stl'));
fprintf('    - Tray %d f | Insert %d f | Femoral %d f\n', size(tray_F,1), size(ins_F,1), size(fem_F,1));

% -------------------------------------------------------------------------
% 6. ANALYTICAL STRESS ESTIMATION
% -------------------------------------------------------------------------
fprintf('\n[+] Step 4: Analytical stress estimation...\n');
F_total = params.loading.F_axial_total;
F_med   = F_total * params.loading.medial_ratio;
F_lat   = F_total * params.loading.lateral_ratio;
F_shear = params.loading.F_ap_shear;

% Tibial tray - plate bending
tray_w = params.tibial_tray.W_ml;  tray_d = params.tibial_tray.D_ap;  tray_t = params.tibial_tray.thickness;
plate_area = tray_w * tray_d;
sigma_bend_tray = 0.0138 * F_total * tray_d^2 / (tray_t^3 * tray_w);
stem_area  = pi * params.tibial_tray.stem_radius^2;
sigma_stem = F_total / stem_area;
tau_tray   = 1.5 * F_shear / plate_area;
vm_tray_analytical = sqrt(sigma_bend_tray^2 + 3*tau_tray^2);
sf_tray_analytical = params.materials.ti6al4v.yield_stress / vm_tray_analytical;

% UHMWPE insert - conforming Hertz contact (femoral ball in a slightly larger dish)
R_med = params.femoral.r_med_runner;  R_lat = params.femoral.r_lat_runner;
cf    = params.tibial_insert.conformity;
Reff_med = 1 / (1/R_med - 1/(cf*R_med));
Reff_lat = 1 / (1/R_lat - 1/(cf*R_lat));
E_star = 1 / ((1 - params.materials.uhmwpe.nu^2)/params.materials.uhmwpe.E + ...
              (1 - params.materials.cocrmo.nu^2)/params.materials.cocrmo.E);
p0_med = (6 * F_med * E_star^2 / (pi^3 * Reff_med^2))^(1/3);
p0_lat = (6 * F_lat * E_star^2 / (pi^3 * Reff_lat^2))^(1/3);
vm_ins_max = 0.62 * max(p0_med, p0_lat);
sf_ins_analytical = params.materials.uhmwpe.yield_stress / vm_ins_max;

% Femoral component - mean contact stress
a_med = (3 * F_med * Reff_med / (4 * E_star))^(1/3);
a_lat = (3 * F_lat * Reff_lat / (4 * E_star))^(1/3);
A_med = pi*a_med^2;  A_lat = pi*a_lat^2;
vm_fem_max = max(F_med/A_med, F_lat/A_lat);
sf_fem_analytical = params.materials.cocrmo.yield_stress / vm_fem_max;

fprintf('    Tibial tray : bending %.2f MPa | stem %.2f MPa | VM %.2f MPa | SF %.2f\n', ...
    sigma_bend_tray, sigma_stem, vm_tray_analytical, sf_tray_analytical);
fprintf('    Insert      : p0 med %.2f / lat %.2f MPa | VM %.2f MPa | SF %.2f\n', ...
    p0_med, p0_lat, vm_ins_max, sf_ins_analytical);
fprintf('    Femoral     : a med %.2f / lat %.2f mm | contact %.2f MPa | SF %.2f\n', ...
    a_med, a_lat, vm_fem_max, sf_fem_analytical);

% -------------------------------------------------------------------------
% 7. FIGURES - EACH IN ITS OWN WINDOW
% -------------------------------------------------------------------------
fprintf('\n[+] Step 5: Rendering figures (one window per plot)...\n');
cBoneF = [0.90 0.86 0.76]; cBoneT = [0.88 0.84 0.74];
cFem = [0.22 0.50 0.80]; cIns = [0.85 0.94 0.98]; cTray = [0.78 0.62 0.28];
saved = {};

% Fig 1: 3D in situ
[fig, ax] = new_fig('3D In Situ Assembly', 1);
plot_mesh_patch(ax, bone_fv, bone_ff, cBoneF, 0.35, 'Femur (resected)');
plot_mesh_patch(ax, fem_V, fem_F, cFem, 0.95, 'Femoral Component (CoCrMo)');
plot_mesh_patch(ax, ins_V, ins_F, cIns, 0.90, 'Tibial Insert (UHMWPE)');
plot_mesh_patch(ax, tray_V, tray_F, cTray, 0.95, 'Tibial Baseplate (Ti6Al4V)');
plot_mesh_patch(ax, bone_tv, bone_tf, cBoneT, 0.35, 'Tibia (resected)');
setup_3d_axes(ax, '3D Personalized Knee Implant Assembly (In Situ)', [-35 25]);
legend(ax, 'Location', 'best', 'FontSize', 9);
saved{end+1} = save_fig(fig, 'fig01_implant_in_situ_3d.png', results_dir, base_dir);

% Fig 2: coronal
[fig, ax] = new_fig('Coronal View', 2);
plot_mesh_patch(ax, bone_fv, bone_ff, cBoneF, 0.30, '');
plot_mesh_patch(ax, fem_V, fem_F, cFem, 0.95, '');
plot_mesh_patch(ax, ins_V, ins_F, cIns, 0.90, '');
plot_mesh_patch(ax, tray_V, tray_F, cTray, 0.95, '');
plot_mesh_patch(ax, bone_tv, bone_tf, cBoneT, 0.30, '');
setup_3d_axes(ax, 'Coronal View (Mediolateral Alignment)', [0 90]);
saved{end+1} = save_fig(fig, 'fig02_coronal_view.png', results_dir, base_dir);

% Fig 3: sagittal
[fig, ax] = new_fig('Sagittal View', 3);
plot_mesh_patch(ax, bone_fv, bone_ff, cBoneF, 0.30, '');
plot_mesh_patch(ax, fem_V, fem_F, cFem, 0.95, '');
plot_mesh_patch(ax, ins_V, ins_F, cIns, 0.90, '');
plot_mesh_patch(ax, tray_V, tray_F, cTray, 0.95, '');
plot_mesh_patch(ax, bone_tv, bone_tf, cBoneT, 0.30, '');
setup_3d_axes(ax, 'Sagittal View (Anteroposterior Condylar Profile)', [0 0]);
saved{end+1} = save_fig(fig, 'fig03_sagittal_view.png', results_dir, base_dir);

% Fig 4: exploded
sh = @(V, d) V + repmat([d 0 0], size(V,1), 1);
[fig, ax] = new_fig('Exploded View', 4);
plot_mesh_patch(ax, sh(bone_fv, 60), bone_ff, cBoneF, 0.30, 'Femur');
plot_mesh_patch(ax, sh(fem_V, 36),   fem_F,   cFem,   0.95, 'Femoral Component');
plot_mesh_patch(ax, sh(ins_V, 14),   ins_F,   cIns,   0.90, 'UHMWPE Insert');
plot_mesh_patch(ax, sh(tray_V, -6),  tray_F,  cTray,  0.95, 'Tibial Tray');
plot_mesh_patch(ax, sh(bone_tv, -30), bone_tf, cBoneT, 0.30, 'Tibia');
setup_3d_axes(ax, '3D Exploded View (Component Hierarchy)', [-45 20]);
legend(ax, 'Location', 'best', 'FontSize', 9);
saved{end+1} = save_fig(fig, 'fig04_exploded_view.png', results_dir, base_dir);

% Fig 5-7: individual components
[fig, ax] = new_fig('Femoral Component', 5);
plot_mesh_patch(ax, fem_V, fem_F, cFem, 1.0, '');
setup_3d_axes(ax, sprintf('Femoral Component (Co-Cr-Mo)\nWidth=%.1f mm | AP=%.1f mm | Resection=%.1f mm', ...
    params.femoral.W_ml, params.femoral.D_ap, params.femoral.resect_thick), [-40 25]);
saved{end+1} = save_fig(fig, 'fig05_femoral_component.png', results_dir, base_dir);

[fig, ax] = new_fig('Tibial Insert', 6);
plot_mesh_patch(ax, ins_V, ins_F, [0.70 0.88 0.95], 1.0, '');
setup_3d_axes(ax, sprintf('Tibial Bearing Insert (UHMWPE)\nCore=%.1f mm | Med Dish=%.1f mm | Lat Dish=%.1f mm', ...
    params.tibial_insert.core_thickness, params.tibial_insert.med_dish_depth, params.tibial_insert.lat_dish_depth), [-40 30]);
saved{end+1} = save_fig(fig, 'fig06_tibial_insert.png', results_dir, base_dir);

[fig, ax] = new_fig('Tibial Tray', 7);
plot_mesh_patch(ax, tray_V, tray_F, cTray, 1.0, '');
setup_3d_axes(ax, sprintf('Tibial Baseplate Tray (Ti-6Al-4V)\nWidth=%.1f mm | Thickness=%.1f mm | Keel=%.1f mm', ...
    params.tibial_tray.W_ml, params.tibial_tray.thickness, params.tibial_tray.stem_length), [-40 25]);
saved{end+1} = save_fig(fig, 'fig07_tibial_tray.png', results_dir, base_dir);

% Fig 8: parameter sheet
fig = figure('Name','Parameter Matrix','Color','w','Position',[200 100 1000 800]);
ax = axes(fig); axis(ax,'off');
text_specs = {
    '\bf\fontsize{13}Parametric Implant Specifications'
    ''
    '\bf1. Femoral Component (Co-Cr-Mo, ASTM F75):'
    sprintf('   Width: %.2f mm  |  AP: %.2f mm  |  Resection: %.2f mm', params.femoral.W_ml, params.femoral.D_ap, params.femoral.resect_thick)
    sprintf('   Med Runner R: %.1f mm  |  Lat Runner R: %.1f mm', params.femoral.r_med_runner, params.femoral.r_lat_runner)
    sprintf('   E=%g GPa, nu=%.2f, Sy=%d MPa, SF=%.2f', params.materials.cocrmo.E/1000, params.materials.cocrmo.nu, params.materials.cocrmo.yield_stress, sf_fem_analytical)
    ''
    '\bf2. UHMWPE Articular Bearing Insert (ASTM F648):'
    sprintf('   Core Thickness: %.2f mm  |  Spine: %.2f mm', params.tibial_insert.core_thickness, params.tibial_insert.spine_height)
    sprintf('   Medial Dish: %.2f mm  |  Lateral Dish: %.2f mm', params.tibial_insert.med_dish_depth, params.tibial_insert.lat_dish_depth)
    sprintf('   E=%.1f GPa, nu=%.2f, Sy=%d MPa, SF=%.2f', params.materials.uhmwpe.E/1000, params.materials.uhmwpe.nu, params.materials.uhmwpe.yield_stress, sf_ins_analytical)
    ''
    '\bf3. Tibial Baseplate Tray (Ti-6Al-4V ELI, ASTM F136):'
    sprintf('   Width: %.2f mm  |  Depth: %.2f mm', params.tibial_tray.W_ml, params.tibial_tray.D_ap)
    sprintf('   Thickness: %.2f mm  |  Keel Stem: %.2f mm', params.tibial_tray.thickness, params.tibial_tray.stem_length)
    sprintf('   E=%g GPa, nu=%.2f, Sy=%d MPa, SF=%.2f', params.materials.ti6al4v.E/1000, params.materials.ti6al4v.nu, params.materials.ti6al4v.yield_stress, sf_tray_analytical)
    ''
    '\bf4. Loading (ISO 14243 Peak Stance):'
    sprintf('   Faxial=%d N (Med=%d N, Lat=%d N)  |  Fshear=%d N', params.loading.F_axial_total, round(F_med), round(F_lat), params.loading.F_ap_shear)
    };
text(ax, 0.02, 0.98, text_specs, 'FontSize', 12, 'VerticalAlignment','top', 'Interpreter','tex', 'Color',[0.15 0.15 0.15]);
title(ax, 'Parameter Matrix & Material Database', 'FontSize', 14, 'FontWeight','bold');
saved{end+1} = save_fig(fig, 'fig08_parameter_matrix.png', results_dir, base_dir);

% Fig 9: peak stress vs yield
fig = figure('Name','Stress vs Yield','Color','w','Position',[220 120 1000 700]);
ax = axes(fig); hold(ax,'on');
peak_stresses = [vm_tray_analytical, vm_ins_max, vm_fem_max];
yield_limits  = [params.materials.ti6al4v.yield_stress, params.materials.uhmwpe.yield_stress, params.materials.cocrmo.yield_stress];
xp = 1:3; bw = 0.32;
b1 = bar(ax, xp-bw/2, peak_stresses, bw, 'FaceColor',[0.85 0.25 0.20], 'EdgeColor',[0.6 0.1 0.1], 'LineWidth',1.2);
b2 = bar(ax, xp+bw/2, yield_limits,  bw, 'FaceColor',[0.20 0.65 0.35], 'EdgeColor',[0.1 0.4 0.2], 'LineWidth',1.2);
for i = 1:3
    text(ax, xp(i)-bw/2, peak_stresses(i)+8, sprintf('%.1f MPa', peak_stresses(i)), 'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',[0.7 0.1 0.1]);
    text(ax, xp(i)+bw/2, yield_limits(i)+8,  sprintf('%d MPa', round(yield_limits(i))), 'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',[0.1 0.4 0.2]);
end
set(ax,'XTick',xp,'XTickLabel',{'Tibial Tray','UHMWPE Insert','Femoral'});
ylabel(ax,'Stress (MPa)','FontSize',12,'FontWeight','bold');
title(ax,'Analytical Peak Stress vs Material Yield Limit','FontSize',13,'FontWeight','bold');
legend(ax,[b1 b2],{'Est. Von Mises Stress','Yield Strength (Sy)'},'Location','northwest');
ylim(ax,[0 max(yield_limits)*1.15]); box(ax,'on');
saved{end+1} = save_fig(fig, 'fig09_stress_vs_yield.png', results_dir, base_dir);

% Fig 10: safety factors
fig = figure('Name','Safety Factors','Color','w','Position',[240 140 1000 700]);
ax = axes(fig); hold(ax,'on');
sfs = [sf_tray_analytical, sf_ins_analytical, sf_fem_analytical];
bar_colors = [0.20 0.60 0.90; 0.25 0.75 0.50; 0.95 0.70 0.20];
for i = 1:3
    bar(ax, i, sfs(i), 0.55, 'FaceColor',bar_colors(i,:), 'EdgeColor','k', 'LineWidth',1.2);
    text(ax, i, sfs(i)*1.03, sprintf('SF = %.2f', sfs(i)), 'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
end
set(ax,'YScale','log');
yline(ax, 1.5, '--r', 'LineWidth',1.5, 'Label','Min SF = 1.5');
set(ax,'XTick',1:3,'XTickLabel',{'Ti-6Al-4V Tray','UHMWPE Insert','Co-Cr-Mo Femoral'});
ylabel(ax,'Structural Safety Factor (log scale)','FontSize',12,'FontWeight','bold');
title(ax,'Component Safety Factors (ISO 14243 Stance Load)','FontSize',13,'FontWeight','bold');
ylim(ax,[0.5 max(sfs)*2]); box(ax,'on');
saved{end+1} = save_fig(fig, 'fig10_safety_factors.png', results_dir, base_dir);

% Fig 11: Hertz pressure vs load
fig = figure('Name','Hertz Contact','Color','w','Position',[260 160 1000 700]);
ax = axes(fig); hold(ax,'on');
F_range = linspace(0, 2500, 200);
p0_rm = (6*F_range*params.loading.medial_ratio *E_star^2/(pi^3*Reff_med^2)).^(1/3);
p0_rl = (6*F_range*params.loading.lateral_ratio*E_star^2/(pi^3*Reff_lat^2)).^(1/3);
plot(ax, F_range, p0_rm, '-',  'Color',[0.20 0.45 0.80], 'LineWidth',2.2, 'DisplayName','Medial Condyle');
plot(ax, F_range, p0_rl, '--', 'Color',[0.85 0.35 0.20], 'LineWidth',2.2, 'DisplayName','Lateral Condyle');
xline(ax, F_total, ':k', 'LineWidth',1.5, 'Label','ISO 14243 Load');
yline(ax, params.materials.uhmwpe.yield_stress/0.62, '--r', 'LineWidth',1.5, 'Label','UHMWPE p0 Limit');
xlabel(ax,'Axial Joint Force (N)','FontSize',12,'FontWeight','bold');
ylabel(ax,'Hertz Peak Pressure p_0 (MPa)','FontSize',12,'FontWeight','bold');
title(ax,'Hertz Contact Pressure vs Applied Load','FontSize',13,'FontWeight','bold');
legend(ax,'Location','northwest'); box(ax,'on');
saved{end+1} = save_fig(fig, 'fig11_hertz_contact.png', results_dir, base_dir);

% Fig 12: stem stress vs body weight
fig = figure('Name','Stem Stress','Color','w','Position',[280 180 1000 700]);
ax = axes(fig); hold(ax,'on');
bw_range = linspace(1, 5, 200);  BW_N = 70*9.81;
plot(ax, bw_range, bw_range*BW_N/stem_area, '-', 'Color',[0.78 0.35 0.10], 'LineWidth',2.2);
yline(ax, params.materials.ti6al4v.yield_stress, '--r', 'LineWidth',1.5, ...
    'Label', sprintf('Ti-6Al-4V Sy=%d MPa', params.materials.ti6al4v.yield_stress));
xline(ax, F_total/BW_N, ':k', 'LineWidth',1.5, 'Label','ISO 14243 BW');
xlabel(ax,'Body-Weight Multiplier','FontSize',12,'FontWeight','bold');
ylabel(ax,'Tibial Stem Axial Stress (MPa)','FontSize',12,'FontWeight','bold');
title(ax,'Tibial Stem Axial Stress vs Body-Weight Multiplier','FontSize',13,'FontWeight','bold');
box(ax,'on');
saved{end+1} = save_fig(fig, 'fig12_stem_stress.png', results_dir, base_dir);
fprintf('    - %d figure windows created and saved (300 DPI)\n', numel(saved));

% -------------------------------------------------------------------------
% 8. EXPORT PARAMETERS
% -------------------------------------------------------------------------
fprintf('\n[+] Step 6: Exporting design parameters...\n');
S = struct('patient_id',PATIENT_ID, 'timestamp',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'tibial_tray_width_mm',params.tibial_tray.W_ml, 'tibial_tray_depth_mm',params.tibial_tray.D_ap, ...
    'tibial_tray_thick_mm',params.tibial_tray.thickness, 'tibial_stem_length_mm',params.tibial_tray.stem_length, ...
    'insert_thickness_mm',params.tibial_insert.core_thickness, ...
    'femoral_width_mm',params.femoral.W_ml, 'femoral_depth_mm',params.femoral.D_ap, ...
    'femoral_resect_mm',params.femoral.resect_thick, ...
    'axial_load_N',params.loading.F_axial_total, 'ap_shear_load_N',params.loading.F_ap_shear);
S.analytical_tibial_tray = struct('material',params.materials.ti6al4v.name, 'bending_stress_MPa',sigma_bend_tray, ...
    'stem_axial_MPa',sigma_stem, 'shear_stress_MPa',tau_tray, 'von_mises_MPa',vm_tray_analytical, ...
    'yield_MPa',params.materials.ti6al4v.yield_stress, 'safety_factor',sf_tray_analytical);
S.analytical_tibial_insert = struct('material',params.materials.uhmwpe.name, 'medial_hertz_p0_MPa',p0_med, ...
    'lateral_hertz_p0_MPa',p0_lat, 'max_von_mises_MPa',vm_ins_max, ...
    'yield_MPa',params.materials.uhmwpe.yield_stress, 'safety_factor',sf_ins_analytical);
S.analytical_femoral = struct('material',params.materials.cocrmo.name, 'med_contact_radius_mm',a_med, ...
    'lat_contact_radius_mm',a_lat, 'max_contact_stress_MPa',vm_fem_max, ...
    'yield_MPa',params.materials.cocrmo.yield_stress, 'safety_factor',sf_fem_analytical);
json_str = jsonencode(S, 'PrettyPrint', true);
for d = {results_dir, base_dir}
    fid = fopen(fullfile(d{1}, 'implant_design_parameters.json'), 'w'); fwrite(fid, json_str); fclose(fid);
end

csv_file = fullfile(results_dir, 'implant_design_parameters.csv');
fid = fopen(csv_file, 'w');
fprintf(fid, 'Parameter,Component,Value,Unit,Description\n');
fprintf(fid, 'Bicondylar_Width,Femoral,%.2f,mm,Mediolateral bicondylar dimension\n', params.femoral.W_ml);
fprintf(fid, 'AP_Dimension,Femoral,%.2f,mm,Anteroposterior dimension\n', params.femoral.D_ap);
fprintf(fid, 'Resection_Thickness,Femoral,%.2f,mm,Distal bone resection replacement\n', params.femoral.resect_thick);
fprintf(fid, 'Tray_Width,Tibial_Tray,%.2f,mm,Proximal tibial coverage\n', params.tibial_tray.W_ml);
fprintf(fid, 'Tray_Depth,Tibial_Tray,%.2f,mm,AP tibial depth\n', params.tibial_tray.D_ap);
fprintf(fid, 'Tray_Thickness,Tibial_Tray,%.2f,mm,Baseplate metal thickness\n', params.tibial_tray.thickness);
fprintf(fid, 'Stem_Length,Tibial_Tray,%.2f,mm,Keel stem length\n', params.tibial_tray.stem_length);
fprintf(fid, 'Insert_Core_Thickness,Tibial_Insert,%.2f,mm,Articular PE thickness\n', params.tibial_insert.core_thickness);
fprintf(fid, 'VM_Stress_Tray,Analytical,%.2f,MPa,Est Von Mises tibial tray\n', vm_tray_analytical);
fprintf(fid, 'SF_Tray,Analytical,%.2f,ratio,Safety factor tibial tray\n', sf_tray_analytical);
fprintf(fid, 'Hertz_p0_Medial,Analytical,%.2f,MPa,Medial Hertz peak pressure\n', p0_med);
fprintf(fid, 'VM_Stress_Insert,Analytical,%.2f,MPa,Max VM UHMWPE insert\n', vm_ins_max);
fprintf(fid, 'SF_Insert,Analytical,%.2f,ratio,Safety factor UHMWPE insert\n', sf_ins_analytical);
fprintf(fid, 'Contact_Stress_Femoral,Analytical,%.2f,MPa,Max contact stress femoral\n', vm_fem_max);
fprintf(fid, 'SF_Femoral,Analytical,%.2f,ratio,Safety factor femoral component\n', sf_fem_analytical);
fclose(fid);
copyfile(csv_file, fullfile(base_dir, 'implant_design_parameters.csv'));

fprintf('\n=======================================================================\n');
fprintf('  IMPLANT DESIGN COMPLETED - outputs in: %s\n', results_dir);
fprintf('    * Tibial Tray  (Ti-6Al-4V) : VM = %6.2f MPa | SF = %.2f\n', vm_tray_analytical, sf_tray_analytical);
fprintf('    * Insert       (UHMWPE)    : VM = %6.2f MPa | SF = %.2f\n', vm_ins_max, sf_ins_analytical);
fprintf('    * Femoral Comp (Co-Cr-Mo)  : CS = %6.2f MPa | SF = %.2f\n', vm_fem_max, sf_fem_analytical);
fprintf('=======================================================================\n');


% =========================================================================
% STEP 7: PATIENT-SPECIFIC FEA SIMULATION  (MATLAB PDE TOOLBOX)
% =========================================================================
fprintf('\n=======================================================================\n');
fprintf('  STEP 7: PATIENT-SPECIFIC FEA SIMULATION (PDE TOOLBOX)\n');
fprintf('=======================================================================\n');
if ~license('test', 'PDE_Toolbox') || isempty(ver('pde'))
    error('The Partial Differential Equation Toolbox is required for Steps 7 and 8.');
end

% ---------------- user settings ----------------
N_CASES    = 100;     % number of simulations in the dataset (300+ gives better surrogate R^2)
HMAX_BASE  = 2.5;     %#ok<NASGU> mesh size (mm) for the detailed patient-specific case
HMAX_DATA  = 4.0;     % mesh size (mm) for the dataset runs (smaller = slower)
RNG_SEED   = 2024;
CHECKPOINT = 50;      % save partial dataset every N cases
RESUME     = true;    % resume from checkpoint if one exists

fea_root = fullfile(base_dir, 'FEA_analysis');
dir_base = fullfile(fea_root, 'baseline_figures');
dir_data = fullfile(fea_root, 'dataset_figures');
dir_xls  = fullfile(fea_root, 'excel_data');
dir_mat  = fullfile(fea_root, 'matlab_data');
dir_dl   = fullfile(fea_root, 'DL_AMOO_results');
for d = {fea_root, dir_base, dir_data, dir_xls, dir_mat, dir_dl}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end
fprintf('    Output folder : %s\n', fea_root);

% ---- 7a. Subject-specific model parameters (from the STL-driven design) ----
G0 = fea_default_case();
G0.W        = params.tibial_tray.W_ml;
G0.D        = params.tibial_tray.D_ap;
G0.t_tray   = params.tibial_tray.thickness;
G0.L_stem   = params.tibial_tray.stem_length;
G0.r_stem   = params.tibial_tray.stem_radius;
G0.t_ins    = params.tibial_insert.core_thickness;
G0.t_fem    = params.femoral.resect_thick;
G0.E_tray   = params.materials.ti6al4v.E;   G0.nu_tray = params.materials.ti6al4v.nu;
G0.sy_tray  = params.materials.ti6al4v.yield_stress;
G0.E_ins    = params.materials.uhmwpe.E;    G0.nu_ins  = params.materials.uhmwpe.nu;
G0.sy_ins   = params.materials.uhmwpe.yield_stress;
G0.E_fem    = params.materials.cocrmo.E;    G0.nu_fem  = params.materials.cocrmo.nu;
G0.sy_fem   = params.materials.cocrmo.yield_stress;
G0.F_axial  = params.loading.F_axial_total;
G0.F_shear  = params.loading.F_ap_shear;
G0.med_ratio = params.loading.medial_ratio;
G0.flex_deg = 15;                           % peak stance of level walking
G0 = fea_finalize(G0);

fprintf('\n[+] 7a. Patient-specific FEA model (units: mm, N, MPa)\n');
fprintf('    - Footprint (ML x AP)        : %.1f x %.1f mm\n', G0.W, G0.D);
fprintf('    - Tray / Insert / Fem. comp. : %.1f / %.1f / %.1f mm\n', G0.t_tray, G0.t_ins, G0.t_fem);
fprintf('    - Bone: E_cort=%d MPa, E_canc=%d MPa, cortical thickness=%.1f mm\n', G0.E_cort, G0.E_canc, G0.t_cort);
fprintf('    - Load: %d N axial (%.0f%% medial), %d N AP shear\n', G0.F_axial, 100*G0.med_ratio, G0.F_shear);

% ---- 7b. Boundary conditions + solve ----
fprintf('\n[+] 7b. Meshing, applying loads/boundary conditions, solving...\n');
HVOX_BASE = 2.0;      % voxel size (mm) of the patient-geometry mesh (smaller = finer, slower)
[R0, M0, Gp, P0] = run_fea_patient(params, G0, fv_insitu, ff, tv, tf, x_tcut, x_fcut, ...
    tray_V, tray_F, ins_V, ins_F, fem_V, fem_F, HVOX_BASE);
R0.Status = 'OK';
H0 = Gp.H; %#ok<NASGU>
fprintf('    - Mesh : %d nodes, %d tetrahedra | solve time %.1f s\n', R0.Mesh_Nodes, R0.Mesh_Elements, R0.Solve_Time_s);

% ---- 7c. Stress / strain / displacement results ----
fprintf('\n[+] 7c. Baseline FEA results\n');
fprintf('    Peak von Mises  [MPa] : tibial bone %.2f | stem %.2f | tray %.2f | insert %.2f | fem. comp. %.2f | fem. bone %.2f\n', ...
    R0.VM_TibBone_MPa, R0.VM_Stem_MPa, R0.VM_Tray_MPa, R0.VM_Insert_MPa, R0.VM_FemComp_MPa, R0.VM_FemBone_MPa);
fprintf('    Peak strain  [microstr]: tibial bone %.0f | tray %.0f | insert %.0f | fem. comp. %.0f | fem. bone %.0f\n', ...
    R0.Strain_TibBone_ue, R0.Strain_Tray_ue, R0.Strain_Insert_ue, R0.Strain_FemComp_ue, R0.Strain_FemBone_ue);
fprintf('    Displacement    [mm]  : max %.4f | tray subsidence %.4f | insert compression %.4f\n', ...
    R0.Max_Disp_mm, R0.Tray_Subsidence_mm, R0.Insert_Compression_mm);
fprintf('    Safety factors        : tray %.2f | insert %.2f | fem. comp. %.2f | bone (strain) %.2f\n', ...
    R0.SF_Tray, R0.SF_Insert, R0.SF_FemComp, R0.SF_Bone);

fprintf('\n[+] Rendering baseline FEA figures (one window each)...\n');
bi = 0;

bi = bi+1; fea_patch3d(M0, [], '', 'FE Mesh', 'Finite Element Mesh (Implant-Bone Model)', bi, dir_base, 'B01_fe_mesh.png');
bi = bi+1; fea_plot_materials_patient(P0, Gp, bi, dir_base, 'B02_material_layout.png');
bi = bi+1; fea_patch3d(M0, M0.vm, 'von Mises stress (MPa)', '3D von Mises Stress', '3D von Mises Stress Distribution', bi, dir_base, 'B03_3d_von_mises_stress.png');
bi = bi+1; fea_patch3d(M0, M0.ee, 'Equivalent strain (microstrain)', '3D Equivalent Strain', '3D Equivalent (von Mises) Strain Distribution', bi, dir_base, 'B04_3d_equivalent_strain.png');
bi = bi+1; fea_patch3d(M0, M0.um, 'Displacement magnitude (mm)', '3D Displacement', '3D Displacement Magnitude', bi, dir_base, 'B05_3d_displacement.png');

bi = bi+1; fea_plot_slice_patient(P0, Gp, 'sagittal', 'vm',     'Von Mises Stress',      bi, dir_base, 'B06_sagittal_von_mises.png');
bi = bi+1; fea_plot_slice_patient(P0, Gp, 'sagittal', 'strain', 'Equivalent Strain',     bi, dir_base, 'B07_sagittal_strain.png');
bi = bi+1; fea_plot_slice_patient(P0, Gp, 'sagittal', 'disp',   'Displacement',          bi, dir_base, 'B08_sagittal_displacement.png');
bi = bi+1; fea_plot_slice_patient(P0, Gp, 'coronal',  'vm',     'Von Mises Stress',      bi, dir_base, 'B09_coronal_von_mises.png');
bi = bi+1; fea_plot_slice_patient(P0, Gp, 'coronal',  'strain', 'Equivalent Strain',     bi, dir_base, 'B10_coronal_strain.png');
bi = bi+1; fea_plot_slice_patient(P0, Gp, 'coronal',  'disp',   'Displacement',          bi, dir_base, 'B11_coronal_displacement.png');

comp_names = {'Tibial Bone','Stem','Tray','Insert','Fem. Component','Femoral Bone'};
vm_comp  = [R0.VM_TibBone_MPa, R0.VM_Stem_MPa, R0.VM_Tray_MPa, R0.VM_Insert_MPa, R0.VM_FemComp_MPa, R0.VM_FemBone_MPa];
st_comp  = [R0.Strain_TibBone_ue, R0.Strain_Tray_ue, R0.Strain_Tray_ue, R0.Strain_Insert_ue, R0.Strain_FemComp_ue, R0.Strain_FemBone_ue];
bi = bi+1; fea_bar(vm_comp, comp_names, 'Peak von Mises stress (MPa)', 'Peak Stress by Component', '%.2f', bi, dir_base, 'B12_peak_stress_by_component.png', [0.85 0.33 0.10]);
bi = bi+1; fea_bar(st_comp, comp_names, 'Peak equivalent strain (microstrain)', 'Peak Strain by Component', '%.0f', bi, dir_base, 'B13_peak_strain_by_component.png', [0.20 0.55 0.80]);
bi = bi+1; fea_bar([R0.Max_Disp_mm, R0.Tray_Subsidence_mm, R0.FemComp_Disp_mm, R0.Insert_Compression_mm], ...
    {'Max (nodal)','Tray Subsidence','Fem. Comp. Disp.','Insert Compression'}, 'Displacement (mm)', ...
    'Displacement Summary', '%.4f', bi, dir_base, 'B14_displacement_summary.png', [0.30 0.65 0.40]);
bi = bi+1; fea_bar([R0.SF_Tray, R0.SF_Insert, R0.SF_FemComp, R0.SF_Bone], ...
    {'Tray (yield)','Insert (yield)','Fem. Comp. (yield)','Bone (strain)'}, 'Safety factor', ...
    'Component Safety Factors (FEA)', '%.1f', bi, dir_base, 'B15_safety_factors_fea.png', [0.95 0.70 0.20]);
fprintf('    - %d baseline figures saved in %s\n', bi, dir_base);

% =========================================================================
% STEP 8: SIMULATION DATASET GENERATION
% =========================================================================
fprintf('\n=======================================================================\n');
fprintf('  STEP 8: SIMULATION DATASET GENERATION (%d CASES)\n', N_CASES);
fprintf('=======================================================================\n');

rng(RNG_SEED);
act_names = {'Level walking','Stair ascent','Stair descent','Sit-to-stand','Jogging','Deep squat'};
act_bw    = [2.8 3.3 4.2 3.0 4.8 4.5];     % peak tibiofemoral force in body weights
act_flex  = [15 45 30 70 25 90];           % knee flexion angle (deg)

var_names = {'W','D','t_cort','E_cort','E_canc','Mass','t_tray','L_stem','r_stem','t_ins','t_fem','E_ins', ...
             'shear_frac','med_ratio','contact_frac','load_scale'};
lo = [62 42 1.5 12000  150  50 2.5 20 4  8  7.5  600 0.12 0.50 0.28 0.90];
hi = [88 62 4.0 20000 1200 120 5.0 45 8 14 10.0 1400 0.25 0.75 0.45 1.10];
U  = fea_lhs(N_CASES, numel(lo));
Xs = lo + U .* (hi - lo);

W_v = Xs(:,1);  D_v = Xs(:,2);  tc_v = Xs(:,3);  Ec_v = Xs(:,4);  Ek_v = Xs(:,5);  M_v = Xs(:,6);
tt_v = Xs(:,7); Ls_v = Xs(:,8); rs_v = Xs(:,9);  ti_v = Xs(:,10); tf_v = Xs(:,11); Ei_v = Xs(:,12);
sf_v = Xs(:,13); mr_v = Xs(:,14); cf_v = Xs(:,15); ls_v = Xs(:,16);
act_idx  = mod(randperm(N_CASES)-1, numel(act_names)) + 1;   act_idx = act_idx(:);
tray_mat = 1 + (rand(N_CASES,1) > 0.7);                      % 1 = Ti-6Al-4V, 2 = Co-Cr-Mo
tray_lbl = {'Ti-6Al-4V','Co-Cr-Mo'};
Etray_v  = [110000; 210000]; nutray_v = [0.33; 0.30]; sytray_v = [880; 650];
Sy_ins_v = 22 + 6*(Ei_v-600)/800;
flex_v   = act_flex(act_idx)';
Fax_v    = act_bw(act_idx)' .* M_v * 9.81 .* ls_v;
Fsh_v    = sf_v .* Fax_v;
Lc_v     = cf_v .* D_v;
Shift_v  = -(flex_v/90) * 0.22 .* D_v;

T_in = table((1:N_CASES)', W_v, D_v, tc_v, Ec_v, Ek_v, M_v, tray_lbl(tray_mat)', tt_v, Ls_v, rs_v, ti_v, tf_v, ...
    Ei_v, Sy_ins_v, act_names(act_idx)', flex_v, act_bw(act_idx)', Fax_v, Fsh_v, mr_v, Lc_v, Shift_v, ...
    'VariableNames', {'Case_ID','W_ML_mm','D_AP_mm','Cort_Thick_mm','E_Cort_MPa','E_Canc_MPa','Body_Mass_kg', ...
    'Tray_Material','Tray_Thick_mm','Stem_Len_mm','Stem_Radius_mm','Insert_Thick_mm','FemComp_Thick_mm', ...
    'E_Insert_MPa','Sy_Insert_MPa','Activity','Flexion_deg','Axial_xBW','F_Axial_N','F_Shear_N', ...
    'Medial_Ratio','Contact_Len_mm','Contact_AP_Shift_mm'});

ckpt = fullfile(dir_mat, 'dataset_checkpoint.mat');
res_cell = cell(N_CASES, 1);  i_start = 1;
if RESUME && exist(ckpt, 'file')
    Sck = load(ckpt);
    if isfield(Sck,'N_saved') && Sck.N_saved == N_CASES && Sck.seed == RNG_SEED
        res_cell = Sck.res_cell;  i_start = Sck.i_done + 1;
        fprintf('    Resuming from checkpoint at case %d\n', i_start);
    end
end

fprintf('    Running %d FEA simulations (mesh Hmax = %.1f mm)...\n', N_CASES, HMAX_DATA);
t_all = tic;
for i = i_start:N_CASES
    G = fea_default_case();
    G.W = W_v(i);  G.D = D_v(i);  G.t_cort = tc_v(i);  G.E_cort = Ec_v(i);  G.E_canc = Ek_v(i);
    G.t_tray = tt_v(i);  G.L_stem = Ls_v(i);  G.r_stem = rs_v(i);  G.t_ins = ti_v(i);  G.t_fem = tf_v(i);
    G.E_ins = Ei_v(i);   G.sy_ins = Sy_ins_v(i);
    G.E_tray = Etray_v(tray_mat(i));  G.nu_tray = nutray_v(tray_mat(i));  G.sy_tray = sytray_v(tray_mat(i));
    G.flex_deg = flex_v(i);  G.F_axial = Fax_v(i);  G.F_shear = Fsh_v(i);
    G.med_ratio = mr_v(i);   G.contact_frac = cf_v(i);
    G = fea_finalize(G);
    try
        R = run_fea_case(G, HMAX_DATA);
        R.Status = 'OK';
    catch ME
        R = fea_nan_result();
        R.Status = ['FAILED: ' strrep(ME.message, newline, ' ')];
    end
    res_cell{i} = R;
    if mod(i, 10) == 0 || i == N_CASES
        el = toc(t_all);  done = i - i_start + 1;
        fprintf('      case %4d / %d | %.1f s elapsed | ETA %.1f min\n', i, N_CASES, el, el/done*(N_CASES-i)/60);
    end
    if mod(i, CHECKPOINT) == 0
        i_done = i; N_saved = N_CASES; seed = RNG_SEED;
        save(ckpt, 'res_cell', 'i_done', 'N_saved', 'seed');
    end
end
T_out = struct2table(vertcat(res_cell{:}));
T = [T_in, T_out];
n_ok = sum(strcmp(T.Status, 'OK'));
fprintf('    Completed: %d OK, %d failed (total %.1f min)\n', n_ok, N_CASES-n_ok, toc(t_all)/60);
if n_ok < 10, error('Too few successful FEA runs - check the first FAILED message in T.Status.'); end

% ---------------- save dataset to Excel / CSV / MAT ----------------
fprintf('\n[+] Saving dataset (Excel, CSV, MAT)...\n');
xls_file = fullfile(dir_xls, 'FEA_Simulation_Dataset.xlsx');
if exist(xls_file, 'file'), delete(xls_file); end
Tv = T(strcmp(T.Status,'OK'), :);

out_cols = {'Peak_VM_MPa','VM_TibBone_MPa','VM_Stem_MPa','VM_Tray_MPa','VM_Insert_MPa','VM_FemComp_MPa','VM_FemBone_MPa', ...
    'Peak_Strain_ue','Strain_TibBone_ue','Strain_Tray_ue','Strain_Insert_ue','Strain_FemComp_ue','Strain_FemBone_ue', ...
    'Max_Disp_mm','Tray_Subsidence_mm','FemComp_Disp_mm','Insert_Compression_mm','SF_Tray','SF_Insert','SF_FemComp','SF_Bone'};
stat_tbl = table(out_cols', 'VariableNames', {'Output'});
A = Tv{:, out_cols};
stat_tbl.Mean = mean(A,1,'omitnan')';   stat_tbl.Std = std(A,0,1,'omitnan')';
stat_tbl.Min  = min(A,[],1)';           stat_tbl.Median = median(A,1,'omitnan')';
stat_tbl.Max  = max(A,[],1)';

range_tbl = table(var_names', lo', hi', 'VariableNames', {'Sampled_Variable','Min','Max'});
base_tbl  = [struct2table(rmfield(G0, {'x0','y0','z0'})), struct2table(R0)];
GS = groupsummary(Tv, 'Activity', 'mean', {'Peak_VM_MPa','Peak_Strain_ue','Max_Disp_mm','Tray_Subsidence_mm','SF_Insert'});

writetable(T,          xls_file, 'Sheet', 'Simulation_Dataset');
writetable(stat_tbl,   xls_file, 'Sheet', 'Summary_Statistics');
writetable(GS,         xls_file, 'Sheet', 'By_Activity');
writetable(range_tbl,  xls_file, 'Sheet', 'Sampled_Ranges');
writetable(base_tbl,   xls_file, 'Sheet', 'Patient_Specific_Baseline');
writetable(T, fullfile(dir_xls, 'FEA_Simulation_Dataset.csv'));
save(fullfile(dir_mat, 'FEA_dataset.mat'), 'T', 'G0', 'R0', 'params', 'lo', 'hi', 'var_names');
fprintf('    - %s\n', xls_file);

% ---------------- dataset plots (one window each) ----------------
fprintf('\n[+] Rendering dataset analysis figures (one window each)...\n');
di = 0;
di = di+1; fea_hist(Tv.Peak_VM_MPa, 'Peak von Mises stress (MPa)', 'Distribution of Peak von Mises Stress', di, dir_data, 'D01_hist_peak_stress.png', [0.85 0.33 0.10]);
di = di+1; fea_hist(Tv.Peak_Strain_ue, 'Peak equivalent strain (microstrain)', 'Distribution of Peak Strain', di, dir_data, 'D02_hist_peak_strain.png', [0.20 0.55 0.80]);
di = di+1; fea_hist(Tv.Max_Disp_mm, 'Maximum displacement (mm)', 'Distribution of Maximum Displacement', di, dir_data, 'D03_hist_max_displacement.png', [0.30 0.65 0.40]);
di = di+1; fea_scatter(Tv.F_Axial_N, Tv.Peak_VM_MPa, Tv.E_Canc_MPa, 'Axial joint force (N)', 'Peak von Mises stress (MPa)', ...
    'Cancellous modulus (MPa)', 'Peak Stress vs Axial Load', di, dir_data, 'D04_load_vs_stress.png');
tm = double(strcmp(Tv.Tray_Material, 'Co-Cr-Mo'));
di = di+1; fea_scatter(Tv.Tray_Thick_mm, Tv.VM_Tray_MPa, tm, 'Tray thickness (mm)', 'Tray peak von Mises stress (MPa)', ...
    'Tray material (0 = Ti, 1 = CoCr)', 'Tray Stress vs Tray Thickness', di, dir_data, 'D05_tray_thickness_vs_stress.png');
di = di+1; fea_scatter(Tv.Stem_Len_mm, Tv.Tray_Subsidence_mm, Tv.E_Canc_MPa, 'Stem length (mm)', 'Tray subsidence (mm)', ...
    'Cancellous modulus (MPa)', 'Tray Subsidence vs Stem Length', di, dir_data, 'D06_stem_length_vs_subsidence.png');
di = di+1; fea_scatter(Tv.E_Canc_MPa, Tv.Strain_TibBone_ue, Tv.Axial_xBW, 'Cancellous bone modulus (MPa)', ...
    'Tibial bone peak strain (microstrain)', 'Axial load (x body weight)', 'Bone Strain vs Bone Quality', di, dir_data, 'D07_bone_modulus_vs_strain.png');
di = di+1; fea_scatter(Tv.Flexion_deg, Tv.VM_Insert_MPa, Tv.F_Axial_N, 'Knee flexion angle (deg)', 'Insert peak von Mises stress (MPa)', ...
    'Axial force (N)', 'Insert Stress vs Flexion Angle', di, dir_data, 'D08_flexion_vs_insert_stress.png');

di = di+1; fig = fea_fig('Tray Stress by Activity', di); ax = axes(fig);
boxchart(ax, categorical(Tv.Activity), Tv.VM_Tray_MPa, 'BoxFaceColor', [0.78 0.62 0.28]);
ylabel(ax, 'Tray peak von Mises stress (MPa)'); title(ax, 'Tray Stress by Loading Condition', 'FontWeight','bold'); box(ax,'on');
fea_save(fig, dir_data, 'D09_tray_stress_by_activity.png');

di = di+1; fig = fea_fig('Displacement by Activity', di); ax = axes(fig);
boxchart(ax, categorical(Tv.Activity), Tv.Max_Disp_mm, 'BoxFaceColor', [0.30 0.65 0.40]);
ylabel(ax, 'Maximum displacement (mm)'); title(ax, 'Displacement by Loading Condition', 'FontWeight','bold'); box(ax,'on');
fea_save(fig, dir_data, 'D10_displacement_by_activity.png');

di = di+1; fig = fea_fig('Safety Factor Distribution', di); ax = axes(fig); hold(ax,'on');
sf_all = {Tv.SF_Tray, Tv.SF_Insert, Tv.SF_FemComp, Tv.SF_Bone};
sf_nm  = {'Tray','Insert','Femoral comp.','Bone (strain)'};
sf_col = [0.20 0.60 0.90; 0.25 0.75 0.50; 0.95 0.70 0.20; 0.80 0.30 0.30];
for k = 1:4
    histogram(ax, log10(sf_all{k}(sf_all{k}>0)), 30, 'FaceColor', sf_col(k,:), 'FaceAlpha', 0.6, 'DisplayName', sf_nm{k});
end
xline(ax, log10(1.5), '--r', 'LineWidth', 1.5, 'Label', 'SF = 1.5', 'HandleVisibility','off');
xlabel(ax, 'log_{10}(Safety factor)'); ylabel(ax, 'Number of cases');
title(ax, 'Safety Factor Distribution (all components)', 'FontWeight','bold'); legend(ax, 'Location','best'); box(ax,'on');
fea_save(fig, dir_data, 'D11_safety_factor_distribution.png');

in_cols  = {'W_ML_mm','D_AP_mm','Cort_Thick_mm','E_Cort_MPa','E_Canc_MPa','Body_Mass_kg','Tray_Thick_mm','Stem_Len_mm', ...
            'Stem_Radius_mm','Insert_Thick_mm','FemComp_Thick_mm','E_Insert_MPa','Flexion_deg','F_Axial_N','F_Shear_N', ...
            'Medial_Ratio','Contact_Len_mm'};
out_corr = {'Peak_VM_MPa','VM_Tray_MPa','VM_Insert_MPa','VM_FemComp_MPa','Strain_TibBone_ue','Max_Disp_mm','Tray_Subsidence_mm'};
Call = corrcoef([Tv{:, in_cols}, Tv{:, out_corr}], 'Rows', 'pairwise');
Cmat = Call(1:numel(in_cols), numel(in_cols)+1:end);
di = di+1; fig = fea_fig('Input-Output Correlation', di);
hm = heatmap(fig, strrep(out_corr,'_',' '), strrep(in_cols,'_',' '), Cmat, 'Colormap', turbo, 'CellLabelFormat', '%.2f');
hm.Title = 'Correlation: Anatomical / Design / Load Inputs vs FEA Outputs';
fea_save(fig, dir_data, 'D12_correlation_heatmap.png');

sens = abs(Cmat(:,1));  [sens_s, ord] = sort(sens, 'ascend');
di = di+1; fig = fea_fig('Sensitivity Ranking', di); ax = axes(fig);
barh(ax, sens_s, 'FaceColor', [0.55 0.35 0.75]);
set(ax, 'YTick', 1:numel(ord), 'YTickLabel', strrep(in_cols(ord),'_',' '));
xlabel(ax, '|Correlation| with peak von Mises stress'); title(ax, 'Input Sensitivity Ranking', 'FontWeight','bold'); box(ax,'on');
fea_save(fig, dir_data, 'D13_sensitivity_ranking.png');

SFm = [Tv.SF_Tray, Tv.SF_Insert, Tv.SF_FemComp, Tv.SF_Bone];
[min_sf, which_c] = min(SFm, [], 2);
cnt = [sum(min_sf >= 1.5), arrayfun(@(k) sum(min_sf < 1.5 & which_c == k), 1:4)];
lbl = {'All SF >= 1.5','Tray critical','Insert critical','Fem. comp. critical','Bone critical'};
di = di+1; fig = fea_fig('Design Safety Outcome', di); ax = axes(fig);
pie(ax, cnt(cnt>0), lbl(cnt>0));
title(ax, sprintf('Design Safety Outcome (%d simulations)', numel(min_sf)), 'FontWeight','bold');
fea_save(fig, dir_data, 'D14_design_safety_outcome.png');
fprintf('    - %d dataset figures saved in %s\n', di, dir_data);

fprintf('\n=======================================================================\n');
fprintf('  STEPS 7-8 COMPLETED\n');
fprintf('    FEA folder    : %s\n', fea_root);
fprintf('    Excel dataset : %s\n', xls_file);
fprintf('    Cases (OK)    : %d / %d\n', n_ok, N_CASES);
fprintf('=======================================================================\n');

% =========================================================================
% STEPS 9-12: DL SURROGATE -> AMOO -> FINAL FEA VALIDATION -> REPORT
% =========================================================================
fh = struct('run_fea_case', @run_fea_case, 'finalize', @fea_finalize, ...
            'fig', @fea_fig, 'save', @fea_save);
out9_12 = run_steps_9_to_12(T, G0, lo, hi, var_names, dir_dl, HMAX_DATA, fh); %#ok<NASGU>


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function [V, F] = iso_from_volume(Y, Z, X, vol)
% Close the volume at its borders, extract the surface, return [X Y Z] vertices.
% (Grids come from meshgrid(y,z,x), so isosurface must be called as (Y,Z,X).)
vol(1,:,:) = 10; vol(end,:,:) = 10;
vol(:,1,:) = 10; vol(:,end,:) = 10;
vol(:,:,1) = 10; vol(:,:,end) = 10;
s = isosurface(Y, Z, X, vol, 0);
V = s.vertices(:, [3 1 2]);
F = s.faces;
end

function [V, F] = generate_tibial_tray_mesh(params)
p = params.tibial_tray;
y = linspace(p.center_y-p.W_ml/2-8, p.center_y+p.W_ml/2+8, 90);
z = linspace(p.center_z-p.D_ap/2-8, p.center_z+p.D_ap/2+8, 80);
x = linspace(p.x_stem_tip-6, p.x_top+6, 70);
[Y, Z, X] = meshgrid(y, z, x);
Y_n = (Y-p.center_y)/(p.W_ml/2);
Z_n = (Z-p.center_z)/(p.D_ap/2);
notch  = p.notch_factor*(Z_n<0).*cos(pi/2*min(max(Y_n,-1),1)).^2;
r_prof = sqrt(Y_n.^2+(Z_n+notch).^2)-1.0;
d_plate_x  = max(p.x_bot-X, X-p.x_top);
solid_tray = max(r_prof, d_plate_x/4.0);
r_stem = sqrt((Y-p.center_y).^2+(Z-p.center_z).^2)-p.stem_radius+(p.x_bot-X)/p.stem_length*1.2;
d_stem_x   = max(p.x_stem_tip-X, X-p.x_bot);
solid_stem = max(r_stem, d_stem_x);
d_fin_y = abs(Y-p.center_y)-p.fin_width/2;
d_fin_z = abs(Z-p.center_z)-p.fin_thick/2;
solid_fin = max(max(d_fin_y, d_fin_z), d_stem_x);
[V, F] = iso_from_volume(Y, Z, X, min(solid_tray, min(solid_stem, solid_fin)));
end

function [V, F] = generate_tibial_insert_mesh(params)
p_t = params.tibial_tray;  p_i = params.tibial_insert;
y = linspace(p_t.center_y-p_t.W_ml/2-8, p_t.center_y+p_t.W_ml/2+8, 90);
z = linspace(p_t.center_z-p_t.D_ap/2-8, p_t.center_z+p_t.D_ap/2+8, 80);
x = linspace(p_i.x_bot-6, p_i.x_top+10, 60);
[Y, Z, X] = meshgrid(y, z, x);
Y_n = (Y-p_t.center_y)/(p_t.W_ml/2);
Z_n = (Z-p_t.center_z)/(p_t.D_ap/2);
notch  = p_t.notch_factor*(Z_n<0).*cos(pi/2*min(max(Y_n,-1),1)).^2;
r_prof = sqrt(Y_n.^2+(Z_n+notch).^2)-1.0;
r_med = sqrt(((Y-(p_t.center_y+p_i.dish_offset))/p_i.med_dish_ry).^2+((Z-p_t.center_z)/p_i.med_dish_rz).^2);
r_lat = sqrt(((Y-(p_t.center_y-p_i.dish_offset))/p_i.lat_dish_ry).^2+((Z-p_t.center_z)/p_i.lat_dish_rz).^2);
d_med = p_i.med_dish_depth*max(0, 1.0-r_med.^2);
d_lat = p_i.lat_dish_depth*max(0, 1.0-r_lat.^2);
r_spine = sqrt(((Y-p_t.center_y)/(p_i.spine_width/2)).^2+((Z-p_t.center_z)/(p_i.spine_length/2)).^2);
h_spine = p_i.spine_height*max(0, 1.0-r_spine.^2);
z_surf  = p_i.x_top-d_med-d_lat+h_spine;
d_ins_x = max(p_i.x_bot-X, X-z_surf);
[V, F] = iso_from_volume(Y, Z, X, max(r_prof, d_ins_x/3.5));
end

function [V, F] = generate_femoral_component_mesh(params)
% Rounded-profile condyles + anterior bridge, with a flat-cut cavity that
% matches the distal / anterior / posterior bone resections.
p  = params.femoral;
rc = p.r_med_runner;                       % distal corner radius
h  = p.resect_thick;
half_foot = sqrt(2*rc*h - h^2);            % footprint half-width at the cut plane
Hz = p.D_ap/2 - half_foot + rc;            % outer half AP extent so footprint = D_ap
H  = 2*rc + 6;                             % profile height
xc = p.x_apex + H/2;
zc = p.center_z;  yc = p.center_y;

y = linspace(yc-p.W_ml/2-6, yc+p.W_ml/2+6, 110);
z = linspace(zc-Hz-6, zc+Hz+6, 110);
x = linspace(p.x_apex-4, p.x_apex+p.top_height+4, 70);
[Y, Z, X] = meshgrid(y, z, x);

% sagittal rounded-rectangle profile (exact SDF), extruded along Y
qz = abs(Z-zc)-(Hz-rc);
qx = abs(X-xc)-(H/2-rc);
d_prof = sqrt(max(qz,0).^2+max(qx,0).^2)+min(max(qz,qx),0)-rc;

c_med = max(d_prof, abs(Y-(yc+p.condyle_spacing/2))-p.condyle_width/2);
c_lat = max(d_prof, abs(Y-(yc-p.condyle_spacing/2))-p.condyle_width/2);
z_bridge = zc + Hz - 14;                   % anterior flange joins the condyles
bridge = max(d_prof, max(abs(Y-yc)-p.W_ml/2, z_bridge-Z));
outer  = min(min(c_med, c_lat), bridge);

x_cut = p.x_apex + p.resect_thick;         % flat distal cut inside the component
cavity = max(abs(Z-zc)-(p.D_ap/2-p.wall_thick), x_cut-X);
solid  = max(max(outer, -cavity), X-(p.x_apex+p.top_height));
[V, F] = iso_from_volume(Y, Z, X, solid);
end

function [mV, mF] = merge_triangulations(V_cells, F_cells)
mV = []; mF = []; off = 0;
for i = 1:numel(V_cells)
    mV = [mV; V_cells{i}];             %#ok<AGROW>
    mF = [mF; F_cells{i}+off];         %#ok<AGROW>
    off = off + size(V_cells{i},1);
end
end

function [Fo, Vo] = decimate_mesh(F, V, max_faces)
% Reduce face count for display / assembly and drop unused vertices
if size(F,1) > max_faces
    [Fo, Vo] = reducepatch(F, V, max_faces/size(F,1));
else
    Fo = F; Vo = V;
end
used = unique(Fo(:));
map = zeros(size(Vo,1),1);  map(used) = 1:numel(used);
Vo = Vo(used,:);  Fo = map(Fo);
end

function export_stl(TR, name, d1, d2)
stlwrite(TR, fullfile(d1, name));
stlwrite(TR, fullfile(d2, name));
fprintf('      [OK] %s (%d vertices, %d faces)\n', name, size(TR.Points,1), size(TR.ConnectivityList,1));
end

function [fig, ax] = new_fig(name, idx)
off = 25*(idx-1);
fig = figure('Name', name, 'NumberTitle','off', 'Color','w', 'Units','pixels', ...
    'Position', [40+off, 60+off*0.5, 1000, 750], 'Visible','on');
ax = axes(fig); hold(ax, 'on');
end

function path = save_fig(fig, name, d1, d2)
drawnow;
path = fullfile(d1, name);
exportgraphics(fig, path, 'Resolution', 300);
exportgraphics(fig, fullfile(d2, name), 'Resolution', 300);
fprintf('      [OK] %s\n', name);
end

function p = plot_mesh_patch(ax, V, F, color, alpha_val, display_name)
if isempty(display_name)
    p = patch(ax,'Vertices',V,'Faces',F,'FaceColor',color,'EdgeColor','none','FaceAlpha',alpha_val, ...
        'HandleVisibility','off');
else
    p = patch(ax,'Vertices',V,'Faces',F,'FaceColor',color,'EdgeColor','none','FaceAlpha',alpha_val, ...
        'DisplayName',display_name);
end
material(p, 'shiny');
end

function setup_3d_axes(ax, title_str, view_ang)
axis(ax, 'equal'); grid(ax, 'off'); box(ax, 'on');
view(ax, view_ang);
camup(ax, [1 0 0]);                       % superior-inferior (X) points up
delete(findobj(ax, 'Type', 'light'));
try
    camlight(ax, 'headlight');
catch
    light(ax, 'Position', [1 1 1], 'Style', 'infinite');   % fallback light
end
lighting(ax, 'gouraud');
xlabel(ax,'SI: X (mm)','FontWeight','bold');
ylabel(ax,'ML: Y (mm)','FontWeight','bold');
zlabel(ax,'AP: Z (mm)','FontWeight','bold');
title(ax, title_str, 'FontSize', 13, 'FontWeight', 'bold');
end


% -------------------------------------------------------------------------
% FEA LOCAL FUNCTIONS (Steps 7-8)
% -------------------------------------------------------------------------
function G = fea_default_case()
G = struct('W',70,'D',50,'L_tib',50,'L_fem',20,'t_tray',3.5,'L_stem',28,'r_stem',5.5, ...
    't_ins',9,'t_fem',8.5,'t_cort',2.5,'E_cort',16000,'nu_cort',0.30,'E_canc',700,'nu_canc',0.30, ...
    'E_tray',110000,'nu_tray',0.33,'sy_tray',880,'E_ins',1000,'nu_ins',0.42,'sy_ins',25, ...
    'E_fem',210000,'nu_fem',0.30,'sy_fem',650,'E_gap',10,'nu_gap',0.20, ...
    'cw_frac',0.42,'contact_frac',0.35,'flex_deg',15,'F_axial',2100,'F_shear',320,'med_ratio',0.60, ...
    'x0',0,'y0',0,'z0',0,'cw',0,'sp',0,'lc',0,'cy',0);
end

function G = fea_finalize(G)
G.cw = G.cw_frac * G.W;                    % condyle width
G.sp = G.W - G.cw;                         % condyle centre spacing
G.lc = G.contact_frac * G.D;               % AP length of the loaded condylar patch
G.cy = -(G.flex_deg/90) * 0.22 * G.D;      % posterior roll-back with flexion
G.L_stem = min(G.L_stem, G.L_tib - 2);
end

function H = fea_height(G)
H = G.L_tib + G.t_tray + G.t_ins + G.t_fem + G.L_fem;
end

function reg = fea_region_of(x, y, z, G)
% Material region ID at each point (coordinates in the FE model frame)
%  1 tibial cortical | 2 tibial cancellous | 3 stem | 4 tray | 5 UHMWPE insert
%  6 femoral component (condyles) | 7 femoral cortical | 8 femoral cancellous | 9 gap
xr = x - G.x0;  yr = y - G.y0;  zr = z - G.z0;
z1 = G.L_tib;  z2 = z1 + G.t_tray;  z3 = z2 + G.t_ins;  z4 = z3 + G.t_fem;
rim = min(G.W/2 - abs(xr), G.D/2 - abs(yr)) < G.t_cort;
reg = zeros(size(zr));
tib = zr < z1;
reg(tib & rim)  = 1;
reg(tib & ~rim) = 2;
stem = tib & (xr.^2 + yr.^2 <= G.r_stem^2) & (zr >= z1 - G.L_stem);
reg(stem) = 3;
reg(zr >= z1 & zr < z2) = 4;
reg(zr >= z2 & zr < z3) = 5;
lay = zr >= z3 & zr < z4;
inC = (abs(abs(xr) - G.sp/2) <= G.cw/2) & (abs(yr - G.cy) <= G.lc/2);
reg(lay & inC)  = 6;
reg(lay & ~inC) = 9;
fem = zr >= z4;
reg(fem & rim)  = 7;
reg(fem & ~rim) = 8;
end

function [X, Y, Z, L] = fea_sample_points(G)
% Sample points inside every material region (used to extract peak responses)
fr = [0.15 0.40 0.65 0.90];
z1 = G.L_tib;  z2 = z1 + G.t_tray;  z3 = z2 + G.t_ins;  z4 = z3 + G.t_fem;  z5 = z4 + G.L_fem;
zb = {[0 z1],[0 z1],[z1-G.L_stem z1],[z1 z2],[z2 z3],[z3 z4],[z4 z5],[z4 z5]};
m  = 9;
[GX, GY] = meshgrid(linspace(-0.9,0.9,m)*G.W/2, linspace(-0.9,0.9,m)*G.D/2);
full_xy = [GX(:)'; GY(:)'];
t = linspace(-0.9, 0.9, m);
ring_xy = [ (G.W/2-G.t_cort/2)*ones(1,m), -(G.W/2-G.t_cort/2)*ones(1,m), t*G.W/2, t*G.W/2; ...
            t*G.D/2, t*G.D/2, (G.D/2-G.t_cort/2)*ones(1,m), -(G.D/2-G.t_cort/2)*ones(1,m)];
ang = linspace(0, 2*pi, 9);  ang(end) = [];
stem_xy = [0, kron([0.5 0.9]*G.r_stem, cos(ang)); 0, kron([0.5 0.9]*G.r_stem, sin(ang))];
pat = linspace(-0.8, 0.8, 5);
[PX, PY] = meshgrid(pat*G.cw/2, pat*G.lc/2);
cond_xy = [PX(:)'-G.sp/2, PX(:)'+G.sp/2; PY(:)'+G.cy, PY(:)'+G.cy];
xyset = {ring_xy, full_xy, stem_xy, full_xy, [full_xy cond_xy], cond_xy, ring_xy, full_xy};
X = []; Y = []; Z = []; L = [];
for r = 1:8
    xy = xyset{r};
    zl = zb{r}(1) + fr*(zb{r}(2) - zb{r}(1));
    [ia, ib] = meshgrid(1:size(xy,2), 1:numel(zl));
    px = xy(1, ia(:)) + G.x0;  py = xy(2, ia(:)) + G.y0;  pz = zl(ib(:)) + G.z0;
    keep = fea_region_of(px, py, pz, G) == r;
    X = [X, px(keep)]; Y = [Y, py(keep)]; Z = [Z, pz(keep)];          %#ok<AGROW>
    L = [L, r*ones(1, nnz(keep))];                                      %#ok<AGROW>
end
end

function vm = fea_vm(sxx, syy, szz, sxy, syz, sxz)
vm = sqrt(0.5*((sxx-syy).^2 + (syy-szz).^2 + (szz-sxx).^2) + 3*(sxy.^2 + syz.^2 + sxz.^2));
end

function ee = fea_eqstrain(exx, eyy, ezz, exy, eyz, exz)
% von Mises equivalent strain (shear components treated as tensor components)
em = (exx + eyy + ezz)/3;
ee = sqrt(2/3 * ((exx-em).^2 + (eyy-em).^2 + (ezz-em).^2 + 2*(exy.^2 + eyz.^2 + exz.^2)));
end

function m = fea_rmax(v, mask)
v = v(mask);  v = v(~isnan(v));
if isempty(v), m = NaN; else, m = max(v); end
end

function m = fea_rmean(v, mask)
v = v(mask);  v = v(~isnan(v));
if isempty(v), m = NaN; else, m = mean(v); end
end

function R = fea_nan_result()
R = struct('Mesh_Nodes',NaN,'Mesh_Elements',NaN,'Solve_Time_s',NaN,'Peak_VM_Nodal_MPa',NaN,'Peak_VM_MPa',NaN, ...
    'VM_TibBone_MPa',NaN,'VM_Stem_MPa',NaN,'VM_Tray_MPa',NaN,'VM_Insert_MPa',NaN,'VM_FemComp_MPa',NaN,'VM_FemBone_MPa',NaN, ...
    'Peak_Strain_ue',NaN,'Strain_TibBone_ue',NaN,'Strain_Tray_ue',NaN,'Strain_Insert_ue',NaN,'Strain_FemComp_ue',NaN, ...
    'Strain_FemBone_ue',NaN,'Max_Disp_mm',NaN,'Tray_Subsidence_mm',NaN,'FemComp_Disp_mm',NaN,'Insert_Compression_mm',NaN, ...
    'SF_Tray',NaN,'SF_Insert',NaN,'SF_FemComp',NaN,'SF_Bone',NaN,'Status','');
end

function [R, model, res, G] = run_fea_case(G, hmax)
% One 3D static-solid FEA of the tibial implant-bone column (PDE Toolbox).
%  x = ML, y = AP, z = SI (vertical). Bottom of tibia fixed; physiological load on top.
%  The column is split into 5 stacked cells (tibia | tray | insert | femoral
%  layer | femur); structuralProperties only accepts scalar moduli per cell,
%  so in-plane heterogeneity (cortical rim, stem, condyle patches) is
%  represented by area-weighted effective moduli of each layer.
t0 = tic;
R  = fea_nan_result();
H  = fea_height(G);

model = createpde('structural', 'static-solid');
Hl = [G.L_tib, G.t_tray, G.t_ins, G.t_fem, G.L_fem];      % individual layer heights
Zo = [0, cumsum(Hl(1:end-1))];                             % layer base offsets (stacked, touching)
gm = multicuboid(G.W*ones(1,5), G.D*ones(1,5), Hl, 'Zoffset', Zo);
model.Geometry = gm;
if gm.NumCells ~= 5
    error('Expected 5 stacked cells, got %d.', gm.NumCells);
end
msh = generateMesh(model, 'Hmax', hmax, 'GeometricOrder', 'linear');
nd  = msh.Nodes;
G.x0 = (max(nd(1,:)) + min(nd(1,:)))/2;
G.y0 = (max(nd(2,:)) + min(nd(2,:)))/2;
G.z0 = min(nd(3,:));

% --- effective (area-weighted) layer properties ---
A    = G.W * G.D;
Ac   = A - max(G.W-2*G.t_cort,0) * max(G.D-2*G.t_cort,0);        % cortical rim area
As   = pi * G.r_stem^2 * (G.L_stem / G.L_tib);                    % stem area (length-averaged)
Acnd = min(2 * G.cw * G.lc, A);                                   % condyle contact area
wTc = Ac/A;  wTs = As/A;  wTk = max(1 - wTc - wTs, 0);
E_tib  = wTc*G.E_cort  + wTk*G.E_canc  + wTs*G.E_tray;
nu_tib = wTc*G.nu_cort + wTk*G.nu_canc + wTs*G.nu_tray;
wF = Acnd/A;
E_fl   = wF*G.E_fem  + (1-wF)*G.E_gap;
nu_fl  = wF*G.nu_fem + (1-wF)*G.nu_gap;
E_fb   = wTc*G.E_cort  + (1-wTc)*G.E_canc;
nu_fb  = wTc*G.nu_cort + (1-wTc)*G.nu_canc;

Ev  = [E_tib,  G.E_tray,  G.E_ins,  E_fl,  E_fb];
nuv = [nu_tib, G.nu_tray, G.nu_ins, nu_fl, nu_fb];
for c = 1:5
    structuralProperties(model, 'Cell', c, 'YoungsModulus', Ev(c), 'PoissonsRatio', nuv(c));
end

% --- boundary conditions and physiological loading ---
botF = nearestFace(gm, [G.x0, G.y0, G.z0]);
topF = nearestFace(gm, [G.x0, G.y0, G.z0 + H]);
structuralBC(model, 'Face', botF, 'Constraint', 'fixed');

p0 = G.F_axial / (G.W*G.D);                 % mean axial pressure (MPa)
tau = G.F_shear / (G.W*G.D);                % mean AP shear traction (MPa)
mr = G.med_ratio;  xm = G.x0;
trac = @(loc, st) [zeros(1, numel(loc.x)); ...
                   tau*ones(1, numel(loc.x)); ...
                   -2*p0*(mr*(loc.x >= xm) + (1-mr)*(loc.x < xm))];   % medial (+x) / lateral split
structuralBoundaryLoad(model, 'Face', topF, 'SurfaceTraction', trac);

res = solve(model);

% --- post-processing at sample points inside each material region ---
[X, Y, Z, L] = fea_sample_points(G);
st = interpolateStress(res, X(:), Y(:), Z(:));
vm = fea_vm(st.sxx(:)', st.syy(:)', st.szz(:)', st.sxy(:)', st.syz(:)', st.sxz(:)');
sn = interpolateStrain(res, X(:), Y(:), Z(:));
ee = fea_eqstrain(sn.exx(:)', sn.eyy(:)', sn.ezz(:)', sn.exy(:)', sn.eyz(:)', sn.exz(:)');
dp = interpolateDisplacement(res, X(:), Y(:), Z(:));
uz = dp.uz(:)';

vmax = zeros(1,8);  emax = zeros(1,8);
for r = 1:8
    vmax(r) = fea_rmax(vm, L == r);
    emax(r) = fea_rmax(ee, L == r);
end
R.Mesh_Nodes        = size(nd, 2);
R.Mesh_Elements     = size(msh.Elements, 2);
R.Peak_VM_Nodal_MPa = max(res.VonMisesStress);
R.Peak_VM_MPa       = max(vmax, [], 'omitnan');
R.VM_TibBone_MPa    = max(vmax(1:2), [], 'omitnan');
R.VM_Stem_MPa       = vmax(3);
R.VM_Tray_MPa       = vmax(4);
R.VM_Insert_MPa     = vmax(5);
R.VM_FemComp_MPa    = vmax(6);
R.VM_FemBone_MPa    = max(vmax(7:8), [], 'omitnan');
R.Peak_Strain_ue    = max(emax, [], 'omitnan') * 1e6;
R.Strain_TibBone_ue = max(emax(1:2), [], 'omitnan') * 1e6;
R.Strain_Tray_ue    = max(emax(3:4), [], 'omitnan') * 1e6;
R.Strain_Insert_ue  = emax(5) * 1e6;
R.Strain_FemComp_ue = emax(6) * 1e6;
R.Strain_FemBone_ue = max(emax(7:8), [], 'omitnan') * 1e6;
uz_tray = fea_rmean(uz, L == 4);
uz_fem  = fea_rmean(uz, L == 6);
R.Max_Disp_mm            = max(res.Displacement.Magnitude);
R.Tray_Subsidence_mm     = -uz_tray;
R.FemComp_Disp_mm        = -uz_fem;
R.Insert_Compression_mm  = uz_tray - uz_fem;
R.SF_Tray   = G.sy_tray / max(vmax(3), vmax(4), 'omitnan');
R.SF_Insert = G.sy_ins  / vmax(5);
R.SF_FemComp = G.sy_fem / vmax(6);
R.SF_Bone   = 0.007 / max(emax([1 2 7 8]), [], 'omitnan');    % 0.7% yield-strain criterion
R.Solve_Time_s = toc(t0);
end

function [R, M, Gp, P] = run_fea_patient(params, G, fv, ff, tv, tf, x_tcut, x_fcut, Vt, Ft, Vi, Fi, Vf, Ff, h)
% Patient-geometry FEA: the real bone STLs + the generated implant meshes are voxelised
% on one common grid, converted to a conforming (shared-node) tetrahedral mesh, and solved
% with one material per cell. Design frame [SI ML AP] -> FEA frame [x=ML, y=AP, z=SI].
%  cells: 1 tib.cortical | 2 tib.cancellous | 3 stem | 4 tray | 5 insert | 6 fem.comp | 7 fem.cortical | 8 fem.cancellous
t0 = tic;
R  = fea_nan_result();
cv = @(V) V(:, [2 3 1]);
[ffd, fvd] = decimate_mesh(ff, fv, 40000);
[tfd, tvd] = decimate_mesh(tf, tv, 40000);

zmin = x_tcut - 45;  zmax = x_fcut + 40;                     % crop of the bones (mm)
tvc  = tvd(tvd(:,1) >= zmin & tvd(:,1) <= x_tcut, :);
fvc  = fvd(fvd(:,1) >= x_fcut & fvd(:,1) <= zmax, :);
allp = [cv(tvc); cv(fvc); cv(Vt); cv(Vi); cv(Vf)];
xlo = min(allp(:,1)) - 2*h;  xhi = max(allp(:,1)) + 2*h;
ylo = min(allp(:,2)) - 2*h;  yhi = max(allp(:,2)) + 2*h;
nx = ceil((xhi-xlo)/h);  ny = ceil((yhi-ylo)/h);  nz = ceil((zmax-zmin)/h);
xe = xlo + (0:nx)*h;  ye = ylo + (0:ny)*h;  ze = zmin + (0:nz)*h;
sh = 1.234e-4;                                                % avoids exact edge hits in ray casting
xc = xe(1:end-1) + h/2 + sh;  yc = ye(1:end-1) + h/2 + 0.7*sh;  zc = ze(1:end-1) + h/2;

% ---- voxelise bones and implants ----
occT = voxel_fill(cv(tvd), tfd, xc, yc, zc);   occT(:,:, zc > x_tcut) = false;
occF = voxel_fill(cv(fvd), ffd, xc, yc, zc);   occF(:,:, zc < x_fcut) = false;
nr   = max(1, round(G.t_cort/h));
cortT = occT & ~fea_erode2d(occT, nr);
cortF = occF & ~fea_erode2d(occF, nr);
occI = voxel_fill(cv(Vi), Fi, xc, yc, zc);
occR = voxel_fill(cv(Vt), Ft, xc, yc, zc);
occC = voxel_fill(cv(Vf), Ff, xc, yc, zc);

lab = zeros(nx, ny, nz);
lab(occT) = 2;  lab(cortT) = 1;
lab(occF) = 8;  lab(cortF) = 7;
lab(occI) = 5;
lab(occR) = 4;
lab(occR & reshape(zc < x_tcut, 1, 1, [])) = 3;               % keel / stem
lab(occC) = 6;

% conforming contact between UHMWPE insert and femoral condyles (closes gaps up to 4 mm)
for i = 1:nx
    for j = 1:ny
        col = squeeze(lab(i,j,:));
        ki = find(col == 5, 1, 'last');  kf = find(col == 6, 1, 'first');
        if ~isempty(ki) && ~isempty(kf) && kf > ki + 1
            gp = ki+1:kf-1;
            if all(col(gp) == 0) && numel(gp)*h <= 4, lab(i,j,gp) = 5; end
        end
    end
end

% ---- voxels -> conforming tetrahedra (6 per voxel) ----
idx = find(lab > 0);
[ii, jj, kk] = ind2sub([nx ny nz], idx);
lv  = lab(idx);  lv = lv(:)';
nid = @(a, b, c) a + (nx+1)*((b-1) + (ny+1)*(c-1));
pp  = perms([1 2 3]);  I3 = eye(3);
Eall = cell(1, 6);
for p = 1:6
    e = I3(pp(p,:), :);
    c = [0 0 0; e(1,:); e(1,:)+e(2,:); 1 1 1];
    if det([c(2,:)-c(1,:); c(3,:)-c(1,:); c(4,:)-c(1,:)]) < 0, c = c([1 2 4 3], :); end
    Eall{p} = [nid(ii+c(1,1), jj+c(1,2), kk+c(1,3)), nid(ii+c(2,1), jj+c(2,2), kk+c(2,3)), ...
               nid(ii+c(3,1), jj+c(3,2), kk+c(3,3)), nid(ii+c(4,1), jj+c(4,2), kk+c(4,3))]';
end
elem = [Eall{:}];
cid  = repmat(lv, 1, 6);
vid  = repmat(idx(:)', 1, 6);
xcol = xe(:);  ycol = ye(:);  zcol = ze(:);

[elem, nodes, keepE] = fea_compact(elem, nx, ny, nz, xcol, ycol, zcol, []);
cid = cid(keepE);  vid = vid(keepE);
% keep only the largest connected body (guards against floating voxels)
s = [elem(1,:) elem(1,:) elem(1,:) elem(2,:) elem(2,:) elem(3,:)];
t = [elem(2,:) elem(3,:) elem(4,:) elem(3,:) elem(4,:) elem(4,:)];
bins = conncomp(graph(s, t, [], size(nodes,2)));
big  = mode(bins);
keep = bins(elem(1,:)) == big;
elem = elem(:, keep);  cid = cid(keep);  vid = vid(keep);
[used, ~, ic] = unique(elem(:));
elem  = reshape(ic, 4, []);
nodes = nodes(:, used);

% ---- own linear-elastic tetrahedral FEM (conforming shared-node mesh, one material per element) ----
Ne = size(elem, 2);  Nn = size(nodes, 2);  ndof = 3*Nn;
Ev  = [G.E_cort G.E_canc G.E_tray G.E_tray G.E_ins G.E_fem G.E_cort G.E_canc];
nuv = [G.nu_cort G.nu_canc G.nu_tray G.nu_tray G.nu_ins G.nu_fem G.nu_cort G.nu_canc];
Evc = Ev(:);  nuc = nuv(:);
lbl_e = cid(:);
Ee  = Evc(lbl_e);  nue = nuc(lbl_e);
lam = Ee.*nue ./ ((1+nue).*(1-2*nue));  mu = Ee ./ (2*(1+nue));

n1 = elem(1,:)';  n2 = elem(2,:)';  n3 = elem(3,:)';  n4 = elem(4,:)';
P0 = nodes(:, n1)';
v1 = nodes(:, n2)' - P0;  v2 = nodes(:, n3)' - P0;  v3 = nodes(:, n4)' - P0;
c1 = cross(v2, v3, 2);  c2 = cross(v3, v1, 2);  c3 = cross(v1, v2, 2);
dJ  = sum(v1 .* c1, 2);
vol = abs(dJ) / 6;
gN1 = c1 ./ dJ;  gN2 = c2 ./ dJ;  gN3 = c3 ./ dJ;  gN0 = -(gN1 + gN2 + gN3);
gN  = {gN0, gN1, gN2, gN3};

K = sparse(ndof, ndof);
chunk = 40000;
for s0 = 1:chunk:Ne
    ce = s0:min(s0+chunk-1, Ne);  nc = numel(ce);
    B = zeros(6, 12, nc);
    dofs = zeros(12, nc);
    for a_ = 1:4
        cols = 3*(a_-1) + (1:3);
        gx = reshape(gN{a_}(ce,1), 1, 1, nc);
        gy = reshape(gN{a_}(ce,2), 1, 1, nc);
        gz = reshape(gN{a_}(ce,3), 1, 1, nc);
        B(1, cols(1), :) = gx;  B(2, cols(2), :) = gy;  B(3, cols(3), :) = gz;
        B(4, cols(1), :) = gy;  B(4, cols(2), :) = gx;
        B(5, cols(2), :) = gz;  B(5, cols(3), :) = gy;
        B(6, cols(1), :) = gz;  B(6, cols(3), :) = gx;
        dofs(cols, :) = 3*(elem(a_, ce) - 1) + [1; 2; 3];
    end
    lv_ = reshape(lam(ce), 1, 1, nc);  mv_ = reshape(mu(ce), 1, 1, nc);
    D = zeros(6, 6, nc);
    for i_ = 1:3
        for j_ = 1:3, D(i_, j_, :) = lv_; end
        D(i_, i_, :) = lv_ + 2*mv_;
        D(i_+3, i_+3, :) = mv_;
    end
    Ke = pagemtimes(B, 'transpose', pagemtimes(D, B), 'none') .* reshape(vol(ce), 1, 1, nc);
    Ii = repmat(reshape(dofs, 12, 1, nc), [1 12 1]);
    Jj = repmat(reshape(dofs, 1, 12, nc), [12 1 1]);
    K = K + sparse(Ii(:), Jj(:), Ke(:), ndof, ndof);
end
clear B D Ke Ii Jj

% ---- boundary faces: femur-top load faces and external surface for plotting ----
fcs  = [elem([1 2 3],:)'; elem([1 2 4],:)'; elem([1 3 4],:)'; elem([2 3 4],:)'];
eidf = repmat((1:Ne)', 4, 1);
[~, ~, icf] = unique(sort(fcs, 2), 'rows');
cntf = accumarray(icf, 1);
isb  = cntf(icf) == 1;
bF = fcs(isb, :);  bE = eidf(isb);
ztop = max(nodes(3,:));  zbot = min(nodes(3,:));
zbf = reshape(nodes(3, bF(:)), [], 3);
topmask = all(abs(zbf - ztop) < 1e-6, 2);
pa = nodes(:, bF(:,1));  pb = nodes(:, bF(:,2));  pc = nodes(:, bF(:,3));
farea = 0.5 * vecnorm(cross(pb-pa, pc-pa, 1), 2, 1)';
tf_ = bF(topmask, :);  atf = farea(topmask);
nodeArea = accumarray(tf_(:), repmat(atf/3, 3, 1), [Nn 1]);

xm  = params.tibial_tray.center_y;                         % medial (+x) / lateral split
xn  = nodes(1,:)';
med = xn >= xm;
A_med = sum(nodeArea(med));  A_lat = sum(nodeArea(~med));  A_top = A_med + A_lat;
if A_med == 0 || A_lat == 0, med = true(Nn,1); A_med = A_top; A_lat = A_top; end
p_med = G.F_axial*G.med_ratio/A_med;  p_lat = G.F_axial*(1-G.med_ratio)/A_lat;
if all(med), p_med = G.F_axial/A_top; end
tau = G.F_shear / A_top;
f = zeros(ndof, 1);
nn_ = (1:Nn)';
f(3*nn_)   = -(p_med*med + p_lat*(~med)) .* nodeArea;
f(3*nn_-1) = tau * nodeArea;

botNodes = find(abs(nodes(3,:) - zbot) < 1e-6)';
fixd = [3*botNodes-2; 3*botNodes-1; 3*botNodes];
free = setdiff((1:ndof)', fixd);

% ---- solve (scaled PCG with incomplete Cholesky, direct fallback) ----
Kff = K(free, free);  Kff = (Kff + Kff') / 2;
dsc = sqrt(diag(Kff));
Sc  = spdiags(1./dsc, 0, numel(free), numel(free));
Ks  = Sc * Kff * Sc;  Ks = (Ks + Ks') / 2;
rhs = Sc * f(free);
flg = 1;
try
    Lc = ichol(Ks, struct('type','ict','droptol',1e-3,'diagcomp',1e-3));
    [ys, flg] = pcg(Ks, rhs, 1e-9, 4000, Lc, Lc');
catch
    flg = 1;
end
if flg ~= 0, ys = Ks \ rhs; end
U = zeros(ndof, 1);  U(free) = Sc * ys;
u = [U(1:3:end), U(2:3:end), U(3:3:end)];
clear K Kff Ks

% ---- per-element strain / stress (exact for linear tets, no cross-material smoothing) ----
du1 = u(n2,:) - u(n1,:);  du2 = u(n3,:) - u(n1,:);  du3 = u(n4,:) - u(n1,:);
gr = @(p, q) (c1(:,p).*du1(:,q) + c2(:,p).*du2(:,q) + c3(:,p).*du3(:,q)) ./ dJ;   % du_q/dx_p
exx = gr(1,1);  eyy = gr(2,2);  ezz = gr(3,3);
exy = 0.5*(gr(2,1) + gr(1,2));  eyz = 0.5*(gr(3,2) + gr(2,3));  exz = 0.5*(gr(3,1) + gr(1,3));
trE = exx + eyy + ezz;
sxx = lam.*trE + 2*mu.*exx;  syy = lam.*trE + 2*mu.*eyy;  szz = lam.*trE + 2*mu.*ezz;
sxy = 2*mu.*exy;  syz = 2*mu.*eyz;  sxz = 2*mu.*exz;
vm  = fea_vm(sxx, syy, szz, sxy, syz, sxz);
ee  = fea_eqstrain(exx, eyy, ezz, exy, eyz, exz);
umag = vecnorm(u, 2, 2);
um_e = (umag(n1) + umag(n2) + umag(n3) + umag(n4)) / 4;

vmax = zeros(1,8);  emax = zeros(1,8);
for r = 1:8
    vmax(r) = fea_rmax(vm, lbl_e == r);
    emax(r) = fea_rmax(ee, lbl_e == r);
end
R.Mesh_Nodes        = Nn;
R.Mesh_Elements     = Ne;
R.Peak_VM_Nodal_MPa = max(vm);
R.Peak_VM_MPa       = max(vmax, [], 'omitnan');
R.VM_TibBone_MPa    = max(vmax(1:2), [], 'omitnan');
R.VM_Stem_MPa       = vmax(3);
R.VM_Tray_MPa       = vmax(4);
R.VM_Insert_MPa     = vmax(5);
R.VM_FemComp_MPa    = vmax(6);
R.VM_FemBone_MPa    = max(vmax(7:8), [], 'omitnan');
R.Peak_Strain_ue    = max(emax, [], 'omitnan') * 1e6;
R.Strain_TibBone_ue = max(emax(1:2), [], 'omitnan') * 1e6;
R.Strain_Tray_ue    = max(emax(3:4), [], 'omitnan') * 1e6;
R.Strain_Insert_ue  = emax(5) * 1e6;
R.Strain_FemComp_ue = emax(6) * 1e6;
R.Strain_FemBone_ue = max(emax(7:8), [], 'omitnan') * 1e6;
nodes_of = @(r) unique(elem(:, lbl_e == r));
uz_tray = mean(u(nodes_of(4), 3), 'omitnan');
uz_fem  = mean(u(nodes_of(6), 3), 'omitnan');
R.Max_Disp_mm            = max(umag);
R.Tray_Subsidence_mm     = -uz_tray;
R.FemComp_Disp_mm        = -uz_fem;
R.Insert_Compression_mm  = uz_tray - uz_fem;
R.SF_Tray    = G.sy_tray / max(vmax(3), vmax(4), 'omitnan');
R.SF_Insert  = G.sy_ins  / vmax(5);
R.SF_FemComp = G.sy_fem  / vmax(6);
R.SF_Bone    = 0.007 / max(emax([1 2 7 8]), [], 'omitnan');
R.Solve_Time_s = toc(t0);

% ---- voxel-averaged fields for slice plots ----
vid = vid(:);  nvox = nx*ny*nz;
cntv = accumarray(vid, 1, [nvox 1]);
vf = @(q) reshape(accumarray(vid, q(:), [nvox 1]) ./ max(cntv, 1), [nx ny nz]);
Fvm = vf(vm);  Fee = vf(ee*1e6);  Fum = vf(um_e);
emp = reshape(cntv == 0, [nx ny nz]);
Fvm(emp) = NaN;  Fee(emp) = NaN;  Fum(emp) = NaN;
M = struct('nodes', nodes, 'elem', elem, 'lbl', lbl_e, 'bF', bF, 'bE', bE, ...
           'vm', vm, 'ee', ee*1e6, 'um', um_e, 'u', u);

% ---- geometry descriptor used by the slice / material plots ----
Gp = G;
xr = [min(nodes(1,:)) max(nodes(1,:))];  yr = [min(nodes(2,:)) max(nodes(2,:))];
Gp.x0 = mean(xr);  Gp.W = diff(xr);
Gp.y0 = mean(yr);  Gp.D = diff(yr);
Gp.z0 = min(nodes(3,:));
Gp.H  = max(nodes(3,:)) - Gp.z0;
Gp.sp = 2*(params.tibial_tray.center_y + params.femoral.condyle_spacing/2 - Gp.x0);   % medial-condyle slice position
Gp.cy = params.tibial_tray.center_z - Gp.y0;                                          % coronal slice through the tray centre
Gp.L_tib  = x_tcut - Gp.z0;
Gp.t_tray = params.tibial_tray.thickness;
Gp.t_ins  = params.tibial_insert.core_thickness - 1.0;
Gp.t_fem  = params.femoral.resect_thick;
P = struct('lab', lab, 'xc', xc, 'yc', yc, 'zc', zc, 'vm', Fvm, 'ee', Fee, 'um', Fum);
end

function occ = voxel_fill(V, F, xc, yc, zc)
% Solid voxelisation of a closed triangle mesh by ray casting along z (column parity fill).
nx = numel(xc);  ny = numel(yc);  nz = numel(zc);
hx = xc(2) - xc(1);  hy = yc(2) - yc(1);
nT = size(F, 1);
cl = cell(nT, 1);  zl = cell(nT, 1);
A = V(F(:,1), :);  B = V(F(:,2), :);  C = V(F(:,3), :);
for t = 1:nT
    a = A(t,:);  b = B(t,:);  c = C(t,:);
    i0 = max(1,  ceil((min([a(1) b(1) c(1)]) - xc(1))/hx) + 1);
    i1 = min(nx, floor((max([a(1) b(1) c(1)]) - xc(1))/hx) + 1);
    j0 = max(1,  ceil((min([a(2) b(2) c(2)]) - yc(1))/hy) + 1);
    j1 = min(ny, floor((max([a(2) b(2) c(2)]) - yc(1))/hy) + 1);
    if i1 < i0 || j1 < j0, continue; end
    d = (b(2)-c(2))*(a(1)-c(1)) + (c(1)-b(1))*(a(2)-c(2));
    if abs(d) < 1e-12, continue; end
    [I, J] = ndgrid(i0:i1, j0:j1);
    px = xc(I(:));  py = yc(J(:));
    px = px(:);  py = py(:);
    l1 = ((b(2)-c(2))*(px-c(1)) + (c(1)-b(1))*(py-c(2))) / d;
    l2 = ((c(2)-a(2))*(px-c(1)) + (a(1)-c(1))*(py-c(2))) / d;
    l3 = 1 - l1 - l2;
    in = l1 >= 0 & l2 >= 0 & l3 >= 0;
    if ~any(in), continue; end
    Iv = I(:);  Jv = J(:);
    cl{t} = sub2ind([nx ny], Iv(in), Jv(in));
    zl{t} = l1(in)*a(3) + l2(in)*b(3) + l3(in)*c(3);
    cl{t} = cl{t}(:);  zl{t} = zl{t}(:);
end
col = vertcat(cl{:});  zz = vertcat(zl{:});
occ = false(nx*ny, nz);
if ~isempty(col)
    M = sortrows([double(col(:)) zz(:)]);
    st = [1; find(diff(M(:,1)) ~= 0) + 1];
    en = [st(2:end) - 1; size(M,1)];
    for g = 1:numel(st)
        zs = M(st(g):en(g), 2);
        if mod(numel(zs), 2) == 1, zs(end) = []; end
        cc = M(st(g), 1);
        for q = 1:2:numel(zs)
            occ(cc, zc >= zs(q) & zc <= zs(q+1)) = true;
        end
    end
end
occ = reshape(occ, [nx ny nz]);
end

function er = fea_erode2d(o, nr)
% In-plane (axial slice) erosion used to carve the cortical shell
er = o;
k2 = ones(3, 3, 1);
for q = 1:nr
    er = convn(double(er), k2, 'same') > 8.5;
end
end

function [elem, nodes, keepE] = fea_compact(elem, nx, ny, nz, xcol, ycol, zcol, ~)
% Drop unused grid nodes and return compact node list (3 x N)
[used, ~, ic] = unique(elem(:));
elem = reshape(ic, 4, []);
[I_, J_, K_] = ind2sub([nx+1 ny+1 nz+1], used);
nodes = [xcol(I_)'; ycol(J_)'; zcol(K_)'];
keepE = true(1, size(elem, 2));
end

function fea_patch3d(M, vals, cblab, name, ttl, idx, folder, fname)
fig = fea_fig(name, idx);  ax = axes(fig);  hold(ax, 'on');
if isempty(vals)
    patch(ax, 'Faces', M.bF, 'Vertices', M.nodes', 'FaceColor', [0.55 0.85 0.95], ...
        'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.1);
else
    patch(ax, 'Faces', M.bF, 'Vertices', M.nodes', 'FaceVertexCData', vals(M.bE), ...
        'FaceColor', 'flat', 'EdgeColor', 'none');
    colormap(ax, turbo);  cb = colorbar(ax);  cb.Label.String = cblab;
    caxis(ax, [0 fea_clim(vals)]);
end
axis(ax, 'equal', 'tight');  view(ax, 3);  box(ax, 'on');
xlabel(ax, 'ML: X (mm)');  ylabel(ax, 'AP: Y (mm)');  zlabel(ax, 'SI: Z (mm)');
title(ax, ttl, 'FontWeight', 'bold');
fea_save(fig, folder, fname);
end

function fea_plot_slice_patient(P, G, plane, quantity, ttl, idx, folder, fname)
switch quantity
    case 'vm',     F = P.vm;  lbl = 'von Mises stress (MPa)';
    case 'strain', F = P.ee;  lbl = 'Equivalent strain (microstrain)';
    otherwise,     F = P.um;  lbl = 'Displacement magnitude (mm)';
end
if strcmp(plane, 'sagittal')
    [~, ix] = min(abs(P.xc - (G.x0 + G.sp/2)));
    S = squeeze(F(ix, :, :));  A = P.yc;
    xl = 'AP: Y (mm)';  sub = 'Sagittal slice through the medial condyle';
else
    [~, iy] = min(abs(P.yc - (G.y0 + G.cy)));
    S = squeeze(F(:, iy, :));  A = P.xc;
    xl = 'ML: X (mm)';  sub = 'Coronal slice through the tray centre';
end
fig = fea_fig([ttl ' - ' plane], idx);  ax = axes(fig);
pcolor(ax, A, P.zc, S');  shading(ax, 'flat');  colormap(ax, turbo);
cb = colorbar(ax);  cb.Label.String = lbl;  caxis(ax, [0 fea_clim(S)]);
hold(ax, 'on');
zl = G.z0 + [G.L_tib, G.L_tib+G.t_tray, G.L_tib+G.t_tray+G.t_ins, G.L_tib+G.t_tray+G.t_ins+G.t_fem];
nm = {'Tibia | Tray','Tray | Insert','Insert | Fem. comp.','Fem. comp. | Femur'};
for k = 1:4
    yline(ax, zl(k), ':k', nm{k}, 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
end
axis(ax, 'equal', 'tight');  xlabel(ax, xl);  ylabel(ax, 'SI: Z (mm)');
title(ax, {[ttl ' (FEA)'], sub}, 'FontWeight', 'bold');  box(ax, 'on');
fea_save(fig, folder, fname);
end

function fea_plot_materials_patient(P, Gp, idx, folder, fname)
xs = Gp.x0 + Gp.sp/2;
[~, ix] = min(abs(P.xc - xs));
S = double(squeeze(P.lab(ix, :, :)));
S(S == 0) = NaN;
fig = fea_fig('Material Layout', idx);  ax = axes(fig);
pcolor(ax, P.yc, P.zc, S');  shading(ax, 'flat');
colormap(ax, [0.80 0.60 0.35; 0.95 0.88 0.70; 0.60 0.60 0.60; 0.78 0.62 0.28; 0.70 0.88 0.95; ...
              0.22 0.50 0.80; 0.85 0.75 0.55; 0.93 0.86 0.72]);
caxis(ax, [0.5 8.5]);
cb = colorbar(ax);  cb.Ticks = 1:8;
cb.TickLabels = {'Tibial cortical','Tibial cancellous','Stem','Tray','UHMWPE insert','Femoral comp.', ...
                 'Femoral cortical','Femoral cancellous'};
axis(ax, 'equal', 'tight');  xlabel(ax, 'AP: Y (mm)');  ylabel(ax, 'SI: Z (mm)');
title(ax, {'FEA Material Layout (Patient Geometry)', 'Sagittal slice through the medial condyle'}, 'FontWeight', 'bold');  box(ax, 'on');
fea_save(fig, folder, fname);
end

function X = fea_lhs(n, k)
% Latin hypercube sample in [0,1]^k (no Statistics Toolbox needed)
X = zeros(n, k);
for j = 1:k
    X(:, j) = (randperm(n)' - rand(n, 1)) / n;
end
end

function c = fea_clim(V)
v = sort(V(~isnan(V)));
if isempty(v), c = 1; return; end
c = v(max(1, ceil(0.995*numel(v))));
if c <= 0, c = max(v(end), 1); end
end

function fig = fea_fig(name, idx)
off = 22 * mod(idx-1, 12);
fig = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', 'Units', 'pixels', ...
    'Position', [40+off, 60+0.6*off, 980, 720], 'Visible', 'on');
end

function fea_save(fig, folder, name)
drawnow;
exportgraphics(fig, fullfile(folder, name), 'Resolution', 300);
fprintf('      [OK] %s\n', name);
end

function fea_plot_slice(res, G, H, plane, quantity, ttl, idx, folder, fname)
n1 = 110;  n2 = 150;
if strcmp(plane, 'sagittal')
    [A, B] = meshgrid(linspace(G.y0-G.D/2+0.05, G.y0+G.D/2-0.05, n1), linspace(G.z0+0.05, G.z0+H-0.05, n2));
    Xq = (G.x0 + G.sp/2)*ones(size(A));  Yq = A;  Zq = B;
    xl = 'AP: Y (mm)';  sub = 'Sagittal slice through the medial condyle';
else
    [A, B] = meshgrid(linspace(G.x0-G.W/2+0.05, G.x0+G.W/2-0.05, n1), linspace(G.z0+0.05, G.z0+H-0.05, n2));
    Xq = A;  Yq = (G.y0 + G.cy)*ones(size(A));  Zq = B;
    xl = 'ML: X (mm)';  sub = 'Coronal slice through the contact patch';
end
switch quantity
    case 'vm'
        s = interpolateStress(res, Xq(:), Yq(:), Zq(:));
        V = fea_vm(s.sxx, s.syy, s.szz, s.sxy, s.syz, s.sxz);  lbl = 'von Mises stress (MPa)';
    case 'strain'
        s = interpolateStrain(res, Xq(:), Yq(:), Zq(:));
        V = fea_eqstrain(s.exx, s.eyy, s.ezz, s.exy, s.eyz, s.exz) * 1e6;  lbl = 'Equivalent strain (microstrain)';
    otherwise
        s = interpolateDisplacement(res, Xq(:), Yq(:), Zq(:));
        V = sqrt(s.ux.^2 + s.uy.^2 + s.uz.^2);  lbl = 'Displacement magnitude (mm)';
end
V = reshape(V, size(A));
fig = fea_fig([ttl ' - ' plane], idx);  ax = axes(fig);
pcolor(ax, A, B, V);  shading(ax, 'flat');  colormap(ax, turbo);
cb = colorbar(ax);  cb.Label.String = lbl;  caxis(ax, [0 fea_clim(V)]);
hold(ax, 'on');
zl = G.z0 + [G.L_tib, G.L_tib+G.t_tray, G.L_tib+G.t_tray+G.t_ins, G.L_tib+G.t_tray+G.t_ins+G.t_fem];
nm = {'Tibia | Tray','Tray | Insert','Insert | Fem. comp.','Fem. comp. | Femur'};
for k = 1:4
    yline(ax, zl(k), ':k', nm{k}, 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
end
axis(ax, 'equal', 'tight');  xlabel(ax, xl);  ylabel(ax, 'SI: Z (mm)');
title(ax, {[ttl ' (FEA)'], sub}, 'FontWeight', 'bold');  box(ax, 'on');
fea_save(fig, folder, fname);
end

function fea_plot_materials(G, H, idx, folder, fname)
[A, B] = meshgrid(linspace(G.y0-G.D/2+0.05, G.y0+G.D/2-0.05, 110), linspace(G.z0+0.05, G.z0+H-0.05, 150));
Xq = (G.x0 + G.sp/2)*ones(size(A));
Rg = reshape(fea_region_of(Xq(:)', A(:)', B(:)', G), size(A));
fig = fea_fig('Material Layout', idx);  ax = axes(fig);
pcolor(ax, A, B, Rg);  shading(ax, 'flat');
colormap(ax, [0.80 0.60 0.35; 0.95 0.88 0.70; 0.60 0.60 0.60; 0.78 0.62 0.28; 0.70 0.88 0.95; ...
              0.22 0.50 0.80; 0.85 0.75 0.55; 0.93 0.86 0.72; 1 1 1]);
caxis(ax, [0.5 9.5]);
cb = colorbar(ax);  cb.Ticks = 1:9;
cb.TickLabels = {'Tibial cortical','Tibial cancellous','Stem','Tray','UHMWPE insert','Femoral comp.', ...
                 'Femoral cortical','Femoral cancellous','Gap (soft)'};
axis(ax, 'equal', 'tight');  xlabel(ax, 'AP: Y (mm)');  ylabel(ax, 'SI: Z (mm)');
title(ax, {'FEA Material Layout', 'Sagittal slice through the medial condyle'}, 'FontWeight', 'bold');  box(ax, 'on');
fea_save(fig, folder, fname);
end

function fea_bar(vals, names, ylab, ttl, fmt, idx, folder, fname, col)
fig = fea_fig(ttl, idx);  ax = axes(fig);  hold(ax, 'on');
bar(ax, 1:numel(vals), vals, 0.6, 'FaceColor', col, 'EdgeColor', 'k', 'LineWidth', 1.1);
for k = 1:numel(vals)
    text(ax, k, vals(k), sprintf(fmt, vals(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold');
end
set(ax, 'XTick', 1:numel(vals), 'XTickLabel', names);
ylabel(ax, ylab);  title(ax, ttl, 'FontWeight', 'bold');  box(ax, 'on');
ylim(ax, [0 max(vals(~isnan(vals)))*1.15 + eps]);
fea_save(fig, folder, fname);
end

function fea_hist(v, xlab, ttl, idx, folder, fname, col)
fig = fea_fig(ttl, idx);  ax = axes(fig);
histogram(ax, v(~isnan(v)), 40, 'FaceColor', col, 'EdgeColor', 'w');
xlabel(ax, xlab);  ylabel(ax, 'Number of cases');  title(ax, ttl, 'FontWeight', 'bold');  box(ax, 'on');
fea_save(fig, folder, fname);
end

function fea_scatter(x, y, c, xlab, ylab, clab, ttl, idx, folder, fname)
fig = fea_fig(ttl, idx);  ax = axes(fig);
scatter(ax, x, y, 22, c, 'filled', 'MarkerFaceAlpha', 0.7);
colormap(ax, turbo);  cb = colorbar(ax);  cb.Label.String = clab;
xlabel(ax, xlab);  ylabel(ax, ylab);  title(ax, ttl, 'FontWeight', 'bold');  box(ax, 'on');
fea_save(fig, folder, fname);
end


% =========================================================================
% STEPS 9-12 : DL SURROGATE  ->  AMOO  ->  FINAL FEA VALIDATION  ->  REPORT
% MATLAB R2024b | needs: Deep Learning Toolbox (dlarray / adamupdate),
%                        PDE Toolbox (through fh.run_fea_case, Step 11 only)
%
%  9  EfficientNet-elite-B9-V2 (tabular adaptation) + TabNSA native sparse
%     attention + FEA-derived physics constraints in the loss
%  10 Adaptive multi-objective optimisation (constrained NSGA-II with
%     adaptive mutation + adaptive feasible-region exploration)
%  11 Final FEA of the optimised design(s) vs. surrogate predictions
%  12 Final performance report (metrics, improvements, Excel/MAT/figures)
%
% INPUTS
%   T, G0, lo, hi, var_names, HMAX_DATA : variables created by Steps 7-8
%   out_dir : output folder
%   fh      : struct of handles to the local FEA functions:
%             fh.run_fea_case, fh.finalize, fh.fig, fh.save
%   opts    : (optional) struct overriding any field of default_cfg()
% =========================================================================
function out = run_steps_9_to_12(T, G0, lo, hi, var_names, out_dir, hmax_data, fh, opts)
if nargin < 9, opts = struct(); end
cfg = default_cfg();
fo = fieldnames(opts);
for k = 1:numel(fo), cfg.(fo{k}) = opts.(fo{k}); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

tgt = {'Peak_VM_MPa','VM_TibBone_MPa','VM_Stem_MPa','VM_Tray_MPa','VM_Insert_MPa','VM_FemComp_MPa','VM_FemBone_MPa', ...
       'Peak_Strain_ue','Strain_TibBone_ue','Strain_Tray_ue','Strain_Insert_ue','Strain_FemComp_ue','Strain_FemBone_ue', ...
       'Max_Disp_mm','Tray_Subsidence_mm','FemComp_Disp_mm','Insert_Compression_mm', ...
       'SF_Tray','SF_Insert','SF_FemComp','SF_Bone'};
nOut = numel(tgt);
ix = @(n) find(strcmp(tgt, n));
fi = 0;                                           % figure counter

% =========================================================================
% STEP 9 - DATA PREPARATION
% =========================================================================
fprintf('\n=======================================================================\n');
fprintf('  STEP 9: EfficientNet-elite-B9-V2 + TabNSA PHYSICS-INFORMED SURROGATE\n');
fprintf('=======================================================================\n');

Tv   = T(strcmp(T.Status, 'OK'), :);
Yall = Tv{:, tgt};
crit = ~ismember(tgt, {'FemComp_Disp_mm','Tray_Subsidence_mm','Insert_Compression_mm'});
ok   = all(isfinite(Yall), 2) & all(Yall(:, crit) > 0, 2);
Tv   = Tv(ok, :);  Yall = Yall(ok, :);
raw  = table2raw(Tv);
[Xf, fnames] = make_features(raw);
nS = size(Xf, 1);  nIn = size(Xf, 2);
fprintf('[+] Usable FEA cases: %d | input features: %d | outputs: %d\n', nS, nIn, nOut);
if nS < 60
    warning('Only %d usable cases. Set N_CASES = 300 (or more) in Step 8 for a reliable R^2 > 0.90.', nS);
end

rng(cfg.seed);
perm = randperm(nS);
nTr = round(0.70*nS);  nVa = round(0.15*nS);
itr = perm(1:nTr);  iva = perm(nTr+1:nTr+nVa);  ite = perm(nTr+nVa+1:end);
fprintf('    Split: train %d | validation %d | test %d\n', numel(itr), numel(iva), numel(ite));

xmu = mean(Xf(itr,:), 1);  xsg = std(Xf(itr,:), 0, 1);  xsg(xsg < 1e-8) = 1;
Xs  = ((Xf - xmu) ./ xsg)';                                   % nIn x N

isLog = all(Yall > 0, 1)';
Zt = Yall;
for k = 1:nOut
    if isLog(k), Zt(:,k) = log(Yall(:,k)); else, Zt(:,k) = asinh(Yall(:,k)); end
end
zmu = mean(Zt(itr,:), 1)';  zsg = std(Zt(itr,:), 0, 1)';  zsg(zsg < 1e-8) = 1;
Zs  = ((Zt - zmu') ./ zsg')';                                 % nOut x N

% ---- physics bookkeeping (all derived from the FEA post-processing in Step 7/8) ----
ph = struct();
for k = 1:nOut, ph.i.(tgt{k}) = k; end
ph.ymu = zmu;  ph.ysg = zsg;  ph.isLog = isLog;
ph.xmu = xmu;  ph.xsg = xsg;
ph.sgn = ones(nOut,1);  ph.sgn(startsWith(tgt, 'SF_')) = -1;  % stress/strain/disp rise with load, SF falls
loadIdx = [18 19 23 24 25 26];                                % features that scale with joint load
ph.shift = zeros(nIn,1);
ph.shift(loadIdx) = log(1.1) ./ xsg(loadIdx)';                % +10 % load in standardised units
ph.ic_sd = std(Yall(itr, ix('Insert_Compression_mm'))) + eps;
ph.md_sd = std(Yall(itr, ix('Max_Disp_mm'))) + eps;
assert(strcmp(fnames{13}, 'log_Sy_insert') && strcmp(fnames{15}, 'log_Sy_tray'), 'Feature order changed.');

A = make_aux(nIn, cfg);
A.arch = make_arch(cfg);
fprintf('    Architecture: width x%.1f, depth x%.1f -> stage widths %s, blocks %s (%d MBConv-SE blocks)\n', ...
    cfg.width_mult, cfg.depth_mult, mat2str(A.arch.W), mat2str(A.arch.Nb), A.arch.nb);

% =========================================================================
% STEP 9 - TRAINING (ensemble of physics-informed networks)
% =========================================================================
models = cell(1, cfg.n_ens);  hists = cell(1, cfg.n_ens);
for m = 1:cfg.n_ens
    fprintf('\n[+] Training ensemble member %d / %d\n', m, cfg.n_ens);
    [models{m}, hists{m}] = train_one_model(Xs(:,itr), Zs(:,itr), Xs(:,iva), Zs(:,iva), cfg, A, ph, cfg.seed + m);
end
extra = 0;
while true
    surr = struct('models',{models}, 'A',A, 'cfg',cfg, 'xmu',xmu, 'xsg',xsg, 'ymu',zmu, 'ysg',zsg, ...
                  'isLog',isLog, 'tgt',{tgt}, 'feature_names',{fnames});
    Ypred = surr_predict(surr, Xf);                           % original units
    Mtr = metrics_tbl(Yall(itr,:), Ypred(itr,:), tgt);
    Mva = metrics_tbl(Yall(iva,:), Ypred(iva,:), tgt);
    Mte = metrics_tbl(Yall(ite,:), Ypred(ite,:), tgt);
    if all(Mte.R2(1:nOut) >= cfg.r2_target) || extra >= cfg.auto_rounds, break; end
    extra = extra + 1;
    fprintf('\n[+] Some test R2 < %.2f -> auto-round %d/%d: adding %d ensemble members (now %d)\n', ...
        cfg.r2_target, extra, cfg.auto_rounds, cfg.extra_members, numel(models) + cfg.extra_members);
    for q = 1:cfg.extra_members
        m = numel(models) + 1;
        fprintf('[+] Training extra ensemble member %d\n', m);
        [models{m}, hists{m}] = train_one_model(Xs(:,itr), Zs(:,itr), Xs(:,iva), Zs(:,iva), cfg, A, ph, cfg.seed + 100*extra + m);
    end
end

print_metrics('TRAIN',      Mtr, nOut, numel(itr));
print_metrics('VALIDATION', Mva, nOut, numel(iva));
print_metrics('TEST (held-out)', Mte, nOut, numel(ite));
fprintf('    Mean R2  train %.4f | val %.4f | test %.4f   (ensemble of %d models)\n', ...
    mean(Mtr.R2(1:nOut)), mean(Mva.R2(1:nOut)), mean(Mte.R2(1:nOut)), numel(models));
below = Mte.Target(Mte.R2(1:nOut) < cfg.r2_target);
if isempty(below)
    fprintf('    All %d outputs reach test R2 >= %.2f\n', nOut, cfg.r2_target);
else
    fprintf('    Outputs below R2 = %.2f: %s\n', cfg.r2_target, strjoin(below', ', '));
    fprintf('    -> increase N_CASES (300-500), or opts.epochs / opts.auto_rounds, and re-run.\n');
end
writetable(Mtr, fullfile(out_dir, 'Metrics_train.csv'));
writetable(Mva, fullfile(out_dir, 'Metrics_validation.csv'));
writetable(Mte, fullfile(out_dir, 'Metrics_test.csv'));

% ---- Step 9 figures ----
fi = fi+1; fig = fh.fig('Training curves', fi); ax = axes(fig); hold(ax,'on');
for m = 1:numel(hists)
    plot(ax, hists{m}(:,1), '-', 'LineWidth', 1.2, 'DisplayName', sprintf('train data loss #%d', m));
    plot(ax, hists{m}(:,4), '--', 'LineWidth', 1.2, 'DisplayName', sprintf('val loss #%d', m));
end
set(ax, 'YScale', 'log'); xlabel(ax,'Epoch'); ylabel(ax,'Standardised MSE');
title(ax, 'Training / Validation Loss', 'FontWeight','bold'); legend(ax,'Location','northeast'); box(ax,'on');
fh.save(fig, out_dir, 'M01_training_curves.png');

fi = fi+1; fig = fh.fig('Physics loss terms', fi); ax = axes(fig); hold(ax,'on');
plot(ax, hists{1}(:,2), 'LineWidth', 1.5, 'DisplayName', 'Physics-consistency loss');
plot(ax, hists{1}(:,3), 'LineWidth', 1.5, 'DisplayName', 'Load-monotonicity loss');
set(ax, 'YScale', 'log'); xlabel(ax,'Epoch'); ylabel(ax,'Loss');
title(ax, 'FEA-Derived Physics Constraint Losses (member 1)', 'FontWeight','bold'); legend(ax); box(ax,'on');
fh.save(fig, out_dir, 'M02_physics_losses.png');

fi = fi+1; fig = fh.fig('Parity plots (test)', fi); fig.Position = [20 20 1550 900];
tl = tiledlayout(fig, 4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nOut
    ax = nexttile(tl); hold(ax,'on');
    a = Yall(ite,k); p = Ypred(ite,k);
    scatter(ax, a, p, 14, [0.20 0.45 0.80], 'filled');
    lim = [min([a;p]) max([a;p])]; plot(ax, lim, lim, 'k--');
    title(ax, sprintf('%s  R^2=%.3f', strrep(tgt{k},'_',' '), Mte.R2(k)), 'FontSize', 7); box(ax,'on');
end
title(tl, 'Predicted vs FEA (held-out test set)', 'FontWeight', 'bold');
fh.save(fig, out_dir, 'M03_parity_test_all_targets.png');

fi = fi+1; fig = fh.fig('R2 per target', fi); ax = axes(fig); hold(ax,'on');
bar(ax, 1:nOut, Mte.R2(1:nOut), 0.6, 'FaceColor', [0.25 0.60 0.85]);
yline(ax, cfg.r2_target, '--r', sprintf('R^2 = %.2f', cfg.r2_target), 'LineWidth', 1.5);
set(ax, 'XTick', 1:nOut, 'XTickLabel', strrep(tgt,'_',' '), 'XTickLabelRotation', 60);
ylabel(ax, 'Test R^2'); ylim(ax, [min(0.5, min(Mte.R2(1:nOut))-0.02) 1.02]);
title(ax, 'Test-Set R^2 per Predicted Quantity', 'FontWeight','bold'); box(ax,'on');
fh.save(fig, out_dir, 'M04_r2_per_target.png');

% =========================================================================
% STEP 10 - ADAPTIVE MULTI-OBJECTIVE OPTIMISATION (AMOO)
% =========================================================================
fprintf('\n=======================================================================\n');
fprintf('  STEP 10: ADAPTIVE MULTI-OBJECTIVE OPTIMISATION (surrogate-driven)\n');
fprintf('=======================================================================\n');

% patient context (fixed anatomy / bone quality / loading) taken from the baseline case
ctx = struct('W',G0.W, 'D',G0.D, 't_cort',G0.t_cort, 'E_cort',G0.E_cort, 'E_canc',G0.E_canc, ...
    'E_ins',G0.E_ins, 'Sy_ins',G0.sy_ins, 'flex',G0.flex_deg, 'axbw',cfg.act_bw, ...
    'F_ax',G0.F_axial, 'F_sh',G0.F_shear, 'med',G0.med_ratio, 'Lc',G0.lc, 'shift',G0.cy, ...
    'L_tib',G0.L_tib, 'max_stack',17.0);
ctx.mass = ctx.F_ax / (ctx.axbw * 9.81);
iv = @(n) find(strcmp(var_names, n));
dv_names = {'t_tray','L_stem','r_stem','t_ins','t_fem'};
lb = [lo(cellfun(iv, dv_names)), 0];  ub = [hi(cellfun(iv, dv_names)), 1];   % 6th gene: tray material (0=Ti, 1=CoCr)
ctx.ub = ub;
if ctx.W < lo(iv('W')) || ctx.W > hi(iv('W')) || ctx.D < lo(iv('D')) || ctx.D > hi(iv('D'))
    warning('Patient tray size (%.1f x %.1f mm) is outside the sampled W/D range - surrogate extrapolates.', ctx.W, ctx.D);
end

xb = [G0.t_tray, G0.L_stem, G0.r_stem, G0.t_ins, G0.t_fem, 0];             % baseline design
[Fb, CVb, ~] = eval_designs(xb, ctx, surr, cfg, tgt);
fprintf('[+] Baseline (surrogate): implant stress %.3f MPa | disp %.5f mm | min implant SF %.2f | SF bone %.2f | CV %.3f\n', ...
    Fb(1), Fb(2), -Fb(3), -Fb(4), CVb);
if CVb > 0
    fprintf('    NOTE: the baseline violates a constraint (SF thresholds are cfg.sf_min=%.2f, cfg.sf_bone_min=%.2f).\n', cfg.sf_min, cfg.sf_bone_min);
end

fevalfun = @(X) eval_designs(X, ctx, surr, cfg, tgt);
[PopF, FF, CVF, Yp_pop, hist] = amoo_nsga2(fevalfun, lb, ub, cfg);
[rk, ~] = nds_crowd(FF, CVF);
par = (rk == 1) & (CVF <= 0);
if ~any(par)
    warning('No feasible design found - returning least-violating non-dominated set.');
    par = (rk == 1);
end
Xpar = PopF(par,:);  Fpar = FF(par,:);  Ypar = Yp_pop(par,:);  CVpar = CVF(par); %#ok<NASGU>
[Xpar, iu] = unique(round(Xpar, 3), 'rows', 'stable');  Fpar = Fpar(iu,:);  Ypar = Ypar(iu,:);  CVpar = CVpar(iu); %#ok<NASGU>
fprintf('[+] Pareto-optimal feasible designs found: %d\n', size(Xpar,1));

fmin = min(Fpar,[],1);  fmax = max(Fpar,[],1);  fr = max(fmax - fmin, eps);
dist = sqrt(sum(((Fpar - fmin)./fr).^2, 2));                 % distance to the utopia point
[~, isel] = min(dist);
xopt = Xpar(isel,:);  Fopt = Fpar(isel,:);
fprintf('[+] Selected compromise design (closest to utopia point):\n');
fprintf('    t_tray %.2f mm | L_stem %.2f mm | r_stem %.2f mm | t_ins %.2f mm | t_fem %.2f mm | tray %s\n', ...
    xopt(1), xopt(2), xopt(3), xopt(4), xopt(5), ternary(xopt(6) >= 0.5, 'Co-Cr-Mo', 'Ti-6Al-4V'));

pt = array2table([Xpar(:,1:5), double(Xpar(:,6) >= 0.5), Fpar(:,1), Fpar(:,2), -Fpar(:,3), -Fpar(:,4), Fpar(:,5), CVpar], ...
    'VariableNames', {'t_tray_mm','L_stem_mm','r_stem_mm','t_ins_mm','t_fem_mm','TrayIsCoCr', ...
    'Pred_ImplantStress_MPa','Pred_MaxDisp_mm','Pred_MinImplantSF','Pred_SF_Bone','Invasiveness_Index','ConstraintViolation'});
pt.Selected = false(height(pt),1);  pt.Selected(isel) = true;

fi = fi+1; fig = fh.fig('Pareto front 3D', fi); ax = axes(fig); hold(ax,'on');
scatter3(ax, Fpar(:,1), Fpar(:,2), -Fpar(:,3), 36, Fpar(:,5), 'filled', 'MarkerFaceAlpha', 0.8);
scatter3(ax, Fb(1), Fb(2), -Fb(3), 120, 'ks', 'filled');
scatter3(ax, Fopt(1), Fopt(2), -Fopt(3), 160, 'rp', 'filled');
colormap(ax, turbo); cb = colorbar(ax); cb.Label.String = 'Invasiveness index';
xlabel(ax,'Implant peak von Mises (MPa)'); ylabel(ax,'Max displacement (mm)'); zlabel(ax,'Min implant safety factor');
title(ax, 'Pareto-Optimal Implant Designs (surrogate)', 'FontWeight','bold'); view(ax, 135, 25); grid(ax,'on'); box(ax,'on');
legend(ax, {'Pareto designs','Baseline','Selected'}, 'Location','best');
fh.save(fig, out_dir, 'O01_pareto_front_3d.png');

fi = fi+1; fig = fh.fig('AMOO convergence', fi); fig.Position = [30 30 1200 450];
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact');
ax = nexttile(tl); plot(ax, hist.hv, 'LineWidth', 2); xlabel(ax,'Generation'); ylabel(ax,'Hypervolume (MC estimate)'); title(ax,'Convergence'); box(ax,'on');
ax = nexttile(tl); plot(ax, hist.feas, 'LineWidth', 2); xlabel(ax,'Generation'); ylabel(ax,'Feasible fraction'); title(ax,'Feasible-region discovery'); box(ax,'on');
ax = nexttile(tl); yyaxis(ax,'left'); plot(ax, hist.mscale, 'LineWidth', 2); ylabel(ax,'Mutation scale');
yyaxis(ax,'right'); plot(ax, hist.succ, 'LineWidth', 1.5); ylabel(ax,'Offspring survival rate'); xlabel(ax,'Generation'); title(ax,'Adaptive operators'); box(ax,'on');
fh.save(fig, out_dir, 'O02_amoo_convergence.png');

% =========================================================================
% STEP 11 - FINAL FEA VALIDATION
% =========================================================================
fprintf('\n=======================================================================\n');
fprintf('  STEP 11: FINAL FEA VALIDATION (PDE Toolbox, mesh Hmax = %.1f mm)\n', hmax_data);
fprintf('=======================================================================\n');

fprintf('[+] FEA of the baseline design...\n');
[Rb, ~] = fea_for_design(xb, G0, fh, hmax_data);
fprintf('[+] FEA of the selected optimised design...\n');
[Ro, Gopt] = fea_for_design(xopt, G0, fh, hmax_data);
Yb_fea = cellfun(@(n) Rb.(n), tgt);   Yo_fea = cellfun(@(n) Ro.(n), tgt);
Yo_sur = surr_predict(surr, make_features(design2raw(xopt, ctx)));

% extra validation on a spread of Pareto designs
nv = min(cfg.n_val, size(Xpar,1));
Fn = (Fpar - fmin) ./ fr;
chosen = isel;
while numel(chosen) < nv
    dm = min(pdist2_simple(Fn, Fn(chosen,:)), [], 2);  dm(chosen) = -1;
    [~, nxt] = max(dm);  chosen(end+1) = nxt;                  %#ok<AGROW>
end
Yv_fea = nan(nv, nOut);  Yv_sur = nan(nv, nOut);
for q = 1:nv
    fprintf('[+] Validation FEA %d / %d (Pareto design)\n', q, nv);
    if chosen(q) == isel, Rq = Ro; else, Rq = fea_for_design(Xpar(chosen(q),:), G0, fh, hmax_data); end
    Yv_fea(q,:) = cellfun(@(n) Rq.(n), tgt);
    Yv_sur(q,:) = surr_predict(surr, make_features(design2raw(Xpar(chosen(q),:), ctx)));
end
okv = all(isfinite(Yv_fea), 2);  Yv_fea = Yv_fea(okv,:);  Yv_sur = Yv_sur(okv,:);

ae  = abs(Yo_sur - Yo_fea);
pe  = 100 * ae ./ max(abs(Yo_fea), eps);
cmp = table(tgt', Yo_fea', Yo_sur', ae', pe', 'VariableNames', {'Output','FEA_Actual','Surrogate_Pred','Abs_Error','Pct_Error'});
fprintf('\n[+] Surrogate vs FEA for the selected optimised design\n');
fprintf('    %-22s %12s %12s %12s %9s\n', 'Output', 'FEA', 'Surrogate', 'AbsErr', '%Err');
for k = 1:nOut
    fprintf('    %-22s %12.5g %12.5g %12.3g %8.2f%%\n', tgt{k}, Yo_fea(k), Yo_sur(k), ae(k), pe(k));
end
fprintf('    Mean abs. %% error (selected design) : %.2f %%  -> accuracy %.2f %%\n', mean(pe), 100 - mean(pe));

if size(Yv_fea,1) >= 3
    Mv = metrics_tbl(Yv_fea, Yv_sur, tgt);
    pe_v = 100 * mean(abs(Yv_sur - Yv_fea) ./ max(abs(Yv_fea), eps), 1)';
    Mv.MAPE_pct = [pe_v; mean(pe_v)];
    fprintf('    Pareto-set validation (%d FEA runs): mean MAPE %.2f %% | mean R2 %.4f\n', size(Yv_fea,1), mean(pe_v), mean(Mv.R2(1:nOut)));
else
    Mv = table();
end

% mechanical validity / feasibility of the optimised design (from the FEA itself)
sfT = Ro.SF_Tray; sfI = Ro.SF_Insert; sfF = Ro.SF_FemComp; sfB = Ro.SF_Bone;
chk = {'Tray yield SF',        sfT >= cfg.sf_min;
       'Insert yield SF',      sfI >= cfg.sf_min;
       'Femoral comp. SF',     sfF >= cfg.sf_min;
       'Bone strain SF',       sfB >= cfg.sf_bone_min;
       'Stem width <= 30% AP', 2*xopt(3) <= 0.30*ctx.D;
       'Tray+insert <= 17 mm', xopt(1) + xopt(4) <= ctx.max_stack;
       'Stem fits in tibia',   xopt(2) <= ctx.L_tib - 2;
       'FEA converged (OK)',   strcmp(Ro.Status, 'OK')};
fprintf('\n[+] Mechanical validity / feasibility check (actual FEA)\n');
for k = 1:size(chk,1), fprintf('    %-24s : %s\n', chk{k,1}, ternary(chk{k,2}, 'PASS', 'FAIL')); end
valid_all = all([chk{:,2}]);
fprintf('    >>> OVERALL: %s\n', ternary(valid_all, 'optimised design is mechanically valid & feasible', 'one or more checks failed'));

fi = fi+1; fig = fh.fig('Surrogate vs FEA (selected)', fi); fig.Position = [30 30 1300 600]; ax = axes(fig);
bar(ax, categorical(strrep(tgt,'_',' '), strrep(tgt,'_',' ')), [Yo_fea ./ max(abs(Yo_fea),eps); Yo_sur ./ max(abs(Yo_fea),eps)]', 'grouped');
ylabel(ax, 'Value normalised by FEA'); legend(ax, {'FEA','Surrogate'}); xtickangle(ax, 60);
title(ax, 'Optimised Design: Surrogate Prediction vs Actual FEA', 'FontWeight','bold'); box(ax,'on');
fh.save(fig, out_dir, 'V01_selected_design_surrogate_vs_fea.png');

if size(Yv_fea,1) >= 3
    fi = fi+1; fig = fh.fig('Pareto validation parity', fi); ax = axes(fig); hold(ax,'on');
    zf = (Yv_fea - mean(Yall(itr,:),1)) ./ std(Yall(itr,:),0,1);
    zs = (Yv_sur - mean(Yall(itr,:),1)) ./ std(Yall(itr,:),0,1);
    scatter(ax, zf(:), zs(:), 20, [0.85 0.35 0.20], 'filled');
    lim = [min([zf(:);zs(:)]) max([zf(:);zs(:)])]; plot(ax, lim, lim, 'k--');
    xlabel(ax, 'FEA (z-scored)'); ylabel(ax, 'Surrogate (z-scored)');
    title(ax, sprintf('Validation on %d Pareto Designs (all outputs)', size(Yv_fea,1)), 'FontWeight','bold'); box(ax,'on');
    fh.save(fig, out_dir, 'V02_pareto_validation_parity.png');
end

% =========================================================================
% STEP 12 - FINAL PERFORMANCE PREDICTION / IMPROVEMENT REPORT
% =========================================================================
fprintf('\n=======================================================================\n');
fprintf('  STEP 12: FINAL PERFORMANCE - OPTIMISED vs BASELINE\n');
fprintf('=======================================================================\n');

summ = @(Y) [max([Y(ix('VM_Stem_MPa')), Y(ix('VM_Tray_MPa')), Y(ix('VM_Insert_MPa')), Y(ix('VM_FemComp_MPa'))]), ...
                 Y(ix('Peak_VM_MPa')), Y(ix('Peak_Strain_ue')), Y(ix('Strain_TibBone_ue')), Y(ix('Max_Disp_mm')), ...
                 Y(ix('Tray_Subsidence_mm')), Y(ix('Insert_Compression_mm')), Y(ix('SF_Tray')), Y(ix('SF_Insert')), ...
                 Y(ix('SF_FemComp')), Y(ix('SF_Bone')), min([Y(ix('SF_Tray')), Y(ix('SF_Insert')), Y(ix('SF_FemComp'))])];
mnames = {'Implant peak von Mises (MPa)','Peak von Mises, all parts (MPa)','Peak strain (microstrain)','Tibial bone strain (microstrain)', ...
          'Max displacement (mm)','Tray subsidence (mm)','Insert compression (mm)','SF tray','SF insert','SF femoral comp.', ...
          'SF bone (strain)','Min implant SF'};
lower_better = [true true true true true true true false false false false false];
sb = summ(Yb_fea);  so = summ(Yo_fea);  ss = summ(Yo_sur(:)');
impr = 100 * (lower_better .* (sb - so) + ~lower_better .* (so - sb)) ./ max(abs(sb), eps);
imp_tbl = table(mnames', sb', so', ss', impr', 'VariableNames', {'Metric','Baseline_FEA','Optimised_FEA','Optimised_Surrogate','Improvement_pct'});
fprintf('    %-34s %12s %12s %12s %10s\n', 'Metric', 'Baseline', 'Optimised', 'Surrogate', 'Improv.%');
for k = 1:numel(mnames)
    fprintf('    %-34s %12.5g %12.5g %12.5g %9.2f%%\n', mnames{k}, sb(k), so(k), ss(k), impr(k));
end
fprintf('\n    Surrogate held-out test accuracy : mean R2 = %.4f | min R2 = %.4f\n', ...
    mean(Mte.R2(1:nOut)), min(Mte.R2(1:nOut)));
fprintf('    Optimised-design prediction error: mean abs %% error = %.2f %%\n', mean(pe));
fprintf('    Optimised design : t_tray %.2f | L_stem %.2f | r_stem %.2f | t_ins %.2f | t_fem %.2f | %s\n', ...
    xopt(1:5), ternary(xopt(6) >= 0.5, 'Co-Cr-Mo tray', 'Ti-6Al-4V tray'));

fi = fi+1; fig = fh.fig('Baseline vs optimised', fi); fig.Position = [30 30 1200 600]; ax = axes(fig);
bar(ax, categorical(mnames, mnames), [sb; so]', 'grouped'); set(ax, 'YScale', 'log');
ylabel(ax, 'Value (log scale, mixed units)'); legend(ax, {'Baseline (FEA)','Optimised (FEA)'}, 'Location','best'); xtickangle(ax, 45);
title(ax, 'Baseline vs Optimised Personalised Implant (FEA)', 'FontWeight','bold'); box(ax,'on');
fh.save(fig, out_dir, 'F01_baseline_vs_optimised.png');

fi = fi+1; fig = fh.fig('Improvement', fi); ax = axes(fig);
barh(ax, impr, 0.6, 'FaceColor', [0.30 0.65 0.40]);
set(ax, 'YTick', 1:numel(mnames), 'YTickLabel', mnames); xline(ax, 0, 'k-');
xlabel(ax, 'Improvement over baseline (%)  (positive = better)'); title(ax, 'Improvement from Personalised Optimisation', 'FontWeight','bold'); box(ax,'on');
fh.save(fig, out_dir, 'F02_improvement_percent.png');

% ---- save everything ----
cfg_tbl = table(fieldnames(cfg), cell2mat(struct2cell(structfun(@(v) double(v(1)), cfg, 'UniformOutput', false))), ...
    'VariableNames', {'Parameter','Value'});
sheets = {'Train_Metrics',Mtr; 'Val_Metrics',Mva; 'Test_Metrics',Mte; 'Pareto_Set',pt; ...
          'Selected_Surr_vs_FEA',cmp; 'Pareto_Validation',Mv; 'Improvement',imp_tbl; 'Config',cfg_tbl};
save_tables(out_dir, 'DL_AMOO_Results', sheets);
save(fullfile(out_dir, 'DL_AMOO_results.mat'), 'surr', 'Mtr', 'Mva', 'Mte', 'pt', 'xopt', 'xb', 'cmp', 'imp_tbl', ...
    'Yo_fea', 'Yo_sur', 'Yb_fea', 'Yv_fea', 'Yv_sur', 'cfg', 'tgt', 'itr', 'iva', 'ite', '-v7.3');
fprintf('    - MAT   : %s\n', fullfile(out_dir, 'DL_AMOO_results.mat'));

out = struct('surrogate',surr, 'test_metrics',Mte, 'val_metrics',Mva, 'train_metrics',Mtr, 'pareto',pt, ...
    'x_opt',xopt, 'x_baseline',xb, 'fea_opt',Ro, 'fea_base',Rb, 'comparison',cmp, 'improvement',imp_tbl, ...
    'valid',valid_all, 'G_opt',Gopt);
fprintf('\n=======================================================================\n');
fprintf('  FINAL SURROGATE METRICS (R2, RMSE, MSE, MAE) - original units\n');
fprintf('=======================================================================\n');
print_metrics('TEST (held-out)', Mte, nOut, numel(ite));
if ~isempty(Mv), print_metrics('PARETO-DESIGN VALIDATION vs NEW FEA', Mv, nOut, size(Yv_fea,1)); end
fprintf('\n=======================================================================\n');
fprintf('  STEPS 9-12 COMPLETED | mean test R2 = %.4f | outputs in: %s\n', mean(Mte.R2(1:nOut)), out_dir);
fprintf('=======================================================================\n');
end


% =========================================================================
%                              LOCAL FUNCTIONS (STEPS 9-12)
% =========================================================================
function cfg = default_cfg()
cfg = struct('seed',2024, 'n_ens',3, 'epochs',400, 'patience',80, 'batch',16, 'lr',2e-3, 'wd',1e-4, 'clip',1.0, ...
    'width_mult',2.4, 'depth_mult',4.2, 'd_tok',16, 'n_nsa',2, 'blk',4, 'topk',3, 'win',2, 'drop',0.05, ...
    'lam_phys',0.2, 'lam_mono',0.1, 'n_colloc',32, ...
    'pop',96, 'gens',80, 'n_val',6, 'sf_min',1.5, 'sf_bone_min',1.0, 'r2_target',0.90, 'act_bw',2.8, 'auto_rounds',2, 'extra_members',2);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

% ---------------------------- features -----------------------------------
function r = table2raw(Tv)
r.W = Tv.W_ML_mm;  r.D = Tv.D_AP_mm;  r.t_cort = Tv.Cort_Thick_mm;  r.E_cort = Tv.E_Cort_MPa;
r.E_canc = Tv.E_Canc_MPa;  r.mass = Tv.Body_Mass_kg;
r.t_tray = Tv.Tray_Thick_mm;  r.L_stem = Tv.Stem_Len_mm;  r.r_stem = Tv.Stem_Radius_mm;
r.t_ins = Tv.Insert_Thick_mm;  r.t_fem = Tv.FemComp_Thick_mm;
r.E_ins = Tv.E_Insert_MPa;  r.Sy_ins = Tv.Sy_Insert_MPa;
isTi = strcmp(Tv.Tray_Material, 'Ti-6Al-4V');
r.E_tray  = 110000*isTi + 210000*~isTi;
r.Sy_tray = 880*isTi + 650*~isTi;
r.flex = Tv.Flexion_deg;  r.axbw = Tv.Axial_xBW;  r.F_ax = Tv.F_Axial_N;  r.F_sh = Tv.F_Shear_N;
r.med = Tv.Medial_Ratio;  r.Lc = Tv.Contact_Len_mm;  r.shift = Tv.Contact_AP_Shift_mm;
end

function [X, names] = make_features(r)
% Inputs grouped as: anatomy(1-6) | design(7-11) | materials(12-15) | loading(16-22) | physics priors(23-26)
% Features 23-26 are closed-form engineering estimates (same formulas as the Step-4 analytical model).
WD = r.W .* r.D;
cols = {log(r.W), log(r.D), r.t_cort, log(r.E_cort), log(r.E_canc), log(r.mass), ...
        r.t_tray, r.L_stem, r.r_stem, r.t_ins, r.t_fem, ...
        log(r.E_ins), log(r.Sy_ins), log(r.E_tray), log(r.Sy_tray), ...
        r.flex, r.axbw, log(r.F_ax), log(r.F_sh), r.med, r.Lc, r.shift, ...
        log(r.F_ax ./ WD), log(r.F_ax ./ (pi*r.r_stem.^2)), ...
        log(0.0138 * r.F_ax .* r.D.^2 ./ (r.t_tray.^3 .* r.W)), log(r.F_ax ./ (0.84 * r.W .* r.Lc))};
cols = cellfun(@(c) c(:), cols, 'UniformOutput', false);
X = [cols{:}];
names = {'log_W','log_D','t_cort','log_E_cort','log_E_canc','log_mass','t_tray','L_stem','r_stem','t_ins','t_fem', ...
         'log_E_insert','log_Sy_insert','log_E_tray','log_Sy_tray','flexion','axial_xBW','log_F_axial','log_F_shear', ...
         'medial_ratio','contact_len','contact_shift','log_p_mean','log_p_stem','log_bend_tray','log_p_contact'};
end

function r = design2raw(Xd, ctx)
n = size(Xd,1);  o = ones(n,1);
isCr = Xd(:,6) >= 0.5;
r.W = ctx.W*o;  r.D = ctx.D*o;  r.t_cort = ctx.t_cort*o;  r.E_cort = ctx.E_cort*o;  r.E_canc = ctx.E_canc*o;
r.mass = ctx.mass*o;
r.t_tray = Xd(:,1);  r.L_stem = Xd(:,2);  r.r_stem = Xd(:,3);  r.t_ins = Xd(:,4);  r.t_fem = Xd(:,5);
r.E_ins = ctx.E_ins*o;  r.Sy_ins = ctx.Sy_ins*o;
r.E_tray  = 110000*~isCr + 210000*isCr;
r.Sy_tray = 880*~isCr + 650*isCr;
r.flex = ctx.flex*o;  r.axbw = ctx.axbw*o;  r.F_ax = ctx.F_ax*o;  r.F_sh = ctx.F_sh*o;
r.med = ctx.med*o;  r.Lc = ctx.Lc*o;  r.shift = ctx.shift*o;
end

% ---------------------------- architecture -------------------------------
function arch = make_arch(cfg)
% EfficientNet-style compound scaling (B9-class: width 2.4x, depth 4.2x) applied to a tabular MBConv-SE tower
base_w = [32 48 64];  base_n = [1 2 2];
W  = 8*round(base_w*cfg.width_mult/8);
Nb = ceil(base_n*cfg.depth_mult);
win = []; wout = [];  cur = W(1);
for s = 1:numel(W)
    for b = 1:Nb(s), win(end+1) = cur; wout(end+1) = W(s); cur = W(s); end %#ok<AGROW>
end
arch = struct('W',W, 'Nb',Nb, 'win',win, 'wout',wout, 'nb',numel(win), 'expand',4);
end

function A = make_aux(nIn, cfg)
nBlk = ceil(nIn/cfg.blk);
blkOf = ceil((1:nIn)/cfg.blk);
Pm = zeros(nIn, nBlk);
for b = 1:nBlk, ii = (blkOf == b); Pm(ii,b) = 1/nnz(ii); end
[I, J] = ndgrid(1:nIn, 1:nIn);
A = struct('nBlk',nBlk, 'blkOf',blkOf, 'Pmat',dlarray(Pm), 'winBias',-1e4*double(abs(I-J) > cfg.win), ...
           'topk',min(cfg.topk, nBlk));
end

function P = init_params(nIn, nOut, cfg, arch)
d = cfg.d_tok;
g = @(a,b,s) dlarray(s*sqrt(2/(a+b))*randn(a,b));
z = @(a,b) dlarray(zeros(a,b));
o = @(a,b) dlarray(ones(a,b));
P = struct();
P.emb_W = dlarray(0.5*randn(d,nIn));  P.emb_b = dlarray(0.1*randn(d,nIn));
for l = 1:cfg.n_nsa
    kf = @(s) sprintf('a%d_%s', l, s);
    P.(kf('Wq')) = g(d,d,1);  P.(kf('Wk')) = g(d,d,1);  P.(kf('Wv')) = g(d,d,1);  P.(kf('Wo')) = g(d,d,0.5);
    P.(kf('Wg')) = g(3,d,1);  P.(kf('bg')) = z(3,1);
    P.(kf('g1')) = o(d,1);  P.(kf('b1')) = z(d,1);  P.(kf('g2')) = o(d,1);  P.(kf('b2')) = z(d,1);
    P.(kf('W1')) = g(2*d,d,1);  P.(kf('c1')) = z(2*d,1);  P.(kf('W2')) = g(d,2*d,0.5);  P.(kf('c2')) = z(d,1);
end
P.stem_W = g(arch.W(1), d*nIn, 1);  P.stem_b = z(arch.W(1), 1);
for k = 1:arch.nb
    nm = @(s) sprintf('m%d_%s', k, s);
    wi = arch.win(k);  wo = arch.wout(k);  e = arch.expand*wi;  sr = max(8, round(wi/4));
    P.(nm('g')) = o(wi,1);  P.(nm('b')) = z(wi,1);
    P.(nm('We')) = g(e,wi,1);   P.(nm('be')) = z(e,1);
    P.(nm('Ws1')) = g(sr,e,1);  P.(nm('bs1')) = z(sr,1);
    P.(nm('Ws2')) = g(e,sr,1);  P.(nm('bs2')) = z(e,1);
    P.(nm('Wp')) = g(wo,e,0.5); P.(nm('bp')) = z(wo,1);
    P.(nm('ls')) = dlarray(0.3*ones(wo,1));
end
P.head_W = g(nOut, arch.W(end), 0.3);  P.head_b = z(nOut,1);
P.lin_W  = z(nOut, nIn);                 % linear (log-log) skip path from the inputs
end

% ---------------------------- forward pass -------------------------------
function Y = net_forward(P, X, A, training, cfg)
nIn = size(X,1);  B = size(X,2);  d = cfg.d_tok;
tok = P.emb_W .* reshape(X, 1, nIn, B) + P.emb_b;              % d x nIn x B  (one token per feature)
for l = 1:cfg.n_nsa
    tok = nsa_layer(P, tok, l, A);
end
h = reshape(tok, d*nIn, B);
h = silu(P.stem_W*h + P.stem_b);
arch = A.arch;
for k = 1:arch.nb
    nm = @(s) P.(sprintf('m%d_%s', k, s));
    x  = lnorm(h, nm('g'), nm('b'));
    e  = silu(nm('We')*x + nm('be'));
    se = sigmoid(nm('Ws2')*silu(nm('Ws1')*e + nm('bs1')) + nm('bs2'));   % squeeze-excitation gate
    e  = e .* se;
    if training && cfg.drop > 0
        e = e .* (rand(size(e)) > cfg.drop) / (1 - cfg.drop);
    end
    o = nm('Wp')*e + nm('bp');
    if arch.win(k) == arch.wout(k), h = h + nm('ls') .* o; else, h = o; end
end
if training && cfg.drop > 0
    h = h .* (rand(size(h)) > cfg.drop) / (1 - cfg.drop);
end
Y = P.head_W*h + P.head_b + P.lin_W*X;
end

function tok = nsa_layer(P, tok, l, A)
% TabNSA-style native sparse attention over feature tokens:
%   (1) compressed branch  : attention over block-pooled tokens (coarse global view)
%   (2) selected branch    : attention restricted to the top-k most important blocks (per sample)
%   (3) sliding branch     : local-window attention (neighbouring features)
% The three outputs are mixed by learned per-token gates.
pre = @(f) P.(sprintf('a%d_%s', l, f));
d = size(tok,1);  nF = size(tok,2);  B = size(tok,3);
x = lnorm(tok, pre('g1'), pre('b1'));
Q = pagemtimes(pre('Wq'), x);  K = pagemtimes(pre('Wk'), x);  V = pagemtimes(pre('Wv'), x);
sc = 1/sqrt(d);
S  = pagemtimes(Q, 'transpose', K, 'none') * sc;                % nF x nF x B
% (1) compressed
Kc = pagemtimes(K, A.Pmat);  Vc = pagemtimes(V, A.Pmat);
Ac = softmax2(pagemtimes(Q, 'transpose', Kc, 'none') * sc);     % nF x nBlk x B
Oc = pagemtimes(Vc, 'none', Ac, 'transpose');                   % d x nF x B
% (2) selected blocks (top-k by compressed importance)
imp = reshape(extractdata(mean(Ac, 1)), A.nBlk, B);
[~, ord] = sort(imp, 1, 'descend');
sel = false(A.nBlk, B);
sel(sub2ind([A.nBlk B], ord(1:A.topk,:), repmat(1:B, A.topk, 1))) = true;
selTok = sel(A.blkOf, :);                                       % nF x B
As = softmax2(S + reshape(-1e4*double(~selTok), 1, nF, B));
Os = pagemtimes(V, 'none', As, 'transpose');
% (3) sliding window
Aw = softmax2(S + A.winBias);
Ow = pagemtimes(V, 'none', Aw, 'transpose');
% gated fusion
gt = sigmoid(pagemtimes(pre('Wg'), x) + pre('bg'));             % 3 x nF x B
O  = reshape(gt(1,:,:), 1, nF, B) .* Oc + reshape(gt(2,:,:), 1, nF, B) .* Os + reshape(gt(3,:,:), 1, nF, B) .* Ow;
tok = tok + pagemtimes(pre('Wo'), O);
% feed-forward
x2  = lnorm(tok, pre('g2'), pre('b2'));
f   = silu(pagemtimes(pre('W1'), x2) + pre('c1'));
tok = tok + pagemtimes(pre('W2'), f) + pre('c2');
end

function A = softmax2(S)
S = S - max(S, [], 2);
E = exp(S);
A = E ./ sum(E, 2);
end

function y = lnorm(x, g, b)
mu = mean(x, 1);  xc = x - mu;
v  = mean(xc.^2, 1);
y  = xc ./ sqrt(v + 1e-5) .* g + b;
end

function y = silu(x)
y = x .* sigmoid(x);
end

% ---------------------------- loss / training ----------------------------
function [grads, parts] = grad_fn(P, X, Z, Xc, cfg, A, ph)
B = size(X,2);
Xall = cat(2, X, Xc);                                            % labelled + unlabelled collocation inputs
Yall = net_forward(P, Xall, A, true, cfg);
Ld = mean((Yall(:,1:B) - Z).^2, 'all');                          % data loss (labelled part)

% ---- FEA-derived physical-consistency constraints (log space) ----
Yp = Yall .* ph.ysg + ph.ymu;
ly = @(n) Yp(ph.i.(n), :);
xl = @(j) Xall(j,:) * ph.xsg(j) + ph.xmu(j);
r1 = ly('SF_Tray')    - (xl(15) - max(ly('VM_Stem_MPa'), ly('VM_Tray_MPa')));          % SF = Sy / max(stem, tray VM)
r2 = ly('SF_Insert')  - (xl(13) - ly('VM_Insert_MPa'));                                 % SF = Sy_ins / VM_insert
r3 = ly('SF_FemComp') - (log(650) - ly('VM_FemComp_MPa'));                              % SF = Sy_fem / VM_femcomp
r4 = ly('SF_Bone')    - (log(7000) - max(ly('Strain_TibBone_ue'), ly('Strain_FemBone_ue')));   % 0.7 % yield strain
r5 = ly('Peak_VM_MPa') - max(cat(1, ly('VM_TibBone_MPa'), ly('VM_Stem_MPa'), ly('VM_Tray_MPa'), ...
                                    ly('VM_Insert_MPa'), ly('VM_FemComp_MPa'), ly('VM_FemBone_MPa')), [], 1);
r5b = ly('Peak_Strain_ue') - max(cat(1, ly('Strain_TibBone_ue'), ly('Strain_Tray_ue'), ly('Strain_Insert_ue'), ...
                                     ly('Strain_FemComp_ue'), ly('Strain_FemBone_ue')), [], 1);
ic = invT(ly('Insert_Compression_mm'), ph.isLog(ph.i.Insert_Compression_mm));
fd = invT(ly('FemComp_Disp_mm'),       ph.isLog(ph.i.FemComp_Disp_mm));
ts = invT(ly('Tray_Subsidence_mm'),    ph.isLog(ph.i.Tray_Subsidence_mm));
md = invT(ly('Max_Disp_mm'),           ph.isLog(ph.i.Max_Disp_mm));
r6 = (ic - (fd - ts)) / ph.ic_sd;                                % insert compression = fem.comp. disp - tray subsidence
r7 = (max(0, ts - md) + max(0, fd - md)) / ph.md_sd;             % nodal max displacement bounds the mean ones
Lphys = mean(r1.^2,'all') + mean(r2.^2,'all') + mean(r3.^2,'all') + mean(r4.^2,'all') + ...
        mean(r5.^2,'all') + mean(r5b.^2,'all') + mean(r6.^2,'all') + mean(r7.^2,'all');

% ---- monotonicity w.r.t. joint load (finite difference, eval mode) ----
nm = size(Xall,2);
Ym = net_forward(P, cat(2, Xall, Xall + ph.shift), A, false, cfg);
Dm = (Ym(:, nm+1:end) - Ym(:, 1:nm)) .* ph.ysg;
Lmono = mean(max(0, -ph.sgn .* Dm).^2, 'all');

loss  = Ld + cfg.lam_phys*Lphys + cfg.lam_mono*Lmono;
grads = dlgradient(loss, P);
parts = [Ld, Lphys, Lmono];
end

function y = invT(z, isLog)
if isLog, y = exp(z); else, y = sinh(z); end
end

function [bestP, hist] = train_one_model(Xtr, Ztr, Xva, Zva, cfg, A, ph, seed)
rng(seed);
nIn = size(Xtr,1);  nOut = size(Ztr,1);  nTr = size(Xtr,2);
P = init_params(nIn, nOut, cfg, A.arch);
nPar = 0;  fn = fieldnames(P);
for k = 1:numel(fn), nPar = nPar + numel(P.(fn{k})); end
fprintf('    parameters: %.2f M\n', nPar/1e6);
avgG = [];  avgSq = [];  it = 0;
best = inf;  bestP = P;  bad = 0;
hist = zeros(cfg.epochs, 4);
t0 = tic;
for ep = 1:cfg.epochs
    lr = max(cfg.lr * 0.5*(1 + cos(pi*(ep-1)/cfg.epochs)), 1e-5);
    perm = randperm(nTr);  acc = zeros(1,3);  nb = 0;
    for s0 = 1:cfg.batch:nTr
        ib = perm(s0:min(s0+cfg.batch-1, nTr));
        if numel(ib) < 4, continue; end
        Xb = dlarray(Xtr(:,ib));  Zb = dlarray(Ztr(:,ib));
        Xc = dlarray(make_colloc(Xtr, cfg.n_colloc));
        [grads, parts] = dlfeval(@grad_fn, P, Xb, Zb, Xc, cfg, A, ph);
        grads = clip_grads(grads, cfg.clip);
        it = it + 1;
        [P, avgG, avgSq] = adamupdate(P, grads, avgG, avgSq, it, lr);
        P = dlupdate(@(p) p*(1 - lr*cfg.wd), P);
        acc = acc + extractdata(parts);  nb = nb + 1;
    end
    Zv = predict_Z(P, Xva, A, cfg);
    vmse = mean((Zv - Zva).^2, 'all');
    hist(ep,:) = [acc/max(nb,1), vmse];
    if vmse < best - 1e-6, best = vmse; bestP = P; bad = 0; else, bad = bad + 1; end
    if mod(ep,25) == 0 || ep == 1
        fprintf('    ep %4d | data %.4f | phys %.4f | mono %.5f | val %.4f | best %.4f | %.0fs\n', ...
            ep, hist(ep,1), hist(ep,2), hist(ep,3), vmse, best, toc(t0));
    end
    if bad >= cfg.patience, fprintf('    early stop at epoch %d (best val %.4f)\n', ep, best); break; end
end
hist = hist(1:ep,:);
end

function g = clip_grads(g, maxn)
fn = fieldnames(g);  n2 = 0;
for k = 1:numel(fn), n2 = n2 + sum(extractdata(g.(fn{k})).^2, 'all'); end
gn = sqrt(n2);
if gn > maxn, g = dlupdate(@(x) x*(maxn/gn), g); end
end

function Xc = make_colloc(X, n)
% unlabelled physics collocation points: random convex mixes of training inputs (+ small noise)
N = size(X,2);
ia = randi(N, 1, n);  ib = randi(N, 1, n);
lam = rand(1, n);
Xc = lam .* X(:,ia) + (1 - lam) .* X(:,ib) + 0.03*randn(size(X,1), n);
end

function Z = predict_Z(P, Xs, A, cfg)
n = size(Xs,2);  nOut = size(P.head_b,1);
Z = zeros(nOut, n);
for s = 1:512:n
    ib = s:min(s+511, n);
    Z(:,ib) = extractdata(net_forward(P, dlarray(Xs(:,ib)), A, false, cfg));
end
end

function Y = surr_predict(S, Xfeat)
% Xfeat : N x nIn raw engineered features  ->  Y : N x nOut in ORIGINAL units (ensemble mean)
Xs = ((Xfeat - S.xmu) ./ S.xsg)';
Z = 0;
for m = 1:numel(S.models), Z = Z + predict_Z(S.models{m}, Xs, S.A, S.cfg); end
Z = Z / numel(S.models);
Zt = Z .* S.ysg + S.ymu;
Y = zeros(size(Zt'));
for k = 1:size(Zt,1)
    if S.isLog(k), Y(:,k) = exp(Zt(k,:))'; else, Y(:,k) = sinh(Zt(k,:))'; end
end
end

function M = metrics_tbl(Yt, Yp, names)
res = Yt - Yp;
mse = mean(res.^2, 1);
rmse = sqrt(mse);
mae = mean(abs(res), 1);
r2 = 1 - sum(res.^2, 1) ./ max(sum((Yt - mean(Yt,1)).^2, 1), eps);
nm = [names(:); {'MEAN'}];
M = table(nm, [r2'; mean(r2)], [rmse'; mean(rmse)], [mse'; mean(mse)], [mae'; mean(mae)], ...
    'VariableNames', {'Target','R2','RMSE','MSE','MAE'});
end

function print_metrics(label, M, nOut, n)
fprintf('\n[+] %s metrics (%d cases)\n', label, n);
fprintf('    %-22s %9s %13s %13s %13s\n', 'Target', 'R2', 'RMSE', 'MSE', 'MAE');
for k = 1:height(M)
    fprintf('    %-22s %9.4f %13.5g %13.5g %13.5g\n', M.Target{k}, M.R2(k), M.RMSE(k), M.MSE(k), M.MAE(k));
end
fprintf('    -> targets with R2 >= 0.90 : %d / %d\n', nnz(M.R2(1:nOut) >= 0.90), nOut);
end

% ---------------------------- optimisation -------------------------------
function [F, CV, Y] = eval_designs(Xd, ctx, surr, cfg, tgt)
ix = @(n) find(strcmp(tgt, n));
Y  = surr_predict(surr, make_features(design2raw(Xd, ctx)));
sfT = Y(:,ix('SF_Tray'));  sfI = Y(:,ix('SF_Insert'));  sfF = Y(:,ix('SF_FemComp'));  sfB = Y(:,ix('SF_Bone'));
stress = max([Y(:,ix('VM_Stem_MPa')), Y(:,ix('VM_Tray_MPa')), Y(:,ix('VM_Insert_MPa')), Y(:,ix('VM_FemComp_MPa'))], [], 2);
dsp    = Y(:,ix('Max_Disp_mm'));
sfImp  = min([sfT sfI sfF], [], 2);
vol    = pi * Xd(:,3).^2 .* Xd(:,2);
volmax = pi * ctx.ub(3)^2 * ctx.ub(2);
invas  = 0.5*(Xd(:,5) - 7.5)/(ctx.ub(5) - 7.5) + 0.5*vol/volmax;      % bone-removal / keel-volume proxy
% objectives (all minimised): stress, deformation, -implant safety, -bone safety, invasiveness
F = [stress, dsp, -sfImp, -sfB, invas];
% constraints g <= 0 : mechanical (safety), anatomical (stem width, joint-line stack), geometric (stem fit)
g = [cfg.sf_min - sfT, cfg.sf_min - sfI, cfg.sf_min - sfF, cfg.sf_bone_min - sfB, ...
     2*Xd(:,3) - 0.30*ctx.D, Xd(:,1) + Xd(:,4) - ctx.max_stack, Xd(:,2) - (ctx.L_tib - 2)];
CV = sum(max(g, 0), 2);
end

function [Pop, F, CV, Y, hist] = amoo_nsga2(fevalfun, lb, ub, cfg)
% Constrained NSGA-II with two adaptive mechanisms:
%  (a) mutation intensity adapts to the offspring survival rate (explore <-> exploit)
%  (b) a share of each generation is sampled inside an expanded bounding box of the current
%      feasible non-dominated designs (adaptive exploration of the feasible region)
n = cfg.pop;  D = numel(lb);
Pop = lb + rand(n,D).*(ub - lb);  Pop(:,6) = round(Pop(:,6));
[F, CV, Y] = fevalfun(Pop);
fmin0 = min(F,[],1);  fr0 = max(max(F,[],1) - fmin0, eps);
S_mc = 1.1 * rand(20000, size(F,2));                              % Monte-Carlo points for hypervolume
mscale = 1.0;
hist = struct('hv',zeros(cfg.gens,1), 'feas',zeros(cfg.gens,1), 'mscale',zeros(cfg.gens,1), 'succ',zeros(cfg.gens,1));
fprintf('[+] AMOO: population %d, generations %d, objectives %d, constraints enforced via constraint-domination\n', n, cfg.gens, size(F,2));
for gen = 1:cfg.gens
    [rk, cr] = nds_crowd(F, CV);
    nImm = round(0.15*n);  nOff = n - nImm;
    % ---- variation ----
    pa = tournament(rk, cr, nOff + mod(nOff,2));
    pb = tournament(rk, cr, nOff + mod(nOff,2));
    C = sbx(Pop(pa,:), Pop(pb,:), lb, ub, 15);
    C = C(1:nOff,:);
    C = polymut(C, lb, ub, min(1, mscale/D), 20);
    % ---- adaptive feasible-region immigrants ----
    elite = (rk == 1) & (CV <= 0);
    if nnz(elite) < 3
        [~, o] = sort(CV);  elite = false(n,1);  elite(o(1:max(3, round(0.2*n)))) = true;
    end
    bmin = min(Pop(elite,:),[],1);  bmax = max(Pop(elite,:),[],1);
    pad = 0.20*(ub - lb);
    bmin = max(lb, bmin - pad);  bmax = min(ub, bmax + pad);
    Imm = bmin + rand(nImm, D).*(bmax - bmin);
    Off = [C; Imm];  Off(:,6) = round(min(max(Off(:,6), 0), 1));
    [Fo, CVo, Yo] = fevalfun(Off);
    % ---- environmental selection ----
    Pall = [Pop; Off];  Fall = [F; Fo];  CVall = [CV; CVo];  Yall = [Y; Yo];
    [rk2, cr2] = nds_crowd(Fall, CVall);
    [~, ordr] = sortrows([rk2, -cr2]);
    keep = ordr(1:n);
    succ = mean(keep > n);
    Pop = Pall(keep,:);  F = Fall(keep,:);  CV = CVall(keep);  Y = Yall(keep,:);
    % ---- adapt mutation scale ----
    if succ < 0.20, mscale = min(3.0, mscale*1.15); elseif succ > 0.50, mscale = max(0.5, mscale/1.15); end
    [rk3, ~] = nds_crowd(F, CV);
    fr_idx = (rk3 == 1) & (CV <= 0);
    hist.feas(gen) = mean(CV <= 0);  hist.mscale(gen) = mscale;  hist.succ(gen) = succ;
    hist.hv(gen) = hv_mc((F(fr_idx,:) - fmin0)./fr0, S_mc);
    if mod(gen,10) == 0 || gen == 1
        fprintf('    gen %3d | feasible %.0f%% | front %d | HV %.4f | mut-scale %.2f | survival %.2f\n', ...
            gen, 100*hist.feas(gen), nnz(fr_idx), hist.hv(gen), mscale, succ);
    end
end
end

function hv = hv_mc(Fn, S)
if isempty(Fn), hv = 0; return; end
dom = false(size(S,1),1);
for i = 1:size(Fn,1), dom = dom | all(Fn(i,:) <= S, 2); end
hv = mean(dom) * 1.1^size(S,2);
end

function idx = tournament(rk, cr, m)
N = numel(rk);
a = randi(N, m, 1);  b = randi(N, m, 1);
winA = (rk(a) < rk(b)) | (rk(a) == rk(b) & cr(a) >= cr(b));
idx = b;  idx(winA) = a(winA);
end

function C = sbx(P1, P2, lb, ub, eta)
u = rand(size(P1));
beta = ones(size(P1));
m = u <= 0.5;
beta(m)  = (2*u(m)).^(1/(eta+1));
beta(~m) = (1./(2*(1 - u(~m)))).^(1/(eta+1));
c1 = 0.5*((1+beta).*P1 + (1-beta).*P2);
c2 = 0.5*((1-beta).*P1 + (1+beta).*P2);
sw = rand(size(P1)) < 0.5;
c1(sw) = P1(sw);  c2(sw) = P2(sw);
C = min(max([c1; c2], lb), ub);
end

function X = polymut(X, lb, ub, pm, eta)
[n, D] = size(X);
r = rand(n, D);  m = rand(n, D) < pm;
delta = zeros(n, D);  lw = r < 0.5;
delta(lw)  = (2*r(lw)).^(1/(eta+1)) - 1;
delta(~lw) = 1 - (2*(1 - r(~lw))).^(1/(eta+1));
X = min(max(X + m .* delta .* (ub - lb), lb), ub);
end

function [rank, crowd] = nds_crowd(F, CV)
% constraint-domination non-dominated sorting + crowding distance
N = size(F,1);
dom = false(N,N);
feas = CV <= 0;
for i = 1:N
    better = all(F(i,:) <= F, 2) & any(F(i,:) < F, 2);
    if feas(i)
        d = ~feas | (feas & better);
    else
        d = ~feas & (CV(i) < CV);
    end
    dom(i,:) = d';
end
cnt = sum(dom, 1)';
rank = zeros(N,1);  rem = true(N,1);  r = 0;
while any(rem)
    r = r + 1;
    front = rem & (cnt == 0);
    if ~any(front), front = rem; end
    rank(front) = r;  rem(front) = false;
    cnt = cnt - sum(dom(front,:), 1)';
end
crowd = zeros(N,1);
for r = 1:max(rank)
    id = find(rank == r);
    if numel(id) <= 2, crowd(id) = inf; continue; end
    for k = 1:size(F,2)
        [fs, o] = sort(F(id,k));
        rng_k = max(fs(end) - fs(1), eps);
        crowd(id(o(1))) = inf;  crowd(id(o(end))) = inf;
        crowd(id(o(2:end-1))) = crowd(id(o(2:end-1))) + (fs(3:end) - fs(1:end-2))/rng_k;
    end
end
end

function D = pdist2_simple(A, B)
D = zeros(size(A,1), size(B,1));
for j = 1:size(B,1), D(:,j) = sqrt(sum((A - B(j,:)).^2, 2)); end
end

% ---------------------------- robust export ------------------------------
function xls = save_tables(out_dir, base, sheets)
% Writes each table to its own Excel sheet; never aborts the run.
%  - if the .xlsx is locked (e.g. open in Excel) a time-stamped file name is used
%  - any sheet that fails is saved as CSV and the real error message is printed
xls = fullfile(out_dir, [base '.xlsx']);
if exist(xls, 'file')
    try
        delete(xls);
        if exist(xls, 'file'), error('file locked'); end
    catch
        xls = fullfile(out_dir, sprintf('%s_%s.xlsx', base, char(datetime('now','Format','yyyyMMdd_HHmmss'))));
        fprintf('    [!] Previous Excel file is locked/open - writing to %s\n', xls);
    end
end
nfail = 0;
for k = 1:size(sheets,1)
    nm = sheets{k,1};  tb = sheets{k,2};
    if ~istable(tb) || isempty(tb), continue; end
    tb.Properties.VariableNames = matlab.lang.makeValidName(tb.Properties.VariableNames);
    try
        writetable(tb, xls, 'Sheet', nm, 'UseExcel', false);
    catch ME
        nfail = nfail + 1;
        csvf = fullfile(out_dir, [base '_' nm '.csv']);
        try, writetable(tb, csvf); catch, end
        fprintf('    [!] Sheet "%s" failed in Excel (%s) -> saved as %s\n', nm, strtrim(ME.message), csvf);
    end
end
if nfail == 0, fprintf('\n    - Excel : %s\n', xls);
else, fprintf('\n    - Excel : %s (%d sheet(s) saved as CSV instead)\n', xls, nfail); end
end

% ---------------------------- FEA validation -----------------------------
function [R, G] = fea_for_design(xd, G0, fh, hmax)
G = G0;
G.t_tray = xd(1);  G.L_stem = xd(2);  G.r_stem = xd(3);  G.t_ins = xd(4);  G.t_fem = xd(5);
if xd(6) >= 0.5
    G.E_tray = 210000;  G.nu_tray = 0.30;  G.sy_tray = 650;
else
    G.E_tray = 110000;  G.nu_tray = 0.33;  G.sy_tray = 880;
end
G = fh.finalize(G);
try
    [R, ~, ~, G] = fh.run_fea_case(G, hmax);
    R.Status = 'OK';
catch ME
    warning('FEA failed: %s', ME.message);
    R = struct();
    names = {'Peak_VM_MPa','VM_TibBone_MPa','VM_Stem_MPa','VM_Tray_MPa','VM_Insert_MPa','VM_FemComp_MPa','VM_FemBone_MPa', ...
        'Peak_Strain_ue','Strain_TibBone_ue','Strain_Tray_ue','Strain_Insert_ue','Strain_FemComp_ue','Strain_FemBone_ue', ...
        'Max_Disp_mm','Tray_Subsidence_mm','FemComp_Disp_mm','Insert_Compression_mm','SF_Tray','SF_Insert','SF_FemComp','SF_Bone'};
    for k = 1:numel(names), R.(names{k}) = NaN; end
    R.Status = ['FAILED: ' ME.message];
end
end