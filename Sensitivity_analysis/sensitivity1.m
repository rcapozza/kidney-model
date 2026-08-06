% SENSITIVITY1  One-parameter-at-a-time (OAT) sensitivity analysis of the
% fisher_evol tip/stalk/collision model around a reference (previously
% optimized) parameter set. For each of the 8 model parameters
% (as, dif, pc, gcv, bs, rb, rg, scale), only that parameter is varied
% over a +/-25% range while the others are held at their reference
% value; the resulting change in the simulated density is used to
% compute a normalized "sensitivity" and a plain relative error, which
% are plotted against the parameter's relative deviation from reference.
%
% Dependencies: preprocess_image.m, fisher_evol.m, s3.png (all must be on
% the MATLAB path / in the same folder).
clear all
close all
% Number of time steps for the simulation
f = 1300;
nbox = 20;
% Preprocess the image and get the initial density matrices
[A, S, C] = preprocess_image('s3.png', nbox, f);

% Reference (previously optimized) parameters
as_ref = 1200;         % saturation threshold (inactivation)
dif_ref = 0.00032;     % diffusion coefficient
pc_ref = 0.28305;      % fraction of saturated tips that collide
gcv_ref = 0.014;       % stall-condition threshold
bs_ref = 1171.476;     % branching-saturation gate centre (already pre-multiplied by as_ref; see fisher_evol.m)
rb_ref = 0.333;        % branching rate
rg_ref = 0.75;         % growth rate
scale_ref = 407.9998;  % scale factor (controls the sigmoid gate width b = as/scale)

[A_ref, S_ref, C_ref] = fisher_evol(f, as_ref, dif_ref, pc_ref, gcv_ref, bs_ref, rb_ref, rg_ref, scale_ref, A, S, C);
    X = linspace(0,1,nbox);% returns a row vector of I evenly spaced points between 0 and 1
    Y = linspace(0,1,nbox);% in our scale this is about 1->1mm
    figure(4)
    surf(X,Y,A_ref(:,:)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar
    figure(5)
    surf(X,Y,S_ref(:,:)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar
    figure(6)
    surf(X,Y,C_ref(:,:)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar
% ref_concentration is a nbox x (3*nbox) matrix formed by horizontally
% concatenating [A_ref | S_ref | C_ref]; conc_mod is its overall
% magnitude (Frobenius norm), used to normalize the errors below.
ref_concentration = cell2mat({A_ref, S_ref, C_ref});
conc_mod = sqrt(sum(sum((ref_concentration).^2)));
ns = 24; % size of the sweep grid for each parameter - must be even, otherwise one sample would land exactly on the reference value (das=0), causing a divide-by-zero in the sensitivity below
Snt = zeros(ns,ns,ns);
Sens = zeros(ns*ns*ns,1);
s_as = zeros(ns,1);
s_dif = zeros(ns,1);
s_pc = zeros(ns,1);
s_gcv = zeros(ns,1);
s_bs = zeros(ns,1);
s_rb = zeros(ns,1);
s_rg = zeros(ns,1);
s_scale = zeros(ns,1);
err_as = zeros(ns,1);
err_dif = zeros(ns,1);
err_pc = zeros(ns,1);
err_gcv = zeros(ns,1);
err_bs = zeros(ns,1);
err_rb = zeros(ns,1);
err_rg = zeros(ns,1);
err_scale = zeros(ns,1);
fra = 4.; % each parameter is varied by +/-(1/fra) = +/-25% around its reference value
as_values = linspace(as_ref - as_ref*(1./fra), as_ref + as_ref*(1./fra), ns);
pc_values = linspace(pc_ref - pc_ref*(1./fra), pc_ref + pc_ref*(1./fra), ns);
dif_values = linspace(dif_ref - dif_ref*(1./fra), dif_ref + dif_ref*(1./fra), ns);
gcv_values = linspace(gcv_ref - gcv_ref*(1./fra), gcv_ref + gcv_ref*(1./fra), ns);
bs_values = linspace(bs_ref - bs_ref*(1./fra), bs_ref + bs_ref*(1./fra), ns);
rb_values = linspace(rb_ref - rb_ref*(1./fra), rb_ref + rb_ref*(1./fra), ns);
rg_values = linspace(rg_ref - rg_ref*(1./fra), rg_ref + rg_ref*(1./fra), ns);
scale_values = linspace(scale_ref - scale_ref*(1./fra), scale_ref + scale_ref*(1./fra), ns);

% --- Sensitivity sweeps -------------------------------------------------
% For each parameter below, only that parameter is varied (all others
% stay at their reference value) and the model is re-simulated. Two
% quantities are computed and stored for each sweep point:
%   err         = ||ref_concentration - sim_concentration|| / ||ref_concentration||
%                 (plain normalized error between the perturbed and
%                 reference run)
%   sensitivity = err * (param_ref / |param - param_ref|)
%                 (the error rescaled by the reference-to-perturbation
%                 ratio, so it approximates a normalized derivative of
%                 the error with respect to the parameter)
% This one-parameter-at-a-time (OAT) sweep is repeated below for each of
% the 8 model parameters. Note err and sensitivity each recompute the
% same sqrt(sum(sum(...))) term independently, so the underlying error
% norm is effectively calculated twice per sweep point.
index = 0; % not used by the OAT sweeps below (left over from an earlier full-grid approach)
nas = 0;
ndif = 0;
npc = 0;
for as = as_values
            nas = nas + 1;
            das = abs(as - as_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as,dif_ref,pc_ref,gcv_ref,bs_ref,rb_ref,rg_ref,scale_ref, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = (sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod)*(as_ref/das);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            s_as(nas) = squarediff
            err_as(nas) = err
end

for dif = dif_values
            ndif = ndif + 1;
            ddif = abs(dif - dif_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as_ref,dif,pc_ref,gcv_ref,bs_ref,rb_ref,rg_ref,scale_ref, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod*(dif_ref/ddif);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            dif_ref/ddif
            s_dif(ndif) = squarediff
            err_dif(ndif) = err
end

for pc = pc_values
            npc = npc + 1;
            dpc = abs(pc - pc_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as_ref,dif_ref,pc,gcv_ref,bs_ref,rb_ref,rg_ref,scale_ref, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod*(pc_ref/dpc);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            s_pc(npc) = squarediff
            err_pc(npc) = err
end

ni = 0;
for gcv = gcv_values
            ni = ni + 1;
            dgcv = abs(gcv - gcv_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as_ref,dif_ref,pc_ref,gcv,bs_ref,rb_ref,rg_ref,scale_ref, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = (sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod)*(gcv_ref/dgcv);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            s_gcv(ni) = squarediff
            err_gcv(ni) = err
end

% Sensitivity to bs
ni = 0;
for bs = bs_values
            ni = ni + 1;
            dbs = abs(bs - bs_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as_ref,dif_ref,pc_ref,gcv_ref,bs,rb_ref,rg_ref,scale_ref, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = (sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod)*(bs_ref/dbs);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            s_bs(ni) = squarediff
            err_bs(ni) = err
end

% Sensitivity to rb
ni = 0;
for rb = rb_values
            ni = ni + 1;
            drb = abs(rb - rb_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as_ref,dif_ref,pc_ref,gcv_ref,bs_ref,rb,rg_ref,scale_ref, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = (sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod)*(rb_ref/drb);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            s_rb(ni) = squarediff
            err_rb(ni) = err
end

% Sensitivity to rg
ni = 0;
for rg = rg_values
            ni = ni + 1;
            drg = abs(rg - rg_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as_ref,dif_ref,pc_ref,gcv_ref,bs_ref,rb_ref,rg,scale_ref, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = (sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod)*(rg_ref/drg);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            s_rg(ni) = squarediff
            err_rg(ni) = err
end

% Sensitivity to scale
ni = 0;
for scale = scale_values
            ni = ni + 1;
            dscale = abs(scale - scale_ref);
            [A_sim, S_sim, C_sim] = fisher_evol(f,as_ref,dif_ref,pc_ref,gcv_ref,bs_ref,rb_ref,rg_ref,scale, A, S, C);
            sim_concentration = cell2mat({A_sim, S_sim, C_sim});
            squarediff = (sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod)*(scale_ref/dscale);
            err = sqrt(sum(sum(((ref_concentration - sim_concentration)).^2)))/conc_mod;
            s_scale(ni) = squarediff
            err_scale(ni) = err
end

%% Plot sensitivity and error for each parameter
% Figure 7: as, dif, pc, gcv
figure(7);
subplot(2,2,1);
yyaxis left
plot((as_values - as_ref) / as_ref, s_as, '-o', 'LineWidth', 1.8);
hold on
yyaxis right
plot((as_values - as_ref) / as_ref, err_as, '-o', 'LineWidth', 1.8);
hold on
xlabel('\Deltaa_{c}/a_{c}^{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{a_{c}}', 'Err_{a_{c}}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;      % tick label size
ax.FontName = 'Arial'; % font type for ticks
ax.LineWidth = 1.2;    % axis line thickness

subplot(2,2,2);
yyaxis left
plot((dif_values - dif_ref) / dif_ref, s_dif, '-o', 'LineWidth', 1.8);
hold on
yyaxis right
plot((dif_values - dif_ref) / dif_ref, err_dif, '-o', 'LineWidth', 1.8);
hold on
xlabel('\DeltaD/D^{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{D}', 'Err_{D}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;

subplot(2,2,3);
yyaxis left
plot((pc_values - pc_ref) / pc_ref, s_pc, '-o', 'LineWidth', 1.8);
hold on
yyaxis right
plot((pc_values - pc_ref) / pc_ref, err_pc, '-o', 'LineWidth', 1.8);
hold on
xlabel('\Deltar_i/r_i^{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{r_i}', 'Err_{r_i}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;

subplot(2,2,4);
yyaxis left
plot((gcv_values - gcv_ref) / gcv_ref, s_gcv, '-o', 'LineWidth', 1.8);
hold on
yyaxis right
plot((gcv_values - gcv_ref) / gcv_ref, err_gcv, '-o', 'LineWidth', 1.8);
hold on
xlabel('\Deltag_{thresh}/g_{thresh}^{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{g_{thresh}}', 'Err_{g_{thresh}}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;

% Figure 8: bs, rb, rg, scale
figure(8)
subplot(2,2,1);
yyaxis left
plot((bs_values - bs_ref) / bs_ref, s_bs, '-o', 'LineWidth', 1.8);
hold on
yyaxis right
plot((bs_values - bs_ref) / bs_ref, err_bs, '-o', 'LineWidth', 1.8);
hold on
xlabel('\Deltaa_{b}/a_{b}^{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{a_{b}}', 'Err_{a_{b}}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;

subplot(2,2,2);
yyaxis left
plot((rb_values - rb_ref) / rb_ref, s_rb, '-o', 'LineWidth', 1.8);
hold on
yyaxis right
plot((rb_values - rb_ref) / rb_ref, err_rb, '-o', 'LineWidth', 1.8);
hold on
xlabel('\Deltar_{b}/r_{b}^{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{r_{b}}', 'Err_{r_{b}}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;

subplot(2,2,3);
yyaxis left
plot((rg_values - rg_ref) / rg_ref, s_rg, '-o', 'LineWidth', 1.8);
hold on
yyaxis right
plot((rg_values - rg_ref) / rg_ref, err_rg, '-o', 'LineWidth', 1.8);
hold on
xlabel('\Deltar_{g}/r_{g}^{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{r_{g}}', 'Err_{r_{g}}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;

subplot(2,2,4);
hold on
yyaxis left
plot((scale_values - scale_ref) / scale_ref, s_scale, '-s', 'LineWidth', 1.8);
yyaxis right
plot((scale_values - scale_ref) / scale_ref, err_scale, '-s', 'LineWidth', 1.8);
xlabel('\Delta_{scale}/scale_{ref}', 'Interpreter', 'tex', 'FontSize', 16, 'FontName', 'Arial');
legend('S_{scale}', 'Err_{scale}', 'Location', 'best', 'FontSize', 16);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;
