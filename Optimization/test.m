clear all
close all

%% Load the target image and extract density grids
% s39.png is the "target" data set: the density pattern that the
% optimization below tries to reproduce by running the reaction-diffusion
% model implemented in fisher4.m.
image = imread('s39.png');

% Classify each pixel by color into the tracked populations:
%   red    -> active tips
%   black  -> stalks
%   yellow -> inactive/collided entities
% (green is also computed but not used further in this script)
black_channel = sum(image, 3) == 0; % black assumed to be RGB [0 0 0]
red_channel = image(:,:,1)-image(:,:,2)==255;
yellow_channel = image(:,:,2)-image(:,:,3)==255;
green_channel = image(:,:,2)-image(:,:,1)>50;
I=size(image)

nbox = 20; % number of boxes to divide the lateral field into (spatial grid resolution)
npix = fix(I(1)/nbox)% pixels per box; assumes the image is square, fix() keeps the integer part
Af = zeros(nbox,nbox,1); % active-tip density grid
Sf = zeros(nbox,nbox,1); % stalk density grid
Cf = zeros(nbox,nbox,1); % inactive/collided density grid
Exc = zeros(nbox,nbox,1);
X = linspace(0,1,nbox);% row vector of nbox evenly spaced points between 0 and 1
Y = linspace(0,1,nbox);% in this model's scale, 0->1 corresponds to about 0->1 mm

% Calculate the density in each box by counting matching pixels
for j=1:nbox
    nr = npix*(j-1);
    for i=1:nbox
        nc = npix*(i-1);
        for k = 1:npix
            for m = 1:npix
            Af(j,i,1) = Af(j,i,1) + red_channel(nr+k,nc+m);
            Sf(j,i,1) = Sf(j,i,1) + black_channel(nr+k,nc+m);
            Cf(j,i,1) = Cf(j,i,1) + yellow_channel(nr+k,nc+m);
            end
        end

    end

end

% Normalize pixel counts into densities: each tip corresponds to 953 pixels
dx = 1/nbox;
Af=Af/(953*dx*dx);
Sf=Sf/(953*dx*dx);
Cf=Cf/(953*dx*dx);

%% Plot the target densities (figures 1-3)
    figure(1)
    surf(X,Y,Af(:,:,1)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

    figure(2)
    surf(X,Y,Sf(:,:,1)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

    figure(3)
    surf(X,Y,Cf(:,:,1)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

%% Target concentration data used by the optimizer
% Only concentration_data_a_6s (active-tip density) is actually compared
% against the model output inside objective_fun; the other two are passed
% through but not currently used there.
concentration_data_s_6s = Sf;  % target stalk density
concentration_data_a_6s = Af;  % target active-tip density
concentration_data_c_6s = Cf;  % target inactive/collided density

%% Optimization setup
% Parameter order: [dif, as, pc, gcv, bs, rb, rg, scale]
%   dif   - diffusion coefficient
%   as    - saturation threshold (inactivation)
%   pc    - fraction of saturated tips that collide
%   gcv   - stall-condition threshold
%   bs    - branching-saturation fraction
%   rb    - branching rate
%   rg    - growth rate
%   scale - scale factor
initial_guesses = [0.00045,1000,0.5,0.014,0.98,1./(9./3.),1./(4./3.),400.];

%     dif     as   pc    gcv    bs   rb  rg   scale
lb = [0.0001, 600, 0.1, 0.001, 0.4, 0.1, 0.2, 5.];   % lower bounds
ub = [0.0009, 1200, 2,   0.9,  0.999, 0.5, 0.8, 500.]; % upper bounds

options = optimoptions('fmincon','Display','iter','MaxFunEvals',5000, ...
    'MaxIter',2000,'Algorithm','sqp');

%% Run the constrained optimization
% For each candidate parameter set, objective_fun runs the fisher4 PDE
% simulation and scores it by squared error against the target density.
[optimized_params, fval] = fmincon(@(params) objective_fun(params, ...
    concentration_data_a_6s, concentration_data_s_6s, concentration_data_c_6s), ...
    initial_guesses, [], [], [], [], lb, ub, [], options);

%% Extract and report the optimized parameters
dif_optimized = optimized_params(1);
as_optimized = optimized_params(2);
p_C_optimized = optimized_params(3);
gcv_optimized = optimized_params(4);
b_S_optimized = optimized_params(5);
rb_optimized = optimized_params(6);
rg_optimized = optimized_params(7);
scale_optimized = optimized_params(8);

disp('Optimized parameters:');
disp(['dif:', num2str(dif_optimized)]);
disp(['as:', num2str(as_optimized)]);
disp(['p_C:', num2str(p_C_optimized)]);
disp(['gcv:', num2str(gcv_optimized)]);
disp(['b_S:', num2str(b_S_optimized)]);
disp(['rb:', num2str(rb_optimized)]);
disp(['rg:', num2str(rg_optimized)]);
disp(['scale:', num2str(scale_optimized)]);

% Re-run the simulation with the optimized parameters and plot the result
% (figures 4-9); see fisher_plot.m for details. Note fisher_plot re-reads
% s3.png rather than reusing the s39.png data loaded above.
fisher_plot(dif_optimized,as_optimized,p_C_optimized,gcv_optimized,b_S_optimized,rb_optimized,rg_optimized,scale_optimized);
