function predicted_concentration=fisher4_restart(as,dif,pc)
% FISHER4_RESTART  Simulate the tip/stalk/collision reaction-diffusion
% model with configurable initial conditions and a mid-run domain
% "restart" (bisection) step, then plot the results at several stages.
%
% Numerical simulation of a Fisher-type equation u_t = u_xx + u^2 solved
% on a square domain. The solution can exhibit finite-time blow-up; u is
% assumed to equal zero along all boundaries.
%
% Inputs (scalars):
%   as  - saturation threshold (inactivation)
%   dif - diffusion coefficient
%   pc  - fraction of saturated tips that collide
% (bs, gcv, rb, rg and the sigmoid width b are hardcoded below rather
% than passed as inputs)
%
% Output:
%   predicted_concentration - nbox x (3*nbox) matrix formed by
%   horizontally concatenating the final active-tip, stalk, and
%   inactive/collided density grids [A | S | C] 

image = imread('s3.png');

% Classify each pixel by color into the tracked populations:
%   red    -> active tips
%   black  -> stalks
%   yellow -> inactive/collided entities
% (green is also computed but not used further)
% NOTE ON ANNOTATION COLOURS: in the manuscript figures, active tips are
% shown in green, to match the Six2-GFP signal in the micrographs. In the
% annotation files read here, active tips are instead marked in red
% for practical reasons: red marks stand out better against the green
% fluorescence of the micrographs, whereas green marks would have been
% hard to see.

black_channel = sum(image, 3) == 0; % black assumed to be RGB [0 0 0]
red_channel = image(:,:,1)-image(:,:,2)==255;
yellow_channel = image(:,:,2)-image(:,:,3)==255;
green_channel = image(:,:,2)-image(:,:,1)>50;

% Initial Variable Setup
% These flags select which initial configuration and which extra rules
% are active; only one of collision/cut/tconf/tcut should be 1 at a time.
collision = 1; % 1 = build two stalks approaching head-on from the image data (see below)
cut = 0;       % 1 = alternative "cut" initial layout (top/bottom shifted oppositely)
print=0;       % 1 = enable console debug output during the stall/restart transition
crow=0;        % 1 = enable the "crowding" rule that suppresses growth in saturated neighborhoods
tconf=0;       % 1 = build a synthetic "T"-shaped stalk from scratch (ignores the image)
tcut=0;        % 1 = build two synthetic stalk segments from scratch (ignores the image)
bisect=0;      % 1 = cut the domain along the vertical midline partway through the run and continue (the "restart" feature)
I=size(image)
f = 900;%900;%930; % number of "time units" simulated in the first run; each unit is split into 10 steps below (f*10 steps total)
rsc=1; % density rescaling factor (currently a no-op since rsc=1)
bs=0.97623; % branching-saturation fraction (gcb is centered on bs*as)
gcv=0.014; % stall-condition threshold
nbox = 20; % number of boxes to divide the lateral field into (spatial grid resolution)
trs = nbox/5; % translation distance used by the collision/cut configurations below
npix = fix(I(1)/nbox)% pixels per box; assumes the image is square, fix() keeps the integer part
A = zeros(nbox,nbox,f*10); % active-tip density over time
Ar = zeros(nbox,nbox,f*10); % rotated copy of A, used during initial-condition setup
S = zeros(nbox,nbox,f*10); % stalk density over time
Sr = zeros(nbox,nbox,f*10); % rotated copy of S, used during initial-condition setup
C = zeros(nbox,nbox,f*10); % inactive/collided density over time
Exc = zeros(nbox,nbox,f*10); % total density (A+S+C), filled in after the first time loop
stall = zeros(nbox,nbox);% stall condition per box: 1 = active, 0 = stalled
stall(:,:)=1;% initially all boxes are active (not stalled)
Cden = zeros(1,f*10); % time series of the total collided (C) density, integrated over space
Der = zeros(1,f*10); % time series of the total A.*Exc product, integrated over space
X = linspace(0,1,nbox);% row vector of nbox evenly spaced points between 0 and 1
Y = linspace(0,1,nbox);% in this model's scale, 0->1 corresponds to about 0->1 mm
spv=0.0000;% floor value used to clamp negative densities (see loop below)
A(:,:,1)=spv;

% Calculate the initial density in each box by counting matching pixels
for j=1:nbox % j is the "row" direction
    nr = npix*(j-1);
    for i=1:nbox
        nc = npix*(i-1);
        for k = 1:npix
            for m = 1:npix
            A(j,i,1) = A(j,i,1) + red_channel(nr+k,nc+m);
            S(j,i,1) = S(j,i,1) + black_channel(nr+k,nc+m);
            C(j,i,1) = C(j,i,1) + yellow_channel(nr+k,nc+m);
            end
        end
    end
end

% "collision" configuration: translate all stalks/tips to the left, then
% mirror them to build a second stalk approaching from the right, so the
% two fronts grow toward each other and collide head-on.
if (collision == 1)
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
        if (A(j,i,1)>spv)
            A(j,i-trs,1) = A(j,i,1);
            A(j,i,1) = 0.;
        end
        if (S(j,i,1)>0)
            S(j,i-trs,1) = S(j,i,1);
            S(j,i,1) = 0.;
        end
        if (C(j,i,1)>0)
            C(j,i-trs,1) = C(j,i,1);
            C(j,i,1) = 0.;
        end
    end
end
% Build the symmetric (mirrored) stalk on the right
for j=1:nbox % j is the "row" direction
    for i=1:nbox/2 % column direction
        if (A(j,i,1)>spv)
            A(j,nbox-i+1,1) = A(j,i,1);
        end
        if (S(j,i,1)>0)
            S(j,nbox-i+1,1) = S(j,i,1);
        end
        if (C(j,i,1)>0)
            C(j,nbox-i+1,1) = 0.;
        end
    end
end
end

% "cut" configuration: alternative initial layout where the top half is
% shifted left and the bottom half is shifted right (disabled by default, cut=0)
if (cut == 1)
for j=1:nbox/2 % j is the "row" direction
    for i=1:nbox % column direction
        if (A(j,i,1)>spv)
            A(j,i-trs,1) = A(j,i,1);
            A(j,i,1) = 0.;
        end
        if (S(j,i,1)>0)
            S(j,i-trs,1) = S(j,i,1);
            S(j,i,1) = 0.;
        end
        if (C(j,i,1)>0)
            C(j,i-trs,1) = C(j,i,1);
            C(j,i,1) = 0.;
        end
    end
end
% Build the shifted stalk on the bottom half, offset the other way
ncut=nbox/2+1;
trs =3;
for j=ncut:nbox % j is the "row" direction
    for i=nbox:-1:1 % column direction
        if (A(j,i,1)>spv)
            A(j,i+trs,1) = A(j,i,1);
            A(j,i,1)=0.;
        end
        if (S(j,i,1)>0)
            S(j,i+trs,1) = S(j,i,1);
            S(j,i,1)=0.;
        end
        if (C(j,i,1)>0)
            C(j,i+trs,1) = 0.;
            C(j,i,1)=0.;
        end
    end
end
end

% "T" configuration: build a synthetic T-shaped stalk with two active
% tips from scratch, ignoring the image data (disabled by default, tconf=0)
if (tconf == 1)
% Initial setup to 0
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
        A(j,i,1)=0.;
        S(j,i,1)=0.;
        C(j,i,1)=0.;
    end
end
for j=1:nbox/2 % j is the "row" direction
    S(j,nbox/2,1)=900.;
end
for i=8:12  % i is the "column" direction
    S(nbox/2+1,i,1)=900.;
end
A(nbox/2+1,7,1)=800.;%left tip
A(nbox/2+1,13,1)=800.;%right tip
end

% "tcut" configuration: two separate synthetic stalk segments, one on
% each side, built from scratch (disabled by default, tcut=0)
if (tcut == 1)
% Initial setup to 0
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
        A(j,i,1)=0.;
        S(j,i,1)=0.;
        C(j,i,1)=0.;
    end
end
% left segment
for i=8:9  % i is the "column" direction
    S(nbox/2+4,i,1)=900.;
end
A(nbox/2+4,7,1)=800.; % left tip
% right segment
for i=11:12  % i is the "column" direction
    S(nbox/2-4,i,1)=900.;
end
A(nbox/2-4,13,1)=800.; % right tip
end

% Numerical Parameters
dx = 1/nbox;
dt = 0.0015;

% Normalize pixel counts into densities: each tip corresponds to 953 pixels
A=A/(953*dx*dx);
S=S/(953*dx*dx);
C=C/(953*dx*dx);

% Rescaling of densities (no-op here since rsc=1)
A=A/rsc;
S=S/rsc;
C=C/rsc;
ma=max(max(A(:,:,1)));
A0 = A(:,:,1); % snapshot of the initial active-tip density (not referenced further below)
S0 = S(:,:,1); % snapshot of the initial stalk density (not referenced further below)
C0 = C(:,:,1); % snapshot of the initial inactive/collided density (not referenced further below)
Exc0 = A(:,:,1)+S(:,:,1)+C(:,:,1); % snapshot of the initial total density (not referenced further below)
A0den = 0.;
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
            A0den =  A0den + A(j,i,1);
    end
end
A0den

% Rotate the initial configuration about the centre of the grid.
% thetar=0 below, so this rotation is currently a no-op; the mechanism
% is left in place to allow rotating the initial conditions if desired.
     thetar =0;
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
     if  A(j,i,1) > 0.0001
     it = i - nbox/2;
     jt = j - nbox/2;
     ir = round(cos(thetar)*it-sin(thetar)*jt)+nbox/2;
     jr = round(sin(thetar)*it+cos(thetar)*jt)+nbox/2;
     Ar(jr,ir,1) = A(j,i,1);
     end
    end
end
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
     if  S(j,i,1) > 0.0001
     it = i - nbox/2;
     jt = j - nbox/2;
     ir = round(cos(thetar)*it-sin(thetar)*jt)+nbox/2;
     jr = round(sin(thetar)*it+cos(thetar)*jt)+nbox/2;
     Sr(jr,ir,1) = S(j,i,1);
     end
    end
end
A(:,:,1) = Ar(:,:,1);
S(:,:,1) = Sr(:,:,1);

%% Plot the initial conditions (figures 1-2)
figure(1)
surf(X, Y, A(:,:,1)*0.04444, 'LineStyle', 'none')
axis([0 1 0 1 0 300])
axis off
view(0, 90)
set(gca, 'Color', 'w', 'XTick', [], 'YTick', [], 'ZTick', [])
shading interp

% --- Custom green colormap: white -> light green -> dark green ---
n = 256; % number of colors
greenMap = [linspace(1, 0, n)', linspace(1, 0.6, n)', linspace(1, 0, n)'];
colormap(greenMap)

hcb = colorbar('eastoutside');
caxis([0, 18]);
hcb.FontSize = 28;
hcb.Label.String = 'Concentration (mm^{-2})';
text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

figure(2)
surf(X, Y, S(:,:,1)*0.0444444, 'LineStyle', 'none')
axis([0 1 0 1 0 300])
axis off
view(0, 90)
set(gca, 'Color', 'w', 'XTick', [], 'YTick', [], 'ZTick', [])
shading interp

% --- Custom white -> black colormap ---
n = 256; % number of color levels
grayMap = [linspace(1, 0, n)', linspace(1, 0, n)', linspace(1, 0, n)'];
colormap(grayMap)

hcb = colorbar('eastoutside');
hcb.FontSize = 28;
hcb.Label.String = 'Concentration (mm^{-2})';
text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

%% First simulation run (t = 2 to f*10)
rb=1./(9./3.); % branching rate; rb must be larger than rg
rg=1./(4./3.); % growth rate
b=as/408.; % width of the sigmoid gates (gc, gcb) used below
gam=0; % coefficient of the "ste" extra term; 0 disables it

for t = 2:f*10
    for i=2:nbox-1
        for j=2:nbox-1
          % Neighboring box values used for finite differences
          sti=stall(i+1,j)*stall(i-1,j); % if a neighboring box in i has stalled, its derivative contribution is 0
          stj=stall(i,j+1)*stall(i,j-1); % same, for the j direction
          At=A(i,j,t-1)+A(i+1,j,t-1)+A(i-1,j,t-1)+A(i,j+1,t-1)+A(i,j-1,t-1);
          cip1j=A(i+1,j,t-1)+S(i+1,j,t-1)+C(i+1,j,t-1);
          cim1j=A(i-1,j,t-1)+S(i-1,j,t-1)+C(i-1,j,t-1);
          cijp1=A(i,j+1,t-1)+S(i,j+1,t-1)+C(i,j+1,t-1);
          cijm1=A(i,j-1,t-1)+S(i,j-1,t-1)+C(i,j-1,t-1);
          cip1jp1=A(i+1,j+1,t-1)+S(i+1,j+1,t-1)+C(i+1,j+1,t-1);
          cim1jp1=A(i-1,j+1,t-1)+S(i-1,j+1,t-1)+C(i-1,j+1,t-1);
          cip1jm1=A(i+1,j-1,t-1)+S(i+1,j-1,t-1)+C(i+1,j-1,t-1);
          cim1jm1=A(i-1,j-1,t-1)+S(i-1,j-1,t-1)+C(i-1,j-1,t-1);
          cij=A(i,j,t-1)+S(i,j,t-1)+C(i,j,t-1); % total local concentration (all 3 species) at box (i,j)
          Ct=cij+cip1j+cim1j+cijp1+cijm1;
          aip1j=A(i+1,j,t-1);
          aim1j=A(i-1,j,t-1);
          aijp1=A(i,j+1,t-1);
          aijm1=A(i,j-1,t-1);
          aip1jp1=A(i+1,j+1,t-1);
          aim1jp1=A(i-1,j+1,t-1);
          aip1jm1=A(i+1,j-1,t-1);
          aim1jm1=A(i-1,j-1,t-1);
          aij=A(i,j,t-1);
          z=(cij-as)/b;
          zip1j=(cip1j-as)/b;
          zim1j=(cim1j-as)/b;
          zijp1=(cijp1-as)/b;
          zijm1=(cijm1-as)/b;
          gc=exp(-z)/(exp(z)+exp(-z));
          % gcip1j..gcijm1 below are computed but not used further in this loop
          gcip1j=exp(-zip1j)/(exp(zip1j)+exp(-zip1j));
          gcim1j=exp(-zim1j)/(exp(zim1j)+exp(-zim1j));
          gcijp1=exp(-zijp1)/(exp(zijp1)+exp(-zijp1));
          gcijm1=exp(-zijm1)/(exp(zijm1)+exp(-zijm1));

          % Sigmoid "gate" functions controlling inactivation (gc) and
          % branching saturation (gcb), both centered on the saturation
          % threshold "as" (gcb centered on bs*as instead) with width "b"
          zb=(cij-bs*as)/b;
          gcb=exp(-zb)/(exp(zb)+exp(-zb));
          dgc=-(2/b)*(1/(exp(z)+exp(-z))^2); % derivative of gc

          % Finite-difference derivatives of the total concentration c and
          % of the active-tip density a (Laplacians, gradients, mixed
          % second derivatives), with stall gating (sti, stj) applied to c
          d2c=(cip1j+cim1j-2*cij)*sti/dx^2+(cijp1+cijm1-2*cij)*stj/dx^2;
          d2a=(aip1j+aim1j+aijp1+aijm1-4*aij)/dx^2;
          dci=(cip1j-cim1j)*sti/(2*dx); %derivative respect to i
          dcj=(cijp1-cijm1)*stj/(2*dx); %derivative respect to j
          ddcii=(cip1j+cim1j-2*cij)/dx^2; %second derivative to i
          ddcjj=(cijp1+cijm1-2*cij)/dx^2; %second derivative to j
          ddcij=(cip1jp1+cim1jm1-cim1jp1-cip1jm1)/(4*dx^2); %mixed derivative to i and j
          dai=(aip1j-aim1j)/(2*dx); %derivative respect to i
          daj=(aijp1-aijm1)/(2*dx); %derivative respect to j
          ddaii=(aip1j+aim1j-2*aij)/dx^2; %second derivative to i
          ddajj=(aijp1+aijm1-2*aij)/dx^2; %second derivative to j
          ddaij=(aip1jp1+aim1jm1-aim1jp1-aip1jm1)/(4*dx^2); %mixed derivative to i and j
          mdc=sqrt(dci^2+dcj^2);%modulus of gradient of c
          m2dc = mdc*mdc;% square modulus of gradient of c
          mda=sqrt(dai^2+daj^2);%modulus of gradient of a
          m2da = mda*mda;% square modulus of gradient of a
          dadc = dci*dai + dcj*daj; %scalar product da . dc
          ste=gam*(mda^2*gc+aij*d2a*gc+aij*dadc*dgc); % extra term, disabled since gam=0

          % "Crowding" rule: if enabled (crow=1), active growth is
          % suppressed when too many of the 8 neighboring boxes already
          % have a high active-tip density, so tips cannot flow freely
          % like a liquid into already-crowded regions.
          if crow==1
          count=0;
          thrs=40.;
          if aip1j>thrs
              count=count+1;
          end
          if aim1j>thrs
              count=count+1;
          end
          if aijp1>thrs
              count=count+1;
          end
          if aijm1>thrs
              count=count+1;
          end
          if aip1jp1>thrs
              count=count+1;
          end
          if aim1jp1>thrs
              count=count+1;
          end
          if aip1jm1>thrs
              count=count+1;
          end
          if aim1jm1>thrs
              count=count+1;
          end

          if count>6
             gc=0;
          end
          end

          % Update the three densities for boxes with enough surrounding
          % concentration that have not stalled; otherwise freeze them
          if Ct > 0.1 & stall(i,j)>0.1 % check if there is some concentration around
            if cij<90. % clamp cij to avoid division blow-up when concentration is very low
                cij=90.;
            end
%            A(i,j,t) = aij + dt*dif*(d2c*gc*aij/cij+m2dc*dgc*aij/cij+dadc*gc/cij-m2dc*gc*aij/cij^2)+dt*rb*aij*gcb-dt*pc*rb*aij*(1-gcb)+dt*ste;
            A(i,j,t) = aij + dt*dif*(d2c*gc*aij/cij+m2dc*dgc*aij/cij+dadc*gc/cij-m2dc*gc*aij/cij^2)+dt*rb*gcb*aij*gc-dt*pc*rb*aij*(1-gc)+dt*ste;
            S(i,j,t) = S(i,j,t-1) + dt*rg*aij*gc;
%            C(i,j,t)=C(i,j,t-1)+dt*rb*pc*aij*(1-gcb);
            C(i,j,t)=C(i,j,t-1)+dt*rb*pc*aij*(1-gc);
          else
            A(i,j,t) = A(i,j,t-1);
            S(i,j,t) = S(i,j,t-1);
            C(i,j,t) = C(i,j,t-1);
          end
          if A(i,j,t) < -0.0001
             A(i,j,t)=spv; % clamp negative density to the floor value
          end

          % Check the stall condition: a box stalls once its gate value
          % drops below gcv
        if gc<gcv
          stall(i,j)=0;
        end
          % Check the restart condition: a box can restart once its gate
          % value rises back above 0.1 
        if gc>0.1
          stall(i,j)=1;
        end

          % Debug: set print=1 to echo state when a box is about to stall
         if gc<0.1 & i~=13 & print==1
             fprintf('t=%d\n',t);
             fprintf('j=%d\n',j);
             fprintf('i=%d\n',i);
             fprintf('gcb=%d\n',gcb);
         end
        end
    end
end

% Snapshot the three densities at the current time t
concentration_s_pred = S(:,:,t);
concentration_a_pred = A(:,:,t);
concentration_c_pred = C(:,:,t);
pred_concentration_data = {concentration_a_pred, concentration_s_pred, concentration_c_pred};
predicted_concentration = cell2mat(pred_concentration_data);
Exc = A+S+C; % total density over time, used below for the KDE/integration steps

%% Plot the state after the first run (figures 3-5)
        figure(3)
        surf(X,Y,(A(:,:,t))*0.044444,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % --- Custom green colormap: white -> light green -> dark green ---
        n = 256; % number of colors
        greenMap = [linspace(1, 0, n)', linspace(1, 0.6, n)', linspace(1, 0, n)'];
        colormap(greenMap)
        hcb = colorbar('eastoutside');
        caxis([0, 18]);
        hcb.FontSize = 28;
        hcb.Label.String = 'Concentration (mm^{-2})';
        text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

        figure(4)
        surf(X,Y,(S(:,:,t))*0.0444444,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % --- Custom white -> black colormap ---
        n = 256; % number of color levels
        grayMap = [linspace(1, 0, n)', linspace(1, 0, n)', linspace(1, 0, n)'];
        colormap(grayMap)
hcb = colorbar('eastoutside');
hcb.FontSize = 28;
hcb.Label.String = 'Concentration (mm^{-2})';
text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

        figure(5)
        surf(X,Y,C(:,:,t)*0.0444444,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % Define white-to-dark-orange-yellow colormap
        n = 256; % number of color levels
        white_to_orangeyellow = [ ...
        linspace(1, 0.8, n)', ...   % Red: stays strong
        linspace(1, 0.4, n)', ...   % Green: reduced for warmer tone
        linspace(1, 0.0, n)'  ...   % Blue: fades to 0 (adds orange depth)
];

        colormap(white_to_orangeyellow)
        hcb = colorbar('eastoutside');
        hcb.FontSize = 20;
        hcb.Label.String = 'Concentration (mm^{-2})';
        text(1.18, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 20);

% Integrate the density of collisions: D is the pointwise product of
% active-tip density and total density; Cden/Der are its spatial sums
% over time, computed for every time step k
D = A.*Exc;
for k=1:f*10
  for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
             Cden(k) =  Cden(k) + C(j,i,k);
             Der(k) = Der(k) + D(j,i,k);
    end
  end
end

% Time-series probes at two fixed boxes: (10,10) and (12,10). In
% "collision" mode, Sc tracks only the active-tip density (the collision
% point) and Ct tracks the full local concentration; otherwise both Sc
% and Sd track the full local concentration.
if (collision == 1)
   for k=1:f*10
     Sc(k) = A(10,10,k);
     Ct(k) = A(10,10,k)+S(10,10,k)+C(10,10,k);
   end
   for k=1:f*10
     Sd(k) = A(12,10,k);
   end
figure(6)
% Convert nsteps (f*10) to hours
hours= dt*f*10*2;
tim = linspace(0, hours, length(Sc));
%tim = linspace(0, 26, length(Sc));
% Plot with the rescaled (hours) x-axis
plot(tim, Sc * 0.04444, 'g-', 'LineWidth', 4, 'MarkerFaceColor', 'g');
hold on
plot(tim, Ct * 0.04444, 'b-', 'LineWidth', 4, 'MarkerFaceColor', 'b');

% Axis formatting
xlim([0 hours]);
ylim([0 60]);
xlabel('t (hours)', 'FontSize', 28);
ylabel('Concentration (mm^{-2})', 'FontSize', 28);
ax = gca;
ax.XAxis.FontSize = 32;
ax.YAxis.FontSize = 32;

legend({'active(a)', 'total(c)'}, ...
       'FontSize', 32, ...
       'Location', 'northwest');
legend boxoff   % removes the legend frame

figure(7)
plot(Sd*0.04444,'ro');
ylim([0 22]);
else
   for k=1:f*10
     Sc(k) = A(10,10,k)+S(10,10,k)+C(10,10,k);
   end
   for k=1:f*10
     Sd(k) = A(12,10,k)+S(12,10,k)+C(12,10,k);
   end
figure(6)
plot(Sc);
figure(7)
plot(Sd);
end

% Integrate the density of active tips over the whole domain at the final time t
Aden = 0.;
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
            Aden =  Aden + A(j,i,t);
    end
end
Aden

% Integrate the density of inactive/collided entities over the whole domain
Ctden = 0.;
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
            Ctden = Ctden + C(j,i,t);
    end
end
Ctden

% Flag boxes with a high active-tip density that are not stalled
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
        stalla(j,i)=0;
        if (A(j,i,t)>80 && stall(j,i)==0)
           stalla(j,i)=1;
        end
    end
end

%% Plot the "high density & not stalled" flag map (figure 8)
         figure(8)
        surf(X,Y,stalla(:,:),'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % ----- Custom colormap: white -> violet -----
n = 256;  % number of colors
white  = [1 1 1];
violet = [0.55 0.0 0.8];   % a strong violet tone (adjustable)
customMap = [linspace(white(1), violet(1), n)', ...
             linspace(white(2), violet(2), n)', ...
             linspace(white(3), violet(3), n)'];
colormap(customMap);
colorbar;

%% Kernel density estimate (KDE) of the active-tip density A
% Note: X is temporarily reused below as the KDE grid-point coordinates
% (overwriting the [0,1] plotting axis defined earlier); it is reset to
% the plotting axis again further down before being used in surf() calls.
Ak=(A(:,:,t))/10.;
% Create X and Y coordinate grids corresponding to A
[x, y] = meshgrid(1: nbox, 1: nbox);

% Flatten the matrices into vectors
X = [x(:), y(:)];         % grid coordinates as points
V = Ak(:);                 % values at each point

% Define resolution for the KDE output
res = 200;  % change for higher/lower detail

% Define the evaluation grid for the KDE
xi = linspace(1, nbox, res);
yi = linspace(1, nbox, res);
[xxi, yyi] = meshgrid(xi, yi);
gridPoints = [xxi(:), yyi(:)];

% Compute the KDE using the matrix values as weights
pdfxy = ksdensity(X, gridPoints, 'Weights', V, 'Bandwidth', 1);

% Reshape the KDE output (a 1D vector of densities) into a 2D array
% matching the evaluation mesh grid
pdfxy = reshape(pdfxy, res, res);

% Plot the KDE surface
figure;
mesh(xxi, yyi, pdfxy);
xlabel('X'); ylabel('Y'); zlabel('Density Estimate');
title('Kernel Density Estimate of Matrix A');
set(gca, 'XLim', [min(xi) max(xi)]);
set(gca, 'YLim', [min(yi) max(yi)]);
colorbar

% Plot two KDE contour levels (70% and 60% of the peak density)
figure;
contourLevel = max(pdfxy(:)) * 0.7;  % adjust threshold as needed - ref value 0.75
contour(xi, yi, pdfxy, [contourLevel contourLevel], 'LineWidth', 2, 'LineColor', 'b');
hold on
contourLevel = max(pdfxy(:)) * 0.6;  % adjust threshold as needed - ref value 0.75
contour(xi, yi, pdfxy, [contourLevel contourLevel], 'LineWidth', 2, 'LineColor', 'b','LineStyle', '--');
xlabel('X'); ylabel('Y');

%% Restart feature: cut the domain along the vertical midline and
% continue the simulation from this partial state
X = linspace(0,1,nbox);% row vector of nbox evenly spaced points between 0 and 1 (reset after the KDE step above)
Y = linspace(0,1,nbox);% in this model's scale, 0->1 corresponds to about 0->1 mm
if (bisect == 1)
dcut=4.; % divisor applied to the column just right of the cut
dcut1=3.; % divisor applied to the centre column
dcut2=2.; % divisor applied to the column just left of the cut
for j=1:nbox   % j is the "row" direction
    for i=nbox/2+2:nbox % zero out everything right of the cut
            A(j,i,t) = 0.;
            S(j,i,t) = 0.;
            C(j,i,t) = 0.;
            stall(j,i) = 1;
    end
            A(j,nbox/2+1,t) = A(j,nbox/2+1,t)/dcut;
            S(j,nbox/2+1,t) = S(j,nbox/2+1,t)/dcut;
            C(j,nbox/2+1,t) = C(j,nbox/2+1,t)/dcut;
%            stall(j,nbox/2+1) = 1;
            A(j,nbox/2,t) = A(j,nbox/2,t)/dcut1;
            S(j,nbox/2,t) = S(j,nbox/2,t)/dcut1;
            C(j,nbox/2,t) = C(j,nbox/2,t)/dcut1;
%            stall(j,nbox/2) = 1;
            A(j,nbox/2-1,t) = A(j,nbox/2-1,t)/dcut2;
            S(j,nbox/2-1,t) = S(j,nbox/2-1,t)/dcut2;
            C(j,nbox/2-1,t) = C(j,nbox/2-1,t)/dcut2;
%            stall(j,nbox/2-1) = 1;
end
        figure(11)
        surf(X,Y,(A(:,:,t))*0.044444,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % --- Custom green colormap: white -> light green -> dark green ---
        n = 256; % number of colors
        greenMap = [linspace(1, 0, n)', linspace(1, 0.6, n)', linspace(1, 0, n)'];
        colormap(greenMap)
        hcb = colorbar('eastoutside');
        caxis([0, 18]);
        hcb.FontSize = 28;
        hcb.Label.String = 'Concentration (mm^{-2})';
        text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

        figure(12)
        surf(X,Y,(S(:,:,t))*0.0444444,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % --- Custom white -> black colormap ---
        n = 256; % number of color levels
        grayMap = [linspace(1, 0, n)', linspace(1, 0, n)', linspace(1, 0, n)'];
        colormap(grayMap)
hcb = colorbar('eastoutside');
hcb.FontSize = 28;
hcb.Label.String = 'Concentration (mm^{-2})';
text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

% Flag boxes with a high active-tip density that are not stalled (after the cut)
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
        stalla(j,i)=0;
        if (A(j,i,t)>80 && stall(j,i)==0)
           stalla(j,i)=1;
        end
    end
end
         figure(13)
        surf(X,Y,stalla(:,:),'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])

        shading(gca,'interp')

        % ----- Custom colormap: white -> violet -----
n = 256;  % number of colors

white  = [1 1 1];
violet = [0.55 0.0 0.8];   % a strong violet tone (adjustable)

customMap = [linspace(white(1), violet(1), n)', ...
             linspace(white(2), violet(2), n)', ...
             linspace(white(3), violet(3), n)'];
colormap(customMap);
colorbar;

%% Continue the time evolution after the cut (t = f*10+1 to f*15)
 for t = f*10+1:f*15
    for i=2:nbox-1
        for j=2:nbox-1
          sti=stall(i+1,j)*stall(i-1,j); % if a neighboring box in i has stalled, its derivative contribution is 0
          stj=stall(i,j+1)*stall(i,j-1); % same, for the j direction
          At=A(i,j,t-1)+A(i+1,j,t-1)+A(i-1,j,t-1)+A(i,j+1,t-1)+A(i,j-1,t-1);
          cip1j=A(i+1,j,t-1)+S(i+1,j,t-1)+C(i+1,j,t-1);
          cim1j=A(i-1,j,t-1)+S(i-1,j,t-1)+C(i-1,j,t-1);
          cijp1=A(i,j+1,t-1)+S(i,j+1,t-1)+C(i,j+1,t-1);
          cijm1=A(i,j-1,t-1)+S(i,j-1,t-1)+C(i,j-1,t-1);
          cip1jp1=A(i+1,j+1,t-1)+S(i+1,j+1,t-1)+C(i+1,j+1,t-1);
          cim1jp1=A(i-1,j+1,t-1)+S(i-1,j+1,t-1)+C(i-1,j+1,t-1);
          cip1jm1=A(i+1,j-1,t-1)+S(i+1,j-1,t-1)+C(i+1,j-1,t-1);
          cim1jm1=A(i-1,j-1,t-1)+S(i-1,j-1,t-1)+C(i-1,j-1,t-1);
          cij=A(i,j,t-1)+S(i,j,t-1)+C(i,j,t-1); % total local concentration (all 3 species) at box (i,j)
          Ct=cij+cip1j+cim1j+cijp1+cijm1;
          aip1j=A(i+1,j,t-1);
          aim1j=A(i-1,j,t-1);
          aijp1=A(i,j+1,t-1);
          aijm1=A(i,j-1,t-1);
          aip1jp1=A(i+1,j+1,t-1);
          aim1jp1=A(i-1,j+1,t-1);
          aip1jm1=A(i+1,j-1,t-1);
          aim1jm1=A(i-1,j-1,t-1);
          aij=A(i,j,t-1);
          z=(cij-as)/b;
          gc=exp(-z)/(exp(z)+exp(-z));

          % Sigmoid "gate" functions controlling inactivation (gc) and
          % branching saturation (gcb), both centered on the saturation
          % threshold "as" (gcb centered on bs*as instead) with width "b"
          zb=(cij-bs*as)/b;
          gcb=exp(-zb)/(exp(zb)+exp(-zb));
          dgc=-(2/b)*(1/(exp(z)+exp(-z))^2); % derivative of gc

          % Finite-difference derivatives of the total concentration c and
          % of the active-tip density a (Laplacians, gradients, mixed
          % second derivatives), with stall gating (sti, stj) applied to c
          d2c=(cip1j+cim1j-2*cij)*sti/dx^2+(cijp1+cijm1-2*cij)*stj/dx^2;
          d2a=(aip1j+aim1j+aijp1+aijm1-4*aij)/dx^2;
          dci=(cip1j-cim1j)*sti/(2*dx); %derivative respect to i
          dcj=(cijp1-cijm1)*stj/(2*dx); %derivative respect to j
          ddcii=(cip1j+cim1j-2*cij)/dx^2; %second derivative to i
          ddcjj=(cijp1+cijm1-2*cij)/dx^2; %second derivative to j
          ddcij=(cip1jp1+cim1jm1-cim1jp1-cip1jm1)/(4*dx^2); %mixed derivative to i and j
          dai=(aip1j-aim1j)/(2*dx); %derivative respect to i
          daj=(aijp1-aijm1)/(2*dx); %derivative respect to j
          ddaii=(aip1j+aim1j-2*aij)/dx^2; %second derivative to i
          ddajj=(aijp1+aijm1-2*aij)/dx^2; %second derivative to j
          ddaij=(aip1jp1+aim1jm1-aim1jp1-aip1jm1)/(4*dx^2); %mixed derivative to i and j
          mdc=sqrt(dci^2+dcj^2);%modulus of gradient of c
          m2dc = mdc*mdc;% square modulus of gradient of c
          mda=sqrt(dai^2+daj^2);%modulus of gradient of a
          m2da = mda*mda;% square modulus of gradient of a
          dadc = dci*dai + dcj*daj; %scalar product da . dc
          ste=gam*(mda^2*gc+aij*d2a*gc+aij*dadc*dgc); % extra term, disabled since gam=0

          % Update the three densities for boxes with enough surrounding
          % concentration that have not stalled; otherwise freeze them
          if Ct > 0.1 & stall(i,j)>0.1 % check if there is some concentration around
            if cij<90. % clamp cij to avoid division blow-up when concentration is very low
                cij=90.;
            end
%            A(i,j,t) = aij + dt*dif*(d2c*gc*aij/cij+m2dc*dgc*aij/cij+dadc*gc/cij-m2dc*gc*aij/cij^2)+dt*rb*aij*gcb-dt*pc*rb*aij*(1-gcb)+dt*ste;
            A(i,j,t) = aij + dt*dif*(d2c*gc*aij/cij+m2dc*dgc*aij/cij+dadc*gc/cij-m2dc*gc*aij/cij^2)+dt*rb*gcb*aij*gc-dt*pc*rb*aij*(1-gc)+dt*ste;
            S(i,j,t) = S(i,j,t-1) + dt*rg*aij*gc;
%            C(i,j,t)=C(i,j,t-1)+dt*rb*pc*aij*(1-gcb);
            C(i,j,t)=C(i,j,t-1)+dt*rb*pc*aij*(1-gc);
          else
            A(i,j,t) = A(i,j,t-1);
            S(i,j,t) = S(i,j,t-1);
            C(i,j,t) = C(i,j,t-1);
          end
          if A(i,j,t) < -0.0001
             A(i,j,t)=spv; % clamp negative density to the floor value
          end

          % Check the stall condition: a box stalls once its gate value
          % drops below gcv
        if gc<gcv
          stall(i,j)=0;

        end
          % Check the restart condition: a box can restart once its gate
          % value rises back above 0.1
        if gc>0.1
          stall(i,j)=1;
        end

          % Debug: set print=1 to echo state when a box is about to stall
         if gc<0.1 & i~=13 & print==1
             fprintf('t=%d\n',t);
             fprintf('j=%d\n',j);
             fprintf('i=%d\n',i);
             fprintf('gcb=%d\n',gcb);
         end

        end
    end

 end

%% Plot the final state after the continued evolution (figures 14-16)
         figure(14)
        surf(X,Y,(A(:,:,t))*0.044444,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % --- Custom green colormap: white -> light green -> dark green ---
        n = 256; % number of colors
        greenMap = [linspace(1, 0, n)', linspace(1, 0.6, n)', linspace(1, 0, n)'];
        colormap(greenMap)
        hcb = colorbar('eastoutside');
        caxis([0, 18]);
        hcb.FontSize = 28;
        hcb.Label.String = 'Concentration (mm^{-2})';
        text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

        figure(15)
        surf(X,Y,(S(:,:,t))*0.0444444,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        shading(gca,'interp')
        % --- Custom white -> black colormap ---
        n = 256; % number of color levels
        grayMap = [linspace(1, 0, n)', linspace(1, 0, n)', linspace(1, 0, n)'];
        colormap(grayMap)
hcb = colorbar('eastoutside');
hcb.FontSize = 28;
hcb.Label.String = 'Concentration (mm^{-2})';
text(1.15, 1.05, '×10', 'Units', 'normalized', ...
     'Rotation', 0, 'HorizontalAlignment', 'center', 'FontSize', 28);

% Flag boxes with a high active-tip density that are not stalled (final state)
for j=1:nbox % j is the "row" direction
    for i=1:nbox % column direction
        stalla(j,i)=0;
        if (A(j,i,t)>80 && stall(j,i)==0)
           stalla(j,i)=1;
        end
    end
end
         figure(16)
        surf(X,Y,stalla(:,:),'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        view(0, 90)
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])

        shading(gca,'interp')

        % ----- Custom colormap: white -> violet -----
n = 256;  % number of colors

white  = [1 1 1];
violet = [0.55 0.0 0.8];   % a strong violet tone (adjustable)

customMap = [linspace(white(1), violet(1), n)', ...
             linspace(white(2), violet(2), n)', ...
             linspace(white(3), violet(3), n)'];
colormap(customMap);
colorbar;
end % end if condition on bisect
end
