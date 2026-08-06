function fisher_plot(dif,as,pc,gcv,bs,rb,rg,scale)
% FISHER_PLOT  Re-run the fisher4 reaction-diffusion simulation with a
% given (typically already-optimized) parameter set and plot the initial
% and final densities. This function does not return a value; see
% fisher4.m for the version used during optimization, which returns the
% final active-tip density instead of plotting.
%
% Numerical simulation of a Fisher-type equation u_t = u_xx + u^2 solved
% on a square domain. The solution can exhibit finite-time blow-up; u is
% assumed to equal zero along all boundaries.
%
% Inputs (all scalars):
%   dif   - diffusion coefficient
%   as    - saturation threshold (inactivation)
%   pc    - fraction of saturated tips that collide
%   gcv   - stall-condition threshold
%   bs    - branching-saturation fraction
%   rb    - branching rate
%   rg    - growth rate
%   scale - scale factor (controls the width "b" of the sigmoid gates below)

image = imread('s3.png');

% Classify each pixel by color into the tracked populations:
%   red    -> active tips
%   black  -> stalks
%   yellow -> inactive/collided entities
% (green is also computed but not used further)
black_channel = sum(image, 3) == 0; % black assumed to be RGB [0 0 0]
red_channel = image(:,:,1)-image(:,:,2)==255;
yellow_channel = image(:,:,2)-image(:,:,3)==255;
green_channel = image(:,:,2)-image(:,:,1)>50;

% Initial Variable Setup
I=size(image)
f=1300; % number of "time units" simulated; each unit is split into 10 steps below (f*10 steps total)
nbox = 20; % number of boxes to divide the lateral field into (spatial grid resolution)
npix = fix(I(1)/nbox)% pixels per box; assumes the image is square, fix() keeps the integer part
A = zeros(nbox,nbox,f*10); % active-tip density over time
S = zeros(nbox,nbox,f*10); % stalk density over time
C = zeros(nbox,nbox,f*10); % inactive/collided density over time
Exc = zeros(nbox,nbox,f*10);
stall = zeros(nbox,nbox);% stall condition per box: 1 = active, 0 = stalled
stall(:,:)=1;% initially all boxes are active (not stalled)
X = linspace(0,1,nbox);% row vector of nbox evenly spaced points between 0 and 1
Y = linspace(0,1,nbox);% in this model's scale, 0->1 corresponds to about 0->1 mm
spv=0.0000;% floor value used to clamp negative densities (see loop below)
A(:,:,1)=spv;

% Calculate the initial density in each box by counting matching pixels
for j=1:nbox
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

% Numerical Parameters
dx = X(2) - X(1);
dt = 0.001;

% Normalize pixel counts into densities: each tip corresponds to 953 pixels
A=A/(953*dx*dx);
S=S/(953*dx*dx);
C=C/(953*dx*dx);

%% Plot the initial densities (figures 4-6)
    figure(4)
    surf(X,Y,A(:,:,1)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

    figure(5)
    surf(X,Y,S(:,:,1)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

    figure(6)
    surf(X,Y,C(:,:,1)/10,'LineStyle','none')
    axis([0 1 0 1 0 300],'off')
    set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
    lightangle(-90,45)
    shading(gca,'interp')
    colorbar

ma=max(max(A(:,:,1)))

b=as/scale; % width of the sigmoid gates (gc, gcb) used below
gam=0; % coefficient of the "ste" extra term; 0 disables it

% Time-stepping loop: simulate the reaction-diffusion system from t=2 to f*10
for t = 2:f*10
    for i=2:nbox-1
        for j=2:nbox-1
            %-----New variables
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

          % Sigmoid "gate" functions controlling inactivation (gc) and
          % branching saturation (gcb), both centered on the saturation
          % threshold "as" (gcb centered on bs*as instead) with width "b"
          z=(cij-as)/b;
          gc=exp(-z)/(exp(z)+exp(-z));
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
          if Ct > 0.1 & stall(i,j)>0.1
            cij=cij+20.;
            A(i,j,t) = aij + dt*dif*(d2c*gc*aij/cij+m2dc*dgc*aij/cij+dadc*gc/cij-m2dc*gc*aij/cij^2)+dt*rb*aij*gcb-dt*pc*rb*aij*(1-gcb)+dt*ste;
            S(i,j,t) = S(i,j,t-1) + dt*rg*aij*gc;
            C(i,j,t)=C(i,j,t-1)+dt*rb*pc*aij*(1-gcb);
          else
            A(i,j,t) = A(i,j,t-1);
            S(i,j,t) = S(i,j,t-1);
            C(i,j,t) = C(i,j,t-1);
          end
          if A(i,j,t) < -0.0001
             A(i,j,t)=spv; % clamp negative density to the floor value
          end

          % Update the stall condition: a box stalls once its gate value
          % drops below gcv, and can restart once gc rises back above 0.9
        if gc<gcv
          stall(i,j)=0;
        end
        if gc>0.9
          stall(i,j)=1;
        end
        end
    end
end

concentration_s_pred = S(:,:,t);
concentration_a_pred = A(:,:,t);
concentration_c_pred = C(:,:,t);

%% Plot the final densities (figures 7-9)
        figure(7)
        surf(X,Y,A(:,:,t)/10,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        lightangle(-90,45)
        shading(gca,'interp')
        colorbar

        figure(8)
        surf(X,Y,S(:,:,t)/10,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        lightangle(-90,45)
        shading(gca,'interp')
        colorbar

        figure(9)
        surf(X,Y,C(:,:,t)/10,'LineStyle','none')
        axis([0 1 0 1 0 300],'off')
        set(gca,'Color','w','XTick',[],'YTick',[],'ZTick',[])
        lightangle(-90,45)
        shading(gca,'interp')
        colorbar
end
