function [At, St, Ct] = fisher_evol(f, as, dif, pc, gcv, bs, rb, rg, scale, A, S, C)
    % FISHER_EVOL  Run the tip/stalk/collision reaction-diffusion model
    % starting from precomputed initial density grids, and return the
    % density grids at the final time step.
    %
    % This function takes the initial active-tip (A), stalk (S), and
    % inactive/collided (C) density grids as inputs, so it can be reused
    % to re-simulate from the same starting point with different
    % parameters (as done by the sensitivity analysis scripts).
    %
    % Inputs:
    %   f            - number of "time units" simulated; each unit is split
    %                  into 10 steps below (f*10 steps total)
    %   as           - saturation threshold (inactivation)
    %   dif          - diffusion coefficient
    %   pc           - inactive-tip production rate
    %   gcv          - stall-condition threshold
    %   bs           - branching-saturation gate centre. Note: unlike
    %                  fisher4.m (which computes zb=(cij-bs*as)/b from a
    %                  fraction bs), this function uses zb=(cij-bs)/b
    %                  directly, so callers must pass bs already
    %                  pre-multiplied by as (e.g. bs_ref = 0.97623*as_ref
    %                  in sensitivity1.m).
    %   rb           - branching rate
    %   rg           - growth rate
    %   scale        - scale factor (controls the width "b" of the sigmoid
    %                  gates below)
    %   A, S, C      - initial active-tip, stalk, and inactive/collided
    %                  density grids (nbox x nbox)
    %
    % Outputs:
    %   At, St, Ct   - active-tip, stalk, and inactive/collided density
    %                  grids at the final time step. Note: inside the time
    %                  loop below, "At" and "Ct" are also used as local
    %                  scalar temporaries (a sum of neighboring densities,
    %                  and the total local concentration, respectively);
    %                  those temporary uses are unrelated to these output
    %                  values, which are only assigned at the very end.

    % Numerical parameters
    [nbox, ~, ~] = size(A);
    X = linspace(0, 1, nbox); % row vector of nbox evenly spaced points between 0 and 1
    Y = linspace(0, 1, nbox); % in this model's scale, 0->1 corresponds to about 0->1 mm
    spv=0.0000; % floor value used to clamp negative densities (see loop below)
    dx = 1/nbox;
    dt = 0.0015;
    Af = zeros(nbox,nbox,f*10); % active-tip density over time
    Sf = zeros(nbox,nbox,f*10); % stalk density over time
    Cf = zeros(nbox,nbox,f*10); % inactive/collided density over time
    Af(:,:,1) = A;
    Sf(:,:,1) = S;
    Cf(:,:,1) = C;
    stall = zeros(nbox,nbox);% stall condition per box: 1 = active, 0 = stalled
    stall(:,:)=1;% initially all boxes are active (not stalled)
    b=as/scale; % width of the sigmoid gates (gc, gcb) used below
    gam=0; % coefficient of the "ste" extra term; 0 disables it
    crow=0; % 1 = enable the "crowding" rule that suppresses growth in saturated neighborhoods

    % Time-stepping loop: simulate the reaction-diffusion system from t=2 to f*10
for t = 2:f*10
    for i=2:nbox-1
        for j=2:nbox-1
          % Neighboring box values used for finite differences
          sti=stall(i+1,j)*stall(i-1,j); % if a neighboring box in i has stalled, its derivative contribution is 0
          stj=stall(i,j+1)*stall(i,j-1); % same, for the j direction
          At=Af(i,j,t-1)+Af(i+1,j,t-1)+Af(i-1,j,t-1)+Af(i,j+1,t-1)+Af(i,j-1,t-1);
          cip1j=Af(i+1,j,t-1)+Sf(i+1,j,t-1)+Cf(i+1,j,t-1);
          cim1j=Af(i-1,j,t-1)+Sf(i-1,j,t-1)+Cf(i-1,j,t-1);
          cijp1=Af(i,j+1,t-1)+Sf(i,j+1,t-1)+Cf(i,j+1,t-1);
          cijm1=Af(i,j-1,t-1)+Sf(i,j-1,t-1)+Cf(i,j-1,t-1);
          cip1jp1=Af(i+1,j+1,t-1)+Sf(i+1,j+1,t-1)+Cf(i+1,j+1,t-1);
          cim1jp1=Af(i-1,j+1,t-1)+Sf(i-1,j+1,t-1)+Cf(i-1,j+1,t-1);
          cip1jm1=Af(i+1,j-1,t-1)+Sf(i+1,j-1,t-1)+Cf(i+1,j-1,t-1);
          cim1jm1=Af(i-1,j-1,t-1)+Sf(i-1,j-1,t-1)+Cf(i-1,j-1,t-1);
          cij=Af(i,j,t-1)+Sf(i,j,t-1)+Cf(i,j,t-1); % total local concentration (all 3 species) at box (i,j)
          Ct=cij+cip1j+cim1j+cijp1+cijm1;
          aip1j=Af(i+1,j,t-1);
          aim1j=Af(i-1,j,t-1);
          aijp1=Af(i,j+1,t-1);
          aijm1=Af(i,j-1,t-1);
          aip1jp1=Af(i+1,j+1,t-1);
          aim1jp1=Af(i-1,j+1,t-1);
          aip1jm1=Af(i+1,j-1,t-1);
          aim1jm1=Af(i-1,j-1,t-1);
          aij=Af(i,j,t-1);
          z=(cij-as)/b;
          gc=exp(-z)/(exp(z)+exp(-z));
          zb=(cij-bs)/b;
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
          mdc=sqrt(dci^2+dcj^2);%modulus of gradient
          m2dc = mdc*mdc;% square modulus of gradient
          mda=sqrt(dai^2+daj^2);%modulus of gradient
          m2da = mda*mda;% square modulus of gradient
          dadc = dci*dai + dcj*daj; %scalar product da * dc

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
 %           Af(i,j,t) = aij + dt*dif*(d2c*gc*aij/cij+m2dc*dgc*aij/cij+dadc*gc/cij-m2dc*gc*aij/cij^2)+dt*rb*aij*gcb-dt*pc*aij*(1-gcb)+dt*ste;
            Af(i,j,t) = aij + dt*dif*(d2c*gc*aij/cij+m2dc*dgc*aij/cij+dadc*gc/cij-m2dc*gc*aij/cij^2)+dt*rb*aij*gcb*gc-dt*pc*aij*(1-gc)+dt*ste;
            Sf(i,j,t) = Sf(i,j,t-1) + dt*rg*aij*gc;
%            Cf(i,j,t) = Cf(i,j,t-1)+dt*pc*aij*(1-gcb);
            Cf(i,j,t) = Cf(i,j,t-1)+dt*pc*aij*(1-gc);
          else
            Af(i,j,t) = Af(i,j,t-1);
            Sf(i,j,t) = Sf(i,j,t-1);
            Cf(i,j,t) = Cf(i,j,t-1);
          end
          if Af(i,j,t) < -0.0001
             Af(i,j,t)=spv; % clamp negative density to the floor value
          end

          % Check the stall condition: a box stalls once its gate value
          % drops below gcv
        if gc<gcv
          stall(i,j)=0;
        end
          % Check the restart condition: a box can restart once its gate
          % value rises back above 0.9
        if gc>0.1
          stall(i,j)=1;
        end

        end
    end
end
    At = Af(:, :, t);
    St = Sf(:, :, t);
    Ct = Cf(:, :, t);
end
