function total_objective = objective_fun(params, concentration_data_a_6s, concentration_data_s_6s, concentration_data_c_6s)
% OBJECTIVE_FUN  Score a candidate parameter set for fmincon (see test.m).
% Runs the fisher4 reaction-diffusion simulation with the given
% parameters and returns the sum of squared differences between the
% simulated active-tip density and the target active-tip density
% (concentration_data_a_6s), plus a large penalty if any parameter falls
% outside a sensible range.

    dif = params(1);
    as = params(2);
    p_C = params(3);
    gcv = params(4);
    b_S= params(5);
    rb= params(6);
    rg= params(7);
    scale= params(8);

    total_objective = 0;

    % Only the active-tip density is currently used as the fitting
    % target; concentration_data_s_6s and concentration_data_c_6s are
    % accepted as inputs but not compared here.
    concentration_data = {concentration_data_a_6s};
    concentration = cell2mat(concentration_data);

    % Run the PDE simulation with the candidate parameters
    predicted_concentration = fisher4(dif, as, p_C, gcv, b_S, rb, rg, scale);

    % Sum of squared differences between predicted and target density
    squared_diff = sum(sum((predicted_concentration - concentration).^2));

    % Penalize parameter combinations outside sensible physical ranges
    if (as < 600) || (rb < 0.1) || (b_S > 1. || (rg < 0.1))
        squared_diff = squared_diff + 1000000000.;
    end

    total_objective = total_objective + squared_diff;
    disp(['total_objective:', num2str(total_objective)]);
end
