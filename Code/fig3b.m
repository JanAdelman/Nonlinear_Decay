function fig3b

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 5; % cell diameter [µm]
CV = 0; % coefficient of variation for the kinetic parameters 
CV_A = 0; % area variability
mu_lambda_unscaled = 20; % mean exponential gradient decay length [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_d = mu_D/mu_lambda_unscaled^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d; % mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 100; % number of cells in the patterning domain
LS = ncS * diameter; % source length
LP = ncP * diameter; % pattern length
final_readout_positions = [5, 10, 15]; % for plotting 
readout_pos_name = ["five_lambda", "ten_lambda", "fifteen_lambda"];
C_ref = 1; % reference concentration 
diameters = linspace(mu_lambda_unscaled/40, 2*mu_lambda_unscaled, 20); % define a range of varying diameters
powers = [1,2,4]; 

% analytical deterministic solution
C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda_unscaled)) + sinh(LS/mu_lambda_unscaled) / sinh((LS+LP)/mu_lambda_unscaled) * cosh((LP-x)/mu_lambda_unscaled));

% standard error 
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

C_0 = NaN(length(powers), 1);

dir = 'fig3b';
if not(isfolder(dir))
    mkdir(dir)
end

for i = 1:length(powers)

    n = powers(i);
    
    % get domain 
    [l_s, l_p] = helper_functions.build_domain(LS, LP, diameter, CV_A);
    
    % initialise the solver
    x0 = [];
    x0 = [x0, -l_s, 0, l_p];
    x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

    nc = length(l_p) + length(l_s);
    ncS = length(l_s);

    sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p, mu_d, mu_D);
    
    % get the concentration at the start of the domain 
    C_0(i) = pchip(unique(sol.x, 'stable'), unique(sol.y(1,:), 'stable'), 0);

end 

% add noise
CV = 0.3;
CV_A = 0.5;

% correction of the mean lambda (see https://doi.org/10.1038/s41467-022-28834-3)
mu_lambda = mu_lambda_unscaled/(1 - 0.011 * CV + 1.355 * CV^2 - 0.179 * CV^3 + 0.0077 * CV^4)^-0.357;

% get readout positions:
readout_position = linspace(0, LP, 100);

% loop over n    
for i = 1:length(powers)

    % get the power of the nonlinear decay 
    n = powers(i);

    % linear decay, get readout concentrations along the domain 
    if n == 1          
       K = C(readout_position);

    % non-linear decay, use steady-state solution for non-linear decay 
    % to find concentrations along the domain 
    else
        K = helper_functions.get_readout_conc_non_linear(readout_position, n, C_0(i), mu_lambda_unscaled, C_ref);
    end 
    
    interp_results_five_lambda = NaN(length(diameters), 3);
    interp_results_ten_lambda = NaN(length(diameters), 3);
    interp_results_fifteen_lambda = NaN(length(diameters), 3);

    % loop over varying diameters 
    for diam = 1:length(diameters)
        
            diam_iter = NaN(nruns, 1);

            % array for saving average readout position per run 
            x_average = NaN(nruns, length(K));

            for j = 1:nruns

                    [l_s, l_p] = helper_functions.build_domain(LS, LP, diameters(diam), CV_A);

                    % initialise the solver
                    x0 = [];
                    x0 = [x0, -l_s, 0, l_p];
                    x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

                    nc = length(l_p) + length(l_s);
                    ncS = length(l_s);

                    domain = [-l_s,0,l_p];

                    sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p, mu_d, mu_D);

                    % array to store numerical integration results
                    y_sol_average = NaN(length(nc), 1);

                    % get the average solution per cell 
                    y_sol_average = helper_functions.average_concentration(sol.x, sol.y(1, :), domain, y_sol_average, nc);
                
                    diam_iter(j, 1) = mean(diff(domain));

                    % find the index where the concentration threshold is passed
                    x_average(j, :) = helper_functions.getindex(y_sol_average, K, domain);

            end

            % set all negative values to NaN (they lie outside of the
            % patterning domain) 
            x_average(x_average<0) = NaN;

            % discard all rows with NaN values 
            x_average = x_average(:,sum(isnan(x_average),1)==0); 

            % normalise by lambda 
            average_diam = mean(diam_iter)/mu_lambda;   
            mean_pos_average = nanmean(x_average)/mu_lambda;
            std_pos_average = nanstd(x_average)/mu_lambda;
            SE_pos_average = nanstd(bootstrp(nboot, SEfun, x_average))/mu_lambda;
            
            % get unique valus for interpolation 
            [unique_positions, ind, ~] = unique(mean_pos_average, 'stable');
            std_pos_average = std_pos_average(ind);
            SE_pos_average = SE_pos_average(ind);
            
            interp_std_pos = pchip(unique_positions, std_pos_average, final_readout_positions);
            interp_SE_pos= pchip(unique_positions, SE_pos_average, final_readout_positions);
            
            interp_results_five_lambda(diam, :) = [average_diam, interp_std_pos(1), interp_SE_pos(1)];
            interp_results_ten_lambda(diam, :) = [average_diam, interp_std_pos(2), interp_SE_pos(2)];
            interp_results_fifteen_lambda(diam, :) = [average_diam, interp_std_pos(3), interp_SE_pos(3)];
            
    end

    filename_five_lambda = [dir '/non_linear_decay_varying_diameters_' num2str(n) '_' num2str(readout_pos_name(1)) '.csv'];
    filename_ten_lambda = [dir '/non_linear_decay_varying_diameters_' num2str(n) '_' num2str(readout_pos_name(2)) '.csv'];
    filename_fifteen_lambda = [dir '/non_linear_decay_varying_diameters_' num2str(n) '_' num2str(readout_pos_name(3)) '.csv'];
    
    writetable(table(interp_results_five_lambda(:, 1), interp_results_five_lambda(:, 2), interp_results_five_lambda(:, 3), 'VariableNames', {'diameter', 'std', 'SE'}), filename_five_lambda);
    writetable(table(interp_results_ten_lambda(:, 1), interp_results_ten_lambda(:, 2), interp_results_ten_lambda(:, 3), 'VariableNames', {'diameter', 'std', 'SE'}), filename_ten_lambda);
    writetable(table(interp_results_fifteen_lambda(:, 1), interp_results_fifteen_lambda(:, 2), interp_results_fifteen_lambda(:, 3), 'VariableNames', {'diameter', 'std', 'SE'}), filename_fifteen_lambda);

end

%% functions for the ODE

function sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p, mu_d, mu_D)

    options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);

    % draw random kinetic parameters for each cell
    p = random(helper_functions.logndist(mu_p, CV), nc, 1);
    d = random(helper_functions.logndist(mu_d, CV), nc, 1);
    D = random(helper_functions.logndist(mu_D, CV), nc, 1);

    % get initial solution 
    sol0 = bvpinit(x0, @helper_functions.y0_non_lin);

    odefun_init = @(x,y,c) helper_functions.odefun_non_lin(x, y, c, n, D, p, d, ncS, C_ref);

    % solve the equation
    sol = bvp4c(odefun_init, @(ya, yb) helper_functions.bcfun(ya, yb, nc), sol0, options);
    
end 


end
