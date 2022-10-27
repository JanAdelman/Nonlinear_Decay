function fig3a

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 5; % cell diameter [µm]
CV = 0; % CV for the kinetic parameters 
CV_area = 0; % area variability
mu_lambda = 20; % gradient decay length 
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p = mu_d;% mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 200; % number of cells in the patterning domain
LS = ncS * diameter; % source length
LP = ncP * diameter; % pattern length
powers = [1,2,4];
c_ref = 1; % reference concentration 
readout_position = linspace(0, LP, 100); % define readout positions along patterning domain 
final_readout_positions = [5,50,150]; % for plotting 

% analytical deterministic solution
C = @(x) mu_p/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));

% standard errror 
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

C_0 = NaN(length(powers), 1);

if not(isfolder('non_linear_decay_varying_cv_area'))
    mkdir('non_linear_decay_varying_cv_area')
end

% get c_0 for the analytical solution in the noise free case
for i = 1:length(powers)

    n = powers(i);

    % get domain 
    [l_s, l_p] = helper_functions.build_domain(LS, LP, diameter, CV_area);
    
    % initialise the solver
    x0 = [];
    x0 = [x0, -l_s, 0, l_p];
    x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

    nc = length(l_p) + length(l_s);
    ncS = length(l_s);

    sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p, mu_d, mu_D);

    % get the concentration at the start of the domain 
    C_0(i) = pchip(unique(sol.x, 'stable'), unique(sol.y(1,:),'stable'), 0);

end 

% add noise
CV = 0.3;
CV_area = logspace(log10(0.01), log10(10), 20);

% loop over all n
for i = 1:length(powers)

    interp_std = NaN(length(CV_area), 3);
    interp_Se = NaN(length(CV_area), 3);

    filename_readout_positions = ['non_linear_decay_varying_cv_area/non_linear_decay_readout_positions_' num2str(n)  '.csv'];
    
    %loop over the area variabilities 
    for area_var = 1:length(CV_area)
        
        % get the power of the nonlinear decay 
        n = powers(i);

        % filename
        filename = ['non_linear_decay_varying_cv_area/non_linear_decay_'  num2str(CV_area(area_var))  '_' num2str(n)  '.csv'];
               
        % linear decay, get readout concentrations along the domain 
        if n == 1          
             K = C(readout_position);

        % non-linear decay, use steady state solution for non-linear decay 
        % to find concentrations along the domain 
        else
           K = helper_functions.get_readout_conc_non_linear(readout_position, n, C_0(i), mu_lambda, c_ref);

        end 

        diam = NaN(nruns, 1);

        % array for saving average readout position per run 
        x_average = NaN(nruns, length(K));

        for j = 1:nruns
            
                [l_s, l_p] = helper_functions.build_domain(LS, LP, diameter, CV_area(area_var));
                
                % initialise the solver
                x0 = [];
                x0 = [x0, -l_s, 0, l_p];
                x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

                nc = length(l_p) + length(l_s);
                ncS = length(l_s);

                domain = [-l_s,0,l_p];

                sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p, mu_d, mu_D);

                % array to store numerical integration results
                y_sol_average =  NaN(length(nc), 1);

                % get the average solution per cell 
                y_sol_average = helper_functions.average_concentration(sol.x, sol.y(1, :), domain, y_sol_average, nc);

                diam(j, 1) = mean(diff(domain));

                % find the position where the threshold concentration is reached                 
                % find the index where the concentration treshold is
                % passed. 

                x_average(j, :) = helper_functions.getindex(y_sol_average, K, domain);
                
        end

        % set all negative values to NaN (they lie outside of the
        % patterning domain) 
        x_average(x_average<0) = NaN;

        % save the average diameter 
        average_diam = mean(diam);  
        
        % save the average readout position 
        mean_pos_average = nanmean(x_average)/average_diam;

        % save the positional error 
        std_pos_average = nanstd(x_average)/average_diam;
        
        % save the standard error of the positional error 
        SE_pos_average = nanstd(bootstrp(nboot, SEfun, x_average))/average_diam;

         % get unique values for interpolation 
        [unique_positions, indx, ~] = unique(mean_pos_average, 'stable'); 
        std_pos_average = std_pos_average(indx);
        SE_pos_average = SE_pos_average(indx);
        
        % interpolate to get final readout positions 
        interp_std(area_var, :) = pchip(unique_positions, std_pos_average, final_readout_positions);
        interp_Se(area_var, :) = pchip(unique_positions, SE_pos_average, final_readout_positions);
        
        writetable(table(mean_pos_average', std_pos_average', SE_pos_average', 'VariableNames', {'mean_pos', 'std_pos', 'SE_std'}), filename);

    end

    writetable(table(CV_area', interp_std(:, 1), interp_Se(:, 1), interp_std(:, 2), interp_Se(:, 2), interp_std(:, 3), interp_Se(:, 3), ...
    'VariableNames', {'CV_area', 'std_five_cells', 'SE_five_cells', 'std_fifty_cells', 'SE_fifty_cells', 'std_hundred_fifty_cells','SE_hundred_fifty_cells'}), ...
    filename_readout_positions);


end 

%% functions for the ODE

function sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p, mu_d, mu_D, k)

    options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);

    % default: all parameters constant
    p = mu_p * ones(nc, 1);
    d = mu_d * ones(nc, 1);
    D = mu_D * ones(nc, 1);
    
    % draw random kinetic parameters for each cell
    p = random(helper_functions.logndist(mu_p, CV), nc, 1);
    d = random(helper_functions.logndist(mu_d, CV), nc, 1);
    D = random(helper_functions.logndist(mu_D,  CV), nc, 1);

    % get initial solution 
    sol0 = bvpinit(x0, @helper_functions.y0_non_lin);

    odefun_init = @(x,y,c) helper_functions.odefun_non_lin(x, y, c, n, D, p, d, ncS);

    % solve the equation
    sol = bvp4c(odefun_init, @(ya, yb) helper_functions.bcfun(ya, yb, nc), sol0, options);
    
end 

end
