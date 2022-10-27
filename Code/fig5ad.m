function fig5ad

% options
write = true;

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 5; % cell diameter [µm]
mu_lambda = 20; % gradient decay length 
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
mu_p  = logspace(log10(0.001*mu_d), log10(10*mu_d), 20); % change source strength,  mean morphogen production rate [substance/(µm^3*s)]
ncS = 5; % number of cells in the source domain
ncP = 200;% number of cells in the patterning domain
LS = ncS * diameter; % source length
LP = ncP * diameter; % pattern length
final_readout_positions = [5, 50, 150]; % for plotting 
readout_position = linspace(0, LP, 100); % readout positions 
c_ref = 1; % reference concentration 
powers = [1,2,4];

% standard error 
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

if not(isfolder('production_change'))
    mkdir('production_change')
end

% loop over all n
for i = 1:length(powers)
    
    % allocate memory 
    interp_readout_positions_5_cells = NaN(length(mu_p), 3); 
    interp_readout_positions_50_cells = NaN(length(mu_p), 3);  
    interp_readout_positions_150_cells = NaN(length(mu_p), 3);
    
    % get the power of the nonlinear decay 
    n = powers(i);
    
    for p_val = 1:length(mu_p)

            % filename
            filename = ['production_change/non_linear_decay_' num2str(mu_p(p_val)/mu_d) '_' num2str(n) '_varying_p.csv'];
           
            CV = 0;
            CV_area = 0;

            % get domain 
            [l_s, l_p] = helper_functions.build_domain(LS, LP, diameter, CV_area);
            
            % initialise the solver
            x0 = [];
            x0 = [x0, -l_s, 0, l_p];
            x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes
        
            nc = length(l_p) + length(l_s);
            ncS = length(l_s);
        
            sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p(p_val), mu_d, mu_D);
            
            % analytical deterministic solution
            C = @(x, LP) mu_p(p_val)/mu_d * ((x<0) .* (1-cosh(x/mu_lambda)) + sinh(LS/mu_lambda) / sinh((LS+LP)/mu_lambda) * cosh((LP-x)/mu_lambda));
            
            % linear decay, get readout concentrations along the domain 
            if n == 1          
                 K = C(readout_position, LP);

            % non-linear decay, use steady state solution for non-linear decay 
            % to find concentrations along the domain 
            else

                % get the concentration at the start of the patterning domain 
                C_0 = pchip(unique(sol.x, 'stable'), unique(sol.y(1,:),'stable'), 0);
               
                K = helper_functions.get_readout_conc_non_linear(readout_position, n, C_0, mu_lambda, c_ref);

            end 

            diam = NaN(nruns, 1);

            % add noise
            CV = 0.3;
            CV_area = 0.5;

            % array for saving average readout position per run 
            x_average = NaN(nruns, length(K));
            
            for j = 1:nruns

                    [l_s, l_p] = helper_functions.build_domain(LS, LP, diameter, CV_area);
                    
                    % initialise the solver
                    x0 = [];
                    x0 = [x0, -l_s, 0, l_p];
                    x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

                    nc = length(l_p) + length(l_s);
                    ncS = length(l_s);

                    domain = [-l_s,0,l_p];

                    sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p(p_val), mu_d, mu_D);
 
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

            % discharge all rows with NaN values 
            x_average=x_average(:,sum(isnan(x_average),1)==0); 

            average_diam = mean(diam);   

            mean_pos_average = nanmean(x_average)/average_diam;
            std_pos_average = nanstd(x_average)/average_diam;
            SE_pos_average = nanstd(bootstrp(nboot, SEfun, x_average))/average_diam;
                        
            if write == true
                writetable(table(mean_pos_average', std_pos_average', SE_pos_average', 'VariableNames', {'mean_pos', 'std_pos', 'SE_std'}), filename);
            end 
            
             % get unique valus for interpolation 
            [unique_positions, ind, ~] = unique(mean_pos_average, 'stable'); 
            std_pos_average = std_pos_average(ind);
            SE_pos_average = SE_pos_average(ind);
            
            interp_std = pchip(unique_positions, std_pos_average, final_readout_positions);
            interp_Se = pchip(unique_positions, SE_pos_average, final_readout_positions);
            
            interp_readout_positions_5_cells(p_val, :) = [mu_p(p_val)/mu_d, interp_std(1), interp_Se(1)];
            interp_readout_positions_50_cells(p_val, :) = [mu_p(p_val)/mu_d, interp_std(2), interp_Se(2)];
            interp_readout_positions_150_cells(p_val, :) = [mu_p(p_val)/mu_d, interp_std(3), interp_Se(3)];
    
    end
     
     filename_interp_5 = ['production_change/non_linear_decay_'  num2str(n) '_readout_five_cells_varying_p.csv'];
     filename_interp_50 = ['production_change/non_linear_decay_'  num2str(n) '_readout_fifty_cells_varying_p.csv'];
     filename_interp_150 = ['production_change/non_linear_decay_'  num2str(n) '_readout_hundred_fifty_cells_varying_p.csv'];
    
     writetable(table(interp_readout_positions_5_cells(:, 1), interp_readout_positions_5_cells(:, 2), interp_readout_positions_5_cells(:, 3), 'VariableNames', {'p/d', 'std', 'SE'}), filename_interp_5); 
     writetable(table(interp_readout_positions_50_cells(:, 1), interp_readout_positions_50_cells(:, 2), interp_readout_positions_50_cells(:, 3), 'VariableNames', {'p/d', 'std', 'SE'}), filename_interp_50); 
     writetable(table(interp_readout_positions_150_cells(:, 1), interp_readout_positions_150_cells(:, 2), interp_readout_positions_150_cells(:, 3), 'VariableNames', {'p/d', 'std', 'SE'}), filename_interp_150); 
    
end 
%% functions for the ODE

function sol = solve_ode(x0, nc, ncS, n, CV, tol, mu_p, mu_d, mu_D)

    options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);

    % default: all parameters constant
    p = mu_p * ones(nc, 1);
    d = mu_d * ones(nc, 1);
    D = mu_D * ones(nc, 1);

    % draw random kinetic parameters for each cell

    p = random(helper_functions.logndist(mu_p, CV), nc, 1);
    d = random(helper_functions.logndist(mu_d, CV), nc, 1); 
    D = random(helper_functions.logndist(mu_D, CV), nc, 1);

    % get initial solution 
    sol0 = bvpinit(x0, @helper_functions.y0_non_lin);

    odefun_init = @(x,y,c) helper_functions.odefun_non_lin(x, y, c, n, D, p, d, ncS);

    % solve the equation
    sol = bvp4c(odefun_init, @(ya, yb) helper_functions.bcfun(ya, yb, nc), sol0, options);
    
end 

end
