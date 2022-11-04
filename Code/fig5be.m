function fig5be

% options 
write = true; % set to false if no output should be written 

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 5; % cell diameter [µm]
mu_lambda = 20; % mean exponential gradient decay length [µm]
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
ncP = 200; % number of cells in the patterning domain
LP = ncP * diameter; % patterning domain length
C_ref = 1; % reference concentration
j_ref = mu_D/mu_lambda * C_ref; % reference influx
final_readout_positions = [5, 50, 150]; % for plotting 
powers = [1,2,4];

% standard error 
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

% deterministic solution
C_flux = @(x, j0) j0*mu_lambda/mu_D*exp(-x/mu_lambda);

dir = 'fig5be';
if not(isfolder(dir))
    mkdir(dir)
end
if not(isfolder([dir '/flux_change/']))
    mkdir([dir '/flux_change/'])
end

%% Dirichlet boundary conditions & Neumann boundary conditions at x=0

% get readout positions:
readout_position = linspace(0, LP, 100);
   
% loop over n
for i = 1:length(powers)
     
     n = powers(i);

     % due to flat gradients, use different parameters in the case n=4
     if n == 4
         mu_j = logspace(-2, 1, 21) * j_ref; 
     else
         mu_j = logspace(-3, 2, 31) * j_ref; 
     end 
       
     % allocate memory 
     interp_readout_positions_5_cells_flux = NaN(length(mu_j), 3);
     interp_readout_positions_50_cells_flux = NaN(length(mu_j), 3); 
     interp_readout_positions_150_cells_flux = NaN(length(mu_j), 3); 
        
     for bc_val = 1:length(mu_j)

            CV_A = 0;
    
            % linear decay, get readout concentrations along the domain 
            if n == 1      
    
                  K_flux = C_flux(readout_position, mu_j(bc_val));
     
            % non-linear decay, use steady-state solution for non-linear decay 
            % to find concentrations along the domain 
            else
    
                % get domain 
                [~, l_p] = helper_functions.build_domain(0, LP, diameter, CV_A);
            
                % initialise the solver
                x0 = [];
                x0 = [x0, 0, l_p];
                x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes
            
                nc = length(l_p);
            
                % get the solution for Neumann BC 
                sol_flux = solve_ode_flux(x0, nc, 0, 0, n, tol, mu_d, mu_D, mu_j(bc_val));
            
                % get the concentration at the start and the end of the domain 
                C_0_flux = pchip(unique(sol_flux.x, 'stable'), unique(sol_flux.y(1,:), 'stable'), 0);
    
                K_flux = helper_functions.get_readout_conc_non_linear(readout_position, n, C_0_flux, mu_lambda, C_ref);
                      
            end 

            % filename
            filename_flux = [dir '/flux_change/non_linear_decay_flux_bc_' num2str(mu_j(bc_val)/j_ref) '_' num2str(n) '.csv'];
            
            % allocate memory to store the readoutu positions 
            x_average_flux = NaN(nruns, length(readout_position));

            % array to store diameters 
            diam = NaN(nruns, 1);

            % add noise
            CV = 0.3;
            CV_A = 0.5;
            CV_bc = 0.3; 

            % loop over independent runs 
            for j = 1:nruns

                    % get domain (only patterning domain needed) 
                    [~, l_p] = helper_functions.build_domain(0, LP, diameter, CV_A);

                    % initialise the solver
                    x0 = [];
                    x0 = [x0, 0, l_p];
                    x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

                    % get the number of cells 
                    nc = length(l_p);

                    domain = [0,l_p];
                   
                    % solve the diffusion eq. with Neumann BC 
                    sol_flux = solve_ode_flux(x0, nc, CV, CV_bc, n, tol, mu_d, mu_D, mu_j(bc_val));

                    % array to store numerical results
                    y_sol_average_flux = NaN(length(nc), 1);
                                    
                    % get the average solution per cell 
                    y_sol_average_flux = helper_functions.average_concentration(sol_flux.x, sol_flux.y(1, :), domain, y_sol_average_flux, nc);

                    % calculate average
                    diam(j, 1) = mean(diff(domain));

                    % find the index where the concentration threshold is passed
                    x_average_flux(j, :) = helper_functions.getindex(y_sol_average_flux, K_flux, domain);

            end

            % set all negative values to NaN (they lie outside of the
            % patterning domain) 
            x_average_flux(x_average_flux<0) = NaN;

            % discard all rows with NaN values 
            x_average_flux = x_average_flux(:,sum(isnan(x_average_flux),1)==0);

            average_diam = mean(diam);
            
            % stats for flux BC            
            mean_pos_average_flux = nanmean(x_average_flux)/average_diam;
            std_pos_average_flux = nanstd(x_average_flux)/average_diam;
            SE_pos_average_flux = nanstd(bootstrp(nboot, SEfun, x_average_flux))/average_diam;

            % write data to file 
            if write == true
                 writetable(table(mean_pos_average_flux', std_pos_average_flux', SE_pos_average_flux', 'VariableNames', {'mean_pos', 'std_pos', 'SE_std'}), filename_flux);
            end 
            
            % get unique valus for interpolation 
            [unique_positions_flux, ind_flux, ~] = unique(mean_pos_average_flux, 'stable');
            std_pos_average_flux = std_pos_average_flux(ind_flux);
            SE_pos_average_flux = SE_pos_average_flux(ind_flux);
              
            interp_std_flux = pchip(unique_positions_flux, std_pos_average_flux, final_readout_positions);
            interp_SE_flux = pchip(unique_positions_flux, SE_pos_average_flux, final_readout_positions);
            
            interp_readout_positions_5_cells_flux(bc_val, :) = [mu_j(bc_val)/j_ref, interp_std_flux(1), interp_SE_flux(1)];
            interp_readout_positions_50_cells_flux(bc_val, :) = [mu_j(bc_val)/j_ref, interp_std_flux(2), interp_SE_flux(2)];
            interp_readout_positions_150_cells_flux(bc_val, :) = [mu_j(bc_val)/j_ref, interp_std_flux(3), interp_SE_flux(3)];
     
     end
      
     filename_flux_interp_5 = [dir '/flux_change_' num2str(n) '_readout_five_cells.csv'];
     filename_flux_interp_50 = [dir '/flux_change_' num2str(n) '_readout_fifty_cells.csv'];
     filename_flux_interp_150 = [dir '/flux_change_' num2str(n) '_readout_hundred_fifty_cells.csv'];
    
     writetable(table(interp_readout_positions_5_cells_flux(:, 1), interp_readout_positions_5_cells_flux(:, 2), interp_readout_positions_5_cells_flux(:, 3), 'VariableNames', {'mu_j', 'std', 'SE'}), filename_flux_interp_5);
     writetable(table(interp_readout_positions_50_cells_flux(:, 1), interp_readout_positions_50_cells_flux(:, 2), interp_readout_positions_50_cells_flux(:, 3), 'VariableNames', {'mu_j', 'std', 'SE'}), filename_flux_interp_50);
     writetable(table(interp_readout_positions_150_cells_flux(:, 1), interp_readout_positions_150_cells_flux(:, 2), interp_readout_positions_150_cells_flux(:, 3), 'VariableNames', {'mu_j', 'std', 'SE'}), filename_flux_interp_150);

end 

%% functions for the ODE

% function to solve the ODE with flux BC 
function sol_flux = solve_ode_flux(x0, nc, CV, CV_bc, n, tol, mu_d, mu_D, mu_j)

    options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);

    % draw random kinetic parameters for each cell
    d = random(helper_functions.logndist(mu_d, CV), nc, 1); 
    D = random(helper_functions.logndist(mu_D, CV), nc, 1);

    dir_j = random(helper_functions.logndist(mu_j, CV_bc), 1, 1);

    % get initial solution 
    sol0 = bvpinit(x0, @helper_functions.y0_non_lin);
    
    % get boundary conditions 
    odefun_bc = @(ya, yb) helper_functions.bcfun_flux(ya, yb, nc, dir_j);
    odefun_init = @(x,y,c) helper_functions.odefun_non_lin_no_source(x, y, c, n, D, d, C_ref);

    % solve the equation
    sol_flux = bvp4c(odefun_init, odefun_bc, sol0, options);
    
end 


end 
