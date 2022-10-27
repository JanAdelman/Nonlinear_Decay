function fig4ef

% options 
write = true; % set to false if no output should be written 

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 5; % cell diameter [µm]
CV = 0; % CV for the kinetic parameters 
CV_area = 0; % CV for the areas 
CV_bc = logspace(log10(0.01), log10(1), 20); % define amout of variation in the flux/ amplitude values 
mu_lambda = 20; % gradient decay lenght in case of linear decay 
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
ncP = 200; % number of cells in the patterning domain
LP = ncP * diameter; % pattern length
powers = [1, 2];
readout_position = linspace(0, LP, 100); % readout positions
final_readout_positions = [5, 50, 150]; % for plotting 

% Dimensionless flux/ amplitude
c_ref = 1;
j_ref = 1;
a = 1;

mu_j = a * sqrt(mu_D*mu_d)/j_ref;
mu_c0 = a/c_ref;

% standard error 
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

% steady state solutions 
C_dirichlet = @(x, c0) c0*exp(-x/mu_lambda);
C_flux = @(x, j_0) (j_0*mu_lambda)/mu_D*exp(-x/mu_lambda);

if not(isfolder('bc_variability'))
    mkdir('bc_variability')
end
if not(isfolder('bc_variability/flux_variability/'))
    mkdir('bc_variability/flux_variability/')
end
if not(isfolder('bc_variability/c0_variability/'))
    mkdir('bc_variability/c0_variability/')
end

C_0_dirichlet = NaN(length(powers), 1);
C_0_flux = NaN(length(powers), 1);

% loop over the different bc 
for i = 1:length(powers)

    n = powers(i);
    
    % get domain 
    [~, l_p] = helper_functions.build_domain(0, LP, diameter, CV_area);

    % initialise the solver
    x0 = [];
    x0 = [x0, 0, l_p];
    x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

    nc = length(l_p);

    % get the solution for dirichlet BC 
    sol_dirichlet = solve_ode_dirichlet(x0, nc, CV, 0, n, tol, mu_d, mu_D, mu_c0);
    sol_flux = solve_ode_flux_bc(x0, nc, CV, 0, n, tol, mu_d, mu_D, mu_j);

    % get the concentration at the start and the end of the domain 
    C_0_dirichlet(i) = pchip(unique(sol_dirichlet.x, 'stable'), unique(sol_dirichlet.y(1,:),'stable'), 0);
    C_0_flux(i) = pchip(unique(sol_flux.x, 'stable'), unique(sol_flux.y(1,:),'stable'), 0);

end
    
% add noise
CV = 0.3;
CV_area = 0.5;
    
% loop over n
for i = 1:length(powers)
     
     % allocate memory 
     interp_readout_positions_5_cells_flux = NaN(length(CV_bc), 3); 
     interp_readout_positions_5_cells_dirichlet = NaN(length(CV_bc), 3);
     interp_readout_positions_50_cells_flux = NaN(length(CV_bc), 3); 
     interp_readout_positions_50_cells_dirichlet = NaN(length(CV_bc), 3);
     interp_readout_positions_150_cells_flux = NaN(length(CV_bc), 3); 
     interp_readout_positions_150_cells_dirichlet = NaN(length(CV_bc), 3);
     
     n = powers(i);
        
     for bc_val = 1:length(CV_bc)
         
            % filename
            filename_dirichlet = ['bc_variability/c0_variability/non_linear_decay_dirichlet_bc_' num2str(mu_c0) '_' num2str(n) '_' num2str(CV_bc(bc_val)) '.csv'];
            filename_flux = ['bc_variability/flux_variability/non_linear_decay_flux_bc_' num2str(mu_j) '_' num2str(n) '_' num2str(CV_bc(bc_val)) '.csv'];

            % linear decay, get readout concentrations along the domain 
            if n == 1          
                 K_dirichlet = C_dirichlet(readout_position, mu_c0);
                 K_flux = C_flux(readout_position, mu_j);

            % non-linear decay, use steady state solution for non-linear decay 
            % to find concentrations along the domain 
            else
               K_dirichlet = helper_functions.get_readout_conc_non_linear(readout_position, n, C_0_dirichlet(i), mu_lambda, c_ref);
               K_flux = helper_functions.get_readout_conc_non_linear(readout_position, n, C_0_flux(i), mu_lambda, c_ref);
            end 

            % allocate memory to store the readoutu positions 
            x_average_dirichlet = NaN(nruns, length(K_dirichlet));
            x_average_flux = NaN(nruns, length(K_flux));

            % array to store diameters 
            diam = NaN(nruns, 1);

            % loop over independent runs 
            for j = 1:nruns

                    % get domain (only patterning domain needed) 
                    [~, l_p] = helper_functions.build_domain(0, LP, diameter, CV_area);

                    % initialise the solver
                    x0 = [];
                    x0 = [x0, 0, l_p];
                    x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes

                    % get the number of cells 
                    nc = length(l_p);

                    domain = [0,l_p];
       
                    % solve the diffusion eq. with dirichlet BC 
                    sol_dirichlet = solve_ode_dirichlet(x0, nc, CV, CV_bc(bc_val), n, tol, mu_d, mu_D, mu_c0);

                    % array to store numerical results
                    y_sol_average_dirichlet =  NaN(length(nc), 1);                   
               
                    % get the average solution per cell 
                    y_sol_average_dirichlet = helper_functions.average_concentration(sol_dirichlet.x, sol_dirichlet.y(1, :), domain, y_sol_average_dirichlet, nc);
                
                    % solve the diffusion eq. with flux BC 
                    sol_flux = solve_ode_flux_bc(x0, nc, CV, CV_bc(bc_val), n, tol, mu_d, mu_D, mu_j);

                    % array to store numerical results
                    y_sol_average_flux = NaN(length(nc), 1);   
                                    
                    % get the average solution per cell 
                    y_sol_average_flux = helper_functions.average_concentration(sol_flux.x, sol_flux.y(1, :), domain, y_sol_average_flux, nc);

                    % calculate average
                    diam(j, 1) = mean(diff(domain));
          
                    % find the position where the threshold concentration is reached                 
                    % find the index where the concentration treshold is
                    % passed.

                    x_average_dirichlet(j, :) = helper_functions.getindex(y_sol_average_dirichlet, K_dirichlet, domain);
                    x_average_flux(j, :) = helper_functions.getindex(y_sol_average_flux, K_flux, domain);


            end

            % set all negative values to NaN (they lie outside of the
            % patterning domain) 
            x_average_dirichlet(x_average_dirichlet<0) = NaN;
            x_average_flux(x_average_flux<0) = NaN;

            % discharge all rows with NaN values 
            x_average_dirichlet=x_average_dirichlet(:,sum(isnan(x_average_dirichlet),1)==0); 
            x_average_flux=x_average_flux(:,sum(isnan(x_average_flux),1)==0); 

            average_diam = mean(diam); 
            
            % stats for dirichlet BC
            mean_pos_average_dirichlet = nanmean(x_average_dirichlet)/average_diam;
            std_pos_average_dirichlet = nanstd(x_average_dirichlet)/average_diam;
            SE_pos_average_dirichlet = nanstd(bootstrp(nboot, SEfun, x_average_dirichlet))/average_diam;

            % stats for flux BC
            mean_pos_average_flux = nanmean(x_average_flux)/average_diam;
            std_pos_average_flux = nanstd(x_average_flux)/average_diam;
            SE_pos_average_flux = nanstd(bootstrp(nboot, SEfun, x_average_flux))/average_diam;

            % write data to file 
            if write == true
                writetable(table(mean_pos_average_dirichlet', std_pos_average_dirichlet', SE_pos_average_dirichlet', 'VariableNames', {'mean_pos', 'std_pos', 'SE_std'}), filename_dirichlet);
                writetable(table(mean_pos_average_flux', std_pos_average_flux', SE_pos_average_flux', 'VariableNames', {'mean_pos', 'std_pos', 'SE_std'}), filename_flux);
            end 
            
             % get unique valus for interpolation 
            [unique_positions_flux, ind_flux, ~] = unique(mean_pos_average_flux, 'stable'); 
            std_pos_average_flux = std_pos_average_flux(ind_flux);
            SE_pos_average_flux = SE_pos_average_flux(ind_flux);
            
            [unique_positions_dirichlet, ind_dirichlet, ~] = unique(mean_pos_average_dirichlet, 'stable'); 
            std_pos_average_dirichlet = std_pos_average_dirichlet(ind_dirichlet);
            SE_pos_average_dirichlet = SE_pos_average_dirichlet(ind_dirichlet);
                        
            interp_std_flux = pchip(unique_positions_flux, std_pos_average_flux, final_readout_positions);
            interp_Se_flux = pchip(unique_positions_flux, SE_pos_average_flux, final_readout_positions);
            interp_std_dirichlet = pchip(unique_positions_dirichlet, std_pos_average_dirichlet, final_readout_positions);
            interp_Se_dirichelt = pchip(unique_positions_dirichlet, SE_pos_average_dirichlet, final_readout_positions);
            
            interp_readout_positions_5_cells_flux(bc_val, :) = [CV_bc(bc_val), interp_std_flux(1), interp_Se_flux(1)];
            interp_readout_positions_50_cells_flux(bc_val, :) = [CV_bc(bc_val), interp_std_flux(2), interp_Se_flux(2)];
            interp_readout_positions_150_cells_flux(bc_val, :) = [CV_bc(bc_val), interp_std_flux(3), interp_Se_flux(3)];
            
            interp_readout_positions_5_cells_dirichlet(bc_val, :) = [CV_bc(bc_val), interp_std_dirichlet(1), interp_Se_dirichelt(1)];
            interp_readout_positions_50_cells_dirichlet(bc_val, :) = [CV_bc(bc_val), interp_std_dirichlet(2), interp_Se_dirichelt(2)];
            interp_readout_positions_150_cells_dirichlet(bc_val, :) = [CV_bc(bc_val), interp_std_dirichlet(3), interp_Se_dirichelt(3)];

     end
   
     filename_dirichlet_interp_5 = ['bc_variability/non_linear_decay_dirichlet_bc_'  num2str(mu_c0) '_' num2str(n) '_readout_five_cells.csv'];
     filename_dirichlet_interp_50 = ['bc_variability/non_linear_decay_dirichlet_bc_'  num2str(mu_c0) '_' num2str(n) '_readout_fifty_cells.csv'];
     filename_dirichlet_interp_150 = ['bc_variability/non_linear_decay_dirichlet_bc_'  num2str(mu_c0) '_' num2str(n) '_readout_hundred_fifty_cells.csv'];
     
     filename_flux_interp_5 = ['bc_variability/non_linear_decay_flux_bc_'  num2str(mu_j) '_' num2str(n) '_readout_five_cells.csv'];
     filename_flux_interp_50 = ['bc_variability/non_linear_decay_flux_bc_'  num2str(mu_j) '_' num2str(n) '_readout_fifty_cells.csv'];
     filename_flux_interp_150 = ['bc_variability/non_linear_decay_flux_bc_'  num2str(mu_j) '_' num2str(n) '_readout_hundred_fifty_cells.csv'];

     writetable(table(interp_readout_positions_5_cells_dirichlet(:, 1), interp_readout_positions_5_cells_dirichlet(:, 2), interp_readout_positions_5_cells_dirichlet(:, 3), 'VariableNames', {'CV', 'std', 'SE'}), filename_dirichlet_interp_5); 
     writetable(table(interp_readout_positions_50_cells_dirichlet(:, 1), interp_readout_positions_50_cells_dirichlet(:, 2), interp_readout_positions_50_cells_dirichlet(:, 3), 'VariableNames', {'CV', 'std', 'SE'}), filename_dirichlet_interp_50); 
     writetable(table(interp_readout_positions_150_cells_dirichlet(:, 1), interp_readout_positions_150_cells_dirichlet(:, 2), interp_readout_positions_150_cells_dirichlet(:, 3), 'VariableNames', {'CV', 'std', 'SE'}), filename_dirichlet_interp_150); 
    
     writetable(table(interp_readout_positions_5_cells_flux(:, 1), interp_readout_positions_5_cells_flux(:, 2), interp_readout_positions_5_cells_flux(:, 3), 'VariableNames', {'CV', 'std', 'SE'}), filename_flux_interp_5); 
     writetable(table(interp_readout_positions_50_cells_flux(:, 1), interp_readout_positions_50_cells_flux(:, 2), interp_readout_positions_50_cells_flux(:, 3), 'VariableNames', {'CV', 'std', 'SE'}), filename_flux_interp_50); 
     writetable(table(interp_readout_positions_150_cells_flux(:, 1), interp_readout_positions_150_cells_flux(:, 2), interp_readout_positions_150_cells_flux(:, 3), 'VariableNames', {'CV', 'std', 'SE'}), filename_flux_interp_150); 

end 

%% functions for the ODE

% function to solve the ODE with flux BC 
function sol_flux = solve_ode_flux_bc(x0, nc, CV, CV_bc, n, tol, mu_d, mu_D, mu_j)

    options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);

    % default: all parameters constant
    d = mu_d * ones(nc, 1);
    D = mu_D * ones(nc, 1);

    % draw random kinetic parameters for each cell
    d = random(helper_functions.logndist(mu_d, CV), nc, 1); 
    D = random(helper_functions.logndist(mu_D, CV), nc, 1);

    % add noise the the flux 
    j_noisy =  random(helper_functions.logndist(mu_j, CV_bc), 1, 1);

    % get initial solution 
    sol0 = bvpinit(x0, @helper_functions.y0_non_lin);
    
    % get boundary conditions 
    odefun_bc = @(ya, yb) helper_functions.bcfun_flux(ya, yb, nc, j_noisy);

    odefun_init = @(x,y,c) helper_functions.odefun_non_lin_no_source(x, y, c, n, D, d);

    % solve the equation
    sol_flux = bvp4c(odefun_init, odefun_bc, sol0, options);
    
end 


% function to solve the ODE with Dirichlet BC 
function sol_dirichlet = solve_ode_dirichlet(x0, nc, CV, CV_bc, n, tol, mu_d, mu_D, mu_c0)

    options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);

    % default: all parameters constant
    d = mu_d * ones(nc, 1);
    D = mu_D * ones(nc, 1);

    % draw random kinetic parameters for each cell
    d = random(helper_functions.logndist(mu_d, CV), nc, 1); 
    D = random(helper_functions.logndist(mu_D, CV), nc, 1);
    
    % add noise to the amplitude 
    dir_c0 = random(helper_functions.logndist(mu_c0, CV_bc), 1, 1);
     
    % get initial solution 
    sol0 = bvpinit(x0, @helper_functions.y0_non_lin);

    odefun_init = @(x,y,c) helper_functions.odefun_non_lin_no_source(x, y, c, n, D, d);    
    odefun_bc = @(ya, yb) helper_functions.bcfun_dirichlet(ya, yb, nc, dir_c0);

    % solve the equation
    sol_dirichlet = bvp4c(odefun_init, odefun_bc, sol0, options);
    
end 

end 
