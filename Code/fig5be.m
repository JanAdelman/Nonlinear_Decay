function fig5be 

% options 
write = true; % set to false if no output should be written 

% parameters
tol = 1e-10; % numerical tolerance for solver and fitting
nruns = 1000; % number of independent simulation runs
nboot = 1e4; % number of bootstrap samples for error estimation
diameter = 5; % cell diameter [µm]
mu_lambda = 20; % gradient decay lenght in case of linear decay 
mu_D = 0.033; % mean morphogen diffusion constant [µm^2/s]
mu_d = mu_D/mu_lambda^2; % mean morphogen degradation rate [1/s]
ncP = 200; % number of cells in the patterning domain
LP = ncP * diameter; % pattern length
final_readout_positions = [5, 50, 150]; % for plotting 
c_ref = 1; % reference concentration 
readout_position = linspace(0, LP, 100); 
powers = [1,2,4];

% standard error 
SEfun = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

% steady state solutions 
C_dirichlet = @(x, c0) c0*exp(-x/mu_lambda);

if not(isfolder('bc_change'))
    mkdir('bc_change')
end

if not(isfolder('bc_change/c0_change/'))
    mkdir('bc_change/c0_change/')
end

% loop over n
for i = 1:length(powers)
     
     n = powers(i);

     % due to flat gradients/ singular jacobians, use different parameters in the case
     % n=4
     if n == 4
         mu_c0 = logspace(log10(0.4), log10(4), 10)/c_ref; % amplitude 
     else
         mu_c0 = logspace(log10(0.01), log10(100), 30)/c_ref;
     end 

     % allocate memory 
     interp_readout_positions_5_cells_dirichlet = NaN(length(mu_c0), 3);
     interp_readout_positions_50_cells_dirichlet = NaN(length(mu_c0), 3);
     interp_readout_positions_150_cells_dirichlet = NaN(length(mu_c0), 3);
        
     for bc_val = 1:length(mu_c0)
                
           CV_area = 0;

           % linear decay, get readout concentrations along the domain 
           if n == 1      
        
                K_dirichlet = C_dirichlet(readout_position, mu_c0(bc_val));
 
           % non-linear decay, use steady state solution for non-linear decay 
           % to find concentrations along the domain 
           else
  
                % get domain 
                [~, l_p] = helper_functions.build_domain(0, LP, diameter, CV_area);
            
                % initialise the solver
                x0 = [];
                x0 = [x0, 0, l_p];
                x0 = sort([x0 x0(2:end-1)]); % duplicate interface nodes
            
                nc = length(l_p);
            
                % get the solution for dirichlet BC 
                sol_dirichlet = solve_ode_dirichlet(x0, nc, 0, 0, n, tol, mu_d, mu_D, mu_c0(bc_val));
               
                % get the concentration at the start and the end of the domain 
                C_0_dirichlet = pchip(unique(sol_dirichlet.x, 'stable'), unique(sol_dirichlet.y(1,:),'stable'), 0);
             
                K_dirichlet = helper_functions.get_readout_conc_non_linear(readout_position, n, C_0_dirichlet, mu_lambda, c_ref);
                   
            end 
            
            % allocate memory to store the readoutu positions 
            x_average_dirichlet = NaN(nruns, length(readout_position));
          
            % array to store diameters 
            diam = NaN(nruns, 1);

            % filename
            filename_dirichlet = ['bc_change/c0_change/non_linear_decay_dirichlet_bc_' num2str(mu_c0(bc_val)) '_' num2str(n)  '.csv'];
            
            % add noise
            CV = 0.3;
            CV_area = 0.5;
            CV_bc = 0.3;

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
                sol_dirichlet = solve_ode_dirichlet(x0, nc, CV, CV_bc, n, tol, mu_d, mu_D, mu_c0(bc_val));
               
                % array to store numerical results
                y_sol_average_dirichlet =  NaN(length(nc), 1);  
                
                % get the average solution per cell 
                y_sol_average_dirichlet = helper_functions.average_concentration(sol_dirichlet.x, sol_dirichlet.y(1, :), domain, y_sol_average_dirichlet, nc);
           
                % calculate average diameter 
                diam(j, 1) = mean(diff(domain));
      
                % find the position where the threshold concentration is reached                 
                % find the index where the concentration treshold is
                % passed. (Beginning of a cell)

                x_average_dirichlet(j, :) = helper_functions.getindex(y_sol_average_dirichlet, K_dirichlet, domain);
             
            end
            
            % set all negative values to NaN (they lie outside of the
            % patterning domain) 
            x_average_dirichlet(x_average_dirichlet<0) = NaN;

            % discharge all rows with NaN values 
            x_average_dirichlet=x_average_dirichlet(:,sum(isnan(x_average_dirichlet),1)==0); 

            average_diam = mean(diam); 
            
            % stats for dirichlet BC
            mean_pos_average_dirichlet = nanmean(x_average_dirichlet)/average_diam;
            std_pos_average_dirichlet = nanstd(x_average_dirichlet)/average_diam;
            SE_pos_average_dirichlet = nanstd(bootstrp(nboot, SEfun, x_average_dirichlet))/average_diam;
         
            % write data to file 
            if write == true
                writetable(table(mean_pos_average_dirichlet', std_pos_average_dirichlet', SE_pos_average_dirichlet', 'VariableNames', {'mean_pos', 'std_pos', 'SE_std'}), filename_dirichlet);
            end 
            
             % get unique values for interpolation 
            [unique_positions_dirichlet, ind_dirichlet, ~] = unique(mean_pos_average_dirichlet, 'stable'); 
            std_pos_average_dirichlet = std_pos_average_dirichlet(ind_dirichlet);
            SE_pos_average_dirichlet = SE_pos_average_dirichlet(ind_dirichlet);
            
            interp_std_dirichlet = pchip(unique_positions_dirichlet, std_pos_average_dirichlet, final_readout_positions);
            interp_Se_dirichelt = pchip(unique_positions_dirichlet, SE_pos_average_dirichlet, final_readout_positions);
 
            interp_readout_positions_5_cells_dirichlet(bc_val, :) = [mu_c0(bc_val), interp_std_dirichlet(1), interp_Se_dirichelt(1)];
            interp_readout_positions_50_cells_dirichlet(bc_val, :) = [mu_c0(bc_val), interp_std_dirichlet(2), interp_Se_dirichelt(2)];
            interp_readout_positions_150_cells_dirichlet(bc_val, :) = [mu_c0(bc_val), interp_std_dirichlet(3), interp_Se_dirichelt(3)];

     end
   
     filename_dirichlet_interp_5 = ['bc_change/c0_change_'  num2str(n) '_readout_five_cells.csv'];
     filename_dirichlet_interp_50 = ['bc_change/c0_change_'   num2str(n) '_readout_fifty_cells.csv'];
     filename_dirichlet_interp_150 = ['bc_change/c0_change_'  num2str(n) '_readout_hundred_fifty_cells.csv'];
     
     writetable(table(interp_readout_positions_5_cells_dirichlet(:, 1), interp_readout_positions_5_cells_dirichlet(:, 2), interp_readout_positions_5_cells_dirichlet(:, 3), 'VariableNames', {'mu_c0', 'std', 'SE'}), filename_dirichlet_interp_5); 
     writetable(table(interp_readout_positions_50_cells_dirichlet(:, 1), interp_readout_positions_50_cells_dirichlet(:, 2), interp_readout_positions_50_cells_dirichlet(:, 3), 'VariableNames', {'mu_c0', 'std', 'SE'}), filename_dirichlet_interp_50); 
     writetable(table(interp_readout_positions_150_cells_dirichlet(:, 1), interp_readout_positions_150_cells_dirichlet(:, 2), interp_readout_positions_150_cells_dirichlet(:, 3), 'VariableNames', {'mu_c0', 'std', 'SE'}), filename_dirichlet_interp_150); 

end 

%% functions for the ODE

% function to solve the ODE with Dirichlet BC 
function sol_dirichlet = solve_ode_dirichlet(x0, nc, CV, CV_bc, n, tol, mu_d, mu_D, c0)

    options = bvpset('Vectorized', 'on', 'NMax', 100*nc, 'RelTol', tol, 'AbsTol', tol);

    % default: all parameters constant
    d = mu_d * ones(nc, 1);
    D = mu_D * ones(nc, 1);

    % draw random kinetic parameters for each cell
    d = random(helper_functions.logndist(mu_d, CV), nc, 1); 
    D = random(helper_functions.logndist(mu_D, CV), nc, 1);

    dir_c0 = random(helper_functions.logndist(c0, CV_bc), 1, 1);
     
    % get initial solution 
    sol0 = bvpinit(x0, @helper_functions.y0_non_lin);

    odefun_init = @(x,y,c) helper_functions.odefun_non_lin_no_source(x, y, c, n, D, d);    
    odefun_bc = @(ya, yb) helper_functions.bcfun_dirichlet(ya, yb, nc, dir_c0);

    % solve the equation
    sol_dirichlet = bvp4c(odefun_init, odefun_bc, sol0, options);    
end 

end 
