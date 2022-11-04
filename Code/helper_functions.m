% helper functions to run the simulations for non-linear decay 
classdef helper_functions

    methods (Static)

        %% functions for the steady-state reaction-diffusion ODE
        % D*C''(x) - p*H(x) + d*C(x) = 0
        % reaction-diffusion equation
        function dydx = odefun(x, y, c, D, p, d, ncS)
            dC = -y(2,:) / D(c); % mass flux: j = -D*grad(C)
            dj = p(c) * (c <= ncS) - d(c) * y(1,:); % conservation of mass: div(j) = p*H(-x) - d*C
            dydx = [dC; dj];
        end

        % initial guess
        function y = y0(x, c)
            y = [0; 0];
        end

        % boundary & cell interface conditions (zero flux at both ends) 
        function res = bcfun(ya, yb, nc)
            res = ya(:);
            res(1) = ya(2, 1); % zero flux at the left end of the source domain
            res(2) = yb(2,nc); % zero flux at right end of the patterning domain

            for c = 1:nc-1
                res(2*c+1) = ya(1,c+1) - yb(1,c); % concentration continuity
                res(2*c+2) = ya(2,c+1) - yb(2,c); % flux continuity
            end
        end


        %% reaction diffusion equation with non-linear decay and source 
        % D*C''(x) - p*H(x) + d*C(x)^n/C_ref^(n-1) = 0
        function dydx = odefun_non_lin(x, y, c, n, D, p, d, ncS, C_ref)
            
            dC = -y(2,:) / D(c); % mass flux: j = -D*grad(C) 
            dj = p(c) * (c <= ncS) - d(c) * y(1,:).^n / C_ref^(n-1); % conservation of mass: div(j) = p*H(-x) - d*C^n/C_ref^(n-1)
            dydx = [dC; dj];
            
        end

        % Dirichlet BC at x=0 and zero flux at x=LP
        function res = bcfun_dirichlet(ya, yb, nc, c0)
            
            res = ya(:);
            res(1) = ya(1, 1) - c0; 
            res(2) = yb(2,nc);
            
            for c = 1:nc-1
                res(2*c+1) = ya(1,c+1) - yb(1,c); % concentration continuity
                res(2*c+2) = ya(2,c+1) - yb(2,c); % flux continuity
            end
            
        end
        
        % Flux BC at x=0 and zero flux at x=LP
        function res = bcfun_flux(ya, yb, nc, j)
            
            res = ya(:);
            res(1) = ya(2, 1) - j; 
            res(2) = yb(2,nc);
            
            for c = 1:nc-1
                res(2*c+1) = ya(1,c+1) - yb(1,c); % concentration continuity
                res(2*c+2) = ya(2,c+1) - yb(2,c); % flux continuity
            end
            
        end

        % initial guess for non-linear decay 
        function y = y0_non_lin(x, c)

            y = [0.5; 0];

        end
        
        %% ODE without source term 
        % D*C''(x) + d*C(x)^n/C_ref^(n-1) = 0
        function dydx = odefun_non_lin_no_source(x, y, c, n, D, d, C_ref)
            
            dC = -y(2,:) / D(c); % mass flux: j = -D*grad(C) 
            dj = - d(c) * y(1,:).^n / C_ref^(n-1); % conservation of mass: div(j) = p*H(-x) - d*C^n/C_ref^(n-1)
            dydx = [dC; dj];
            
        end

        %% log-normal distribution with prescribed mean and coefficient of variation

        function pd = logndist(mu, CV)
            pd = makedist('Lognormal', 'mu', log(mu/sqrt(1+CV^2)), 'sigma', sqrt(log(1+CV^2)));
        end


        %% function to set up the spatial domain

        % builds a one-dimensional domain consisting of variably-sized cells
        function [l_s, l_p] = build_domain(LS, LP, mu_d, CV_A)

            % Inupt parameters:
            % LS =      length of source domain
            % LP =      length of patterning dmoain
            % mu_d =    mean cell diameter
            % CV_A =    coefficient of variation for cell areas

            % mean cell area
            mu_A = pi * (mu_d/2)^2 * (CV_A^2+1)^(1/4);

            % Build the source domain
            l_s = [];
            while sum(l_s) < LS
                A = random(helper_functions.logndist(mu_A, CV_A), 1, 1);
                diam = 2*sqrt(A/pi);
                l_s = [l_s, diam];
            end

            % remove the last cell if the domain length is closer to the target without it
            if length(l_s) > 1 && abs(sum(l_s) - LS) >= abs(sum(l_s) - l_s(end) - LS)
                l_s = l_s(1:end-1);
            end

            % calculate the normalised diameter for each cell
            l_s = fliplr(cumsum(l_s / sum(l_s) * LS));

            % Build the patterning domain
            l_p = [];
            while sum(l_p) < LP
                A = random(helper_functions.logndist(mu_A, CV_A), 1, 1);
                diam = 2*sqrt(A/pi);
                l_p = [l_p, diam];
            end

            % remove the last cell if the domain length is closer to the target without it
            if length(l_p) > 1 && abs(sum(l_p) - LP) >= abs(sum(l_p) - l_p(end) - LP)
                l_p = l_p(1:end-1);
            end

            % calculate the normalised diameter for each cell
            l_p = cumsum(l_p / sum(l_p) * LP);

        end
        
        %% analytical solution to the equation:
        % D*C''(x) + d*C(x)^n = 0
        function c = get_readout_conc_non_linear(x, n, c0, mu_lambda, C_ref)
    
            % get the concentration profile for non-linear decay. 
    
            % C(0)= c_0
            m = 2/(n - 1);
            lambda_m = mu_lambda*sqrt(1 + 1/m)*(C_ref/c0)^(1/m);
            c = c0*(1 + x/(m*lambda_m)).^(-m);  
        
        end

        %% function that returns the index where threshold conc. is reached 
        function pos = getindex(array, K, domain)
    
            % allocate memory 
            index = NaN(1, length(K));
            pos = NaN(1, length(K));
        
            % loop over concentrations and retrieve x position 
            for conc = 1:length(K)   
                if ~isnan(find(array <= K(conc)))
                    index(1, conc) = find((array <= K(conc)), 1);
                    pos(1, conc) = domain(index(1, conc));            
                end
            end
            
        end
        
        %% get average concentration of the morphogen in each cell along the patterning domain 
        function y_sol_average = average_concentration(sol_x, sol_y, patterning_domain, y_sol_average, nc)
    
            % initialise the start location at the beginning of ther
            % patterning domain 
            cell_beginning = patterning_domain(1);
        
            % loop through the  cell and extract the solutions for each
            % cell separately. Then use the trapezoid method to
            % numerically integrate and find the mean morphogen
            % gradient concentration over one cell.
        
            for cell_loc = 1:nc
                
                % set the upper interval as the end of a cell 
                cell_end = patterning_domain(cell_loc+1);
                
                % define interval where to extract solutions 
                logical_indexes = (sol_x <= cell_end) & (sol_x >= cell_beginning);
        
                % extract indices of the desired solutions 
                interval = find(logical_indexes);
        
                % get lenght of the cell for normalisation
                cell_length = cell_end - cell_beginning;
        
                % get the x and y solution 
                X = sol_x(interval);            
                Y = sol_y(1, interval);
        
                % get unique x, y values for interpolation solver
                x_unique = unique(X,'stable');
                y_unique = unique(Y,'stable');
        
                % increase the resoultion of points in each cell 
                x_high_res = linspace(cell_beginning, cell_end, 100);
        
                % get interpolated solutions for better resolution
                y_high_res = pchip(x_unique, y_unique, x_high_res);
        
                % get the average concentration per cell 
                trapz_sol = trapz(x_high_res,y_high_res)/cell_length;
        
                % append the solution for each cell to the solution
                % array 
                y_sol_average(cell_loc, 1) = trapz_sol;
        
                % set the lower interval for the next iteration as the
                % current end of the cell 
                cell_beginning = cell_end;
        
        
            end
        
        end

    end
end
