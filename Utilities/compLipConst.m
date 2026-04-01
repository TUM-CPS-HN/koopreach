function [gamma, L] = compLipConst(fun,U,R0,steps,initpoints,dim_x)
% compLipConst3 - Computes an approximation of the Lipschitz constant L and
% the covering radius gamma (delta) for a given dataset.
%
% This function correctly computes a per-dimension L using the user-specified
% calculation method and a single, global gamma value for the full input space.
%
% Syntax:
%     [gamma, L] = compLipConst3(fun,U,R0,steps,initpoints,dim_x)
%
% Inputs:
%    fun          - A handle to the function to analyze (e.g., @cstrDiscr2).
%    U            - The set or domain for the input vector 'u'.
%    R0           - The initial state set for the trajectories.
%    steps        - The number of steps for each trajectory.
%    initpoints   - The number of initial points for trajectories.
%    dim_x        - The dimension of the state vector 'x'.
%
% Outputs:
%    gamma      - The approximated covering radius for the entire input space.
%    L          - A 1x(dim_x) vector containing the approximated Lipschitz
%                 constants for each output dimension.

    rand('seed',1);
    totalsamples = initpoints*steps;

    % input random sample points
    dim_u = size(randPointExtreme(U), 1);
    u = zeros(dim_u, totalsamples);
    for i=1:totalsamples
        u(:,i) = randPointExtreme(U);
    end

    % get state trajectories
    x_free = zeros(dim_x*initpoints, steps+1);
    index = 1;
    for j=1:dim_x:initpoints*dim_x
        x_free(j:j+dim_x-1,1) = randPoint(R0*0.1);
        for i=1:steps
            x_free(j:j+dim_x-1,i+1) = fun(x_free(j:j+dim_x-1,i),u(:,index));
            index=index+1;
        end
    end

    % combine trajectories
    x_free_vec_0 = zeros(dim_x, totalsamples);
    x_free_vec_1 = zeros(dim_x, totalsamples);
    index_0 = 1;
    index_1 = 1;
    for j=1:dim_x:initpoints*dim_x
        for i=2:steps+1        
            x_free_vec_1(:,index_1) = x_free(j:j+dim_x-1,i);
            index_1 = index_1 + 1;
        end
        for i=1:steps
            x_free_vec_0(:,index_0) = x_free(j:j+dim_x-1,i);
            index_0 = index_0 + 1;
        end
    end
    
    % --- Step 1: Compute the approximated Lipschitz constant (L) ---
    % This section implements the user's specific logic for L.
    normType = 2; % Using the L2-norm as specified
    L = zeros(1, dim_x);

    for idx = 1:dim_x
        for i = 1:totalsamples
            % z vector is constructed with only a single state dimension
            z1 = [x_free_vec_0(:, i); u(:, i)];
            f1 = x_free_vec_1(idx, i);
            for j = 1:totalsamples
                z2 = [x_free_vec_0(:, j); u(:, j)];
                f2 = x_free_vec_1(idx, j);
                if norm(z1 - z2, 2) > 1e-10
                    newnorm = norm(f1 - f2, normType) ./ norm(z1 - z2, normType);
                    if newnorm > L(idx)
                        L(idx) = newnorm;
                    end
                end
            end
        end
    end
    
    % --- Step 2: Compute the approximated covering radius (gamma) ---
    % This section computes gamma using the standard, global definition.
    % Combine state and input to form the 'z' vectors
    z_data = [x_free_vec_0; u];
    min_distances = inf(1, totalsamples); % Store the min distance for each point
    
    % Use nested loops to find the nearest neighbor for each point.
    for i = 1:totalsamples
        for j = 1:totalsamples
            % Avoid self-comparison
            if i ~= j
                z_i = z_data(:, i);
                z_j = z_data(:, j);
                
                % Calculate the distance in the full input space
                dist = norm(z_i - z_j, 2);
                
                % Update the minimum distance for point i
                if dist < min_distances(i)
                    min_distances(i) = dist;
                end
            end
        end
    end
    
    % The covering radius is the maximum of the minimum distances
    gamma = max(min_distances);
end
