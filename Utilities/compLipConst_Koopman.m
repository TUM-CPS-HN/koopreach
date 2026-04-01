function [gamma, L] = compLipConst_Koopman(fun,U,R0,steps,initpoints,dim_x)

    rand('seed',1);
    totalsamples = initpoints*steps;

    % input random sample points
    dim_u = size(randPointExtreme(U), 1);
    u = zeros(dim_u, totalsamples);
    for i = 1:totalsamples
        u(:,i) = randPointExtreme(U);
    end

    % get state trajectories
    x_free = zeros(dim_x*initpoints, steps+1);
    index = 1;
    for j = 1:dim_x:initpoints*dim_x
        x_free(j:j+dim_x-1,1) = randPoint(R0);
        for i = 1:steps
            x_free(j:j+dim_x-1,i+1) = fun(x_free(j:j+dim_x-1,i), u(:,index));
            index = index + 1;
        end
    end

    % combine trajectories
    x_free_vec_0 = zeros(dim_x, totalsamples);
    x_free_vec_1 = zeros(dim_x, totalsamples);
    index_0 = 1;
    index_1 = 1;
    for j = 1:dim_x:initpoints*dim_x
        for i = 2:steps+1
            x_free_vec_1(:,index_1) = x_free(j:j+dim_x-1,i);
            index_1 = index_1 + 1;
        end
        for i = 1:steps
            x_free_vec_0(:,index_0) = x_free(j:j+dim_x-1,i);
            index_0 = index_0 + 1;
        end
    end

    % --- Step 1: Lipschitz constant (state only, input frozen) ---
    normType = 2;
    L = zeros(1, dim_x);

    for idx = 1:dim_x
        for i = 1:totalsamples
            z1 = x_free_vec_0(:,i);
            u_fixed = u(:,i);                    % ✅ FIX: freeze input
            f1 = fun(z1, u_fixed);

            for j = 1:totalsamples
                z2 = x_free_vec_0(:,j);
                f2 = fun(z2, u_fixed);           % ✅ same input

                dz = norm(z1 - z2, normType);
                if dz > 1e-10
                    newnorm = abs(f1(idx) - f2(idx)) / dz;
                    if newnorm > L(idx)
                        L(idx) = newnorm;
                    end
                end
            end
        end
    end

    % --- Step 2: Covering radius gamma (state space only) ---
    z_data = x_free_vec_0;
    min_distances = inf(1, totalsamples);

    for i = 1:totalsamples
        for j = 1:totalsamples
            if i ~= j
                dist = norm(z_data(:,i) - z_data(:,j), 2);
                if dist < min_distances(i)
                    min_distances(i) = dist;
                end
            end
        end
    end

    gamma = max(min_distances);
end
