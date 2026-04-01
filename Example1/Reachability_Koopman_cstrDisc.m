% t_nonlinearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = f(x(k),u(k)) + w(k)
%------------- BEGIN CODE --------------
clear
close all
warning off
rand('seed',1);
% Simulation & Timing
params.dt = 0.01;
NN = 30;
params.tFinal = params.dt*NN;
params.tStart = 0;
% Reachability Sets
params.U = zonotope([[-0.4; 1.4], diag([.1, .2])]); % Input set
params.R0 = zonotope([[0.5; 0.4], diag([0.05; 0.3])]); % Initial set
RR0 = params.R0;
% Dimensions and Sampling
options.dim_x = 2;
initpoints = 70;
steps = 100;
totalsamples = steps * initpoints;
% Disturbance / Noise
sigma_v = 0.0001;
options.W = zonotope(zeros(options.dim_x,1), sigma_v*diag(ones(options.dim_x,1)));
% 1. System Dynamics (CSTR)
fun = @(x,u) cstrDiscr(x, u, params.dt);
% Reachability Settings  --------------------------------------------------
options.zonotopeOrder = 100;
options.tensorOrder = 3;
options.errorOrder = 5;
% System Dynamics  --------------------------------------------------------
sysDisc = nonlinearSysDT('stirredTankReactor',fun,params.dt);
% Reachability Analysis ---------------------------------------------------
R_model = reach(sysDisc,params,options);

% Data Generation
u = zeros(2, totalsamples);
for i=1:totalsamples
    u(:,i) = randPointExtreme(params.U);
end

%get state trajectories
index=1;
for j=1:initpoints
    x_traj_temp(:,1) = randPoint(params.R0);
    for i=1:steps
        x_traj_temp(:,i+1) = fun(x_traj_temp(:,i),u(:,index)) + randPoint(options.W);
        index=index+1;
    end
    X0_traj{j,1} = x_traj_temp(:,1:steps);
    X1_traj{j,1} = x_traj_temp(:,2:steps+1);
end

% Combine for ID
options.X_0T = cell2mat(X0_traj');
options.X_1T = cell2mat(X1_traj');
options.U_full = u;

sigma_v = 0.0001;
options.W = zonotope(zeros(options.dim_x,1),sigma_v*diag(ones(options.dim_x,1))); % disturbance

dim_x = options.dim_x;
GW = cell(dim_x * totalsamples, 1);
index = 1;
for row = 1:dim_x
    for col = 1:totalsamples
        G = zeros(dim_x, totalsamples);
        G(row, col) = sigma_v;
        GW{index} = G;
        index = index + 1;
    end
end
options.Wmatzono = matZonotope(zeros(dim_x, totalsamples), GW);

% data-driven(LS) based Reachability Analysis
stepsLip=1;
initpointsLip=1000;
[gamma,L]= compLipConst(fun,params.U,params.R0,stepsLip,initpointsLip,options.dim_x);
eps(1)= L(1) .* gamma/2;
eps(2)= L(2) .* gamma/2;
options.Zeps = zonotope([zeros(2,1),diag(eps)]);

 R_data = reach_LS(params,options);

%%  Koopman Identification
% State Lifting: [1; x1; x2; x1^2; x1*x2; x2^2] (6x1)
obs_edmd = @(x) [ones(1, size(x,2)); x(1,:); x(2,:); x(1,:).^2; x(1,:).*x(2,:); x(2,:).^2];
options.dim_z = 6;

% Symbolic for Jacobian
syms x1 x2 real
phi_x_sym = [1; x1; x2; x1^2; x1*x2; x2^2];
Jacobian_phi_sym = jacobian(phi_x_sym, [x1; x2]);

% Interaction Lifting V(x,u): (12x1)
% V_xu_edmd = @(x,u) [u(1,:); u(2,:); ...
%     x(1,:).*u(1,:); x(1,:).*u(2,:); ...
%     x(2,:).*u(1,:); x(2,:).*u(2,:); ...
%     u(1,:).^2; u(1,:).*u(2,:); u(2,:).^2; ...
%     (x(1,:).^2).*u(1,:); (x(2,:).^2).*u(2,:); ...
%     ones(1, size(x,2))];
V_xu_edmd = @(x,u) [u(1,:); u(2,:)];

Theta_curr = obs_edmd(options.X_0T);
Theta_next = obs_edmd(options.X_1T);
Upsilon    = V_xu_edmd(options.X_0T, options.U_full); % FIX: Pass physical X, not Theta
K = Theta_next * pinv([Theta_curr; Upsilon]);
A = K(:, 1:options.dim_z);
B = K(:, options.dim_z+1:end);
C_proj = [0 1 0 0 0 0; 0 0 1 0 0 0]; 

%% ================= 3. CONSERVATIVE RESIDUAL BOUNDING =============
for k = 1:initpoints
    X0_ = X0_traj{k,1};
    X1_ = X1_traj{k,1};
    U_traj = u(:,(k-1)*steps+1:k*steps);
    
    for j = 1:steps-NN
        Z_initial = obs_edmd(X0_(:,j));
        Z_sim = Z_initial;
        Z_actual = obs_edmd(X1_(:,j:j+NN-1));
        
        for m = 1:NN
            u_curr = U_traj(:,j+m-1);
            x_phys = Z_sim(2:3); 
            
            [A_lin, B_lin] = linearizeBilinear_updated(A, B, Z_sim, u_curr);
            V_nominal = V_xu_edmd(x_phys, u_curr);
            
            Z_sim = (A*Z_sim + B*V_nominal) ;
            
            Residuals_x{m,1}(:, (k-1)*(steps-NN) + j) = C_proj * (Z_actual(:,m) - Z_sim);
        end
    end
end

for i = 1:NN
    VInt = intervalMatrix(Residuals_x{i,1});
    Error_Zono{i} = zonotope(interval(min(infimum(VInt)')', max(supremum(VInt)')'));
end

%% ================= 5. GUARANTEED REACHABILITY LOOP ==============
x_star = RR0.center;
Jacobian_phi = double(subs(Jacobian_phi_sym, [x1 x2], x_star'));
G = Jacobian_phi * RR0.generators;
R = zonotope(obs_edmd(x_star), G);
z_nom = R.center;

params.R0 = R;
fun_Koopman = @(z,u) A * z + B*V_xu_edmd(z,u);
stepsLip=1;
initpointsLip=1000;
[gamma,L_psi]= compLipConst_Koopman(fun_Koopman,params.U,params.R0,stepsLip,initpointsLip,options.dim_z);
G_W = diag(L_psi) *sigma_v ;
options.W1 = zonotope(zeros(options.dim_z,1),G_W); % disturbance
options.U_full = ones(2,length(options.U_full)).*params.U.center;
sysDisc_koopman = nonlinearSysDT('nonlinearSysDiscrNonAffineBilinear',fun_Koopman,params.dt,options.dim_z,2);

R_Koopman0= reach_Koopman(sysDisc_koopman,params,options);

for i = 1:NN
  
     R_Koopman1{i} = C_proj * R_Koopman0.timePoint.set{i+1} + Error_Zono{i}+options.Zeps;
    t_vec(k) = k;

end
% Data driven reach set object
time_idx = 1:NN;
timePoint.set = R_Koopman1(time_idx);
timePoint.time = num2cell(t_vec(time_idx));
R_Koopman = reachSet(timePoint);

%% Visualization
figure('Renderer', 'painters', 'Position', [10 10 700 600]); hold on; box on;
handleModel=plot(R_model, [1 2], 'b', 'FaceAlpha', 0.2);
handleDataKoopman=plot(R_Koopman, [1 2], 'g', 'FaceColor', 'none', 'LineWidth', 1.5);
handleData =  plot(R_data,[1 2], ...
    'EdgeColor','r', ...
    'LineWidth',1, ...
    'FaceColor','none');   % outline only, no fill

% handleKoopman2 =  plot(R_Koopman2,[1 2], ...
%     'EdgeColor','y', ...
%     'LineWidth',1, ...
%     'FaceColor','none');   % outline only, no fill

% Simulation --------------------------------------------------------------
% --- 1. Setup CORA Simulation Options ---
simOpt.points = 100;
simOpt.type = 'standard';
simOpt.fracVert = 0;      % 0 = interior of X0, 1 = vertices of X0
simOpt.fracInpVert = 0;   % 1 = matches your randPointExtreme(params.U)
simOpt.nrConstInp = NN;   % Changes input at every step (like your loop)
params.R0=RR0;
% Ensure Noise is zero in the system
% (Assuming sysDisc was created with options.W = zonotope(zeros(dim_x,1)))
simRes = simulateRandom(sysDisc, params, simOpt);
% plot simulation
handleTraj=plot(simRes,[1,2],'Marker','.','LineStyle','none','color','k');
handleX0 = plot(RR0, [1 2], 'c', 'LineWidth', 1);

% --- 2. Setup Manual Loop to match ---
% for j=1:200
%     xt = randPoint(RR0); % Match simOpt.fracVert = 0
%     path = [];
%     for i=1:NN
%         % Match simOpt.fracInpVert = 1
%         u_random = randPointExtreme(params.U); 
% 
%         % No noise (+ 0)
%         xt = fun(xt, u_random); 
%         path = [path xt];
%     end
%     plot(path(1,:), path(2,:), '.k');
% end
xlabel('$x_1$ ','Interpreter', 'latex');
ylabel('$x_2$ ','Interpreter', 'latex');
legend([handleX0,handleTraj,handleModel,handleData,handleDataKoopman],...
    'Initial set ($\mathcal{X}_0$)','Random Trajectories','Reachable Sets ($\mathcal{R}_k$) from Model (CORA)','Reachable Sets ($\tilde{\mathcal{R}}_k$) from data using LS','Reachable Sets ($\mathcal{R}^\prime_k$) from data using Koopman','Location','northwest','Interpreter','latex');
legend boxoff

ax = gca;
ax.FontSize = 20;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3) - 0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left-0.01 bottom+0.03 ax_width+0.01 ax_height-0.03];
ylim([-10 53]); % Example: from -2 to 5 in x_2
 xlim([-0.4 0.6]); % Example: from -2 to 5 in x_2


%% ------------- END OF CODE --------------

function f = cstrDiscr(x,u,T)
    rho = 1000; Cp = 0.239; deltaH = -5e4; E_R = 8750; k0 = 7.2e10;
    UA = 5e4; q = 100; Tf = 350; V = 100; C_Af = 1;
    C_A0 = 0.5; T_0 = 350; T_c0 = 300;
    U_ctrl = [-3 -6.9] * x + T_c0;
    xs = x + [C_A0; T_0]; 
    f(1,1) = ((1-(q*T)/(2*V) - k0*T*exp(-E_R/xs(2)))*xs(1) + q/V * C_Af * T) / (1 + (q*T)/(2*V)) + u(1)*T;
    f(2,1) = (xs(2)*(1-0.5*T*q/V - (T*UA)/(2*V*rho*Cp)) + T*(Tf*q/V + (UA*U_ctrl)/(V*rho*Cp)) ...
              - xs(1)*(deltaH*k0*T)/(rho*Cp) * exp(-E_R/xs(2))) / (1+0.5*T*q/V+(T*UA)/(2*V*rho*Cp)) + u(2)*T;
    f = f - [C_A0; T_0];
end
% 
% function [A_lin, B_lin] = linearizeBilinear_updated(A, B, z0, u0)
%     u1 = u0(1); u2 = u0(2);
%     x1 = z0(2); x2 = z0(3);
%     dVdz = zeros(12, 6);
%     dVdz(3, 2) = u1; dVdz(4, 2) = u2; 
%     dVdz(5, 3) = u1; dVdz(6, 3) = u2; 
%     dVdz(10, 2) = 2*x1*u1; dVdz(11, 3) = 2*x2*u2; 
%     dVdu = zeros(12, 2);
%     dVdu(1, 1) = 1; dVdu(2, 2) = 1;
%     dVdu(3, 1) = x1; dVdu(4, 2) = x1;
%     dVdu(5, 1) = x2; dVdu(6, 2) = x2;
%     dVdu(7, 1) = 2*u1; dVdu(8, 1) = u2; dVdu(8, 2) = u1; dVdu(9, 2) = 2*u2;
%     dVdu(10, 1) = x1^2; dVdu(11, 2) = x2^2;
%     A_lin = A + B * dVdz;
%     B_lin = B * dVdu;
% end


function [A_lin, B_lin] = linearizeBilinear_updated(A, B, z0, u0)
    %#ok<INUSD> % z0 and u0 are unused

    A_lin = A;
    B_lin = B;
end