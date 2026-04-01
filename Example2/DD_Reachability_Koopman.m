% t_nonlinearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = f(x(k),u(k)) + w(k)
% The approach is based on [1].
% Author:       Alireza Naderi Akhormeh
% Written:      20-Feb-2026
% Last update:
% Last revision:---

%------------- BEGIN CODE --------------

clear
close all
warning off
rand('seed',1);
params.dt = 0.01;
NN = 30;
params.tFinal = params.dt*NN;
params.tStart = 0;
%input set
params.U = zonotope([[-0.4;1.4],diag([.1 .2])]);
%initial set
params.R0 = zonotope([[0.5;0.4],diag([0.05;0.3])]);
RR0 = params.R0;
% dimension of x
options.dim_x=2;

%Number of trajectories
initpoints=30;
%Number of time steps
steps=100;

%Totoal number of samples
totalsamples = steps*initpoints;

%noise zonotope
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

% Reachability Settings
options.zonotopeOrder = 100;
options.tensorOrder = 2;
options.errorOrder = 15;

% System Dynamics
fun = @(x,u) nonlinearSysDiscrNonAffine(x,u,params.dt);

%input random sample points
for i=1:totalsamples
    u(:,i) = randPointExtreme(params.U);

end

%get state trajectories
index=1;
for j=1:options.dim_x:initpoints*options.dim_x

    x(j:j+options.dim_x-1,1) = randPoint(params.R0);

    x_free(j:j+options.dim_x-1,1) = x(j:j+options.dim_x-1,1);
    for i=1:steps
        x_free(j:j+options.dim_x-1,i+1) = fun(x(j:j+options.dim_x-1,i),u(:,index));
        x(j:j+options.dim_x-1,i+1) = fun(x(j:j+options.dim_x-1,i),u(:,index)) +randPoint(options.W);
        index=index+1;
    end
end


%combine trajectories
index_0 =1;
index_1 =1;
for j=1:options.dim_x:initpoints*options.dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x(j:j+options.dim_x-1,i);
        x_free_vec_1(:,index_1) = x_free(j:j+options.dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_free_vec_0(:,index_0) = x_free(j:j+options.dim_x-1,i);
        x_meas_vec_0(:,index_0) = x(j:j+options.dim_x-1,i);
        index_0 = index_0 +1;
    end
end


index_0 =1;
index_1 =1;
k = 1;
for j=1:options.dim_x:initpoints*options.dim_x
    for i=2:steps+1
        x_vec_1(:,index_1) = x(j:j+options.dim_x-1,i);
        index_1 = index_1 +1;
    end
    X1_traj{k,1} = x_vec_1;
    index_1 = 1;

    for i=1:steps
        x_vec_0(:,index_0) = x(j:j+options.dim_x-1,i);
        index_0 = index_0 +1;
    end
    X0_traj{k,1} = x_vec_0;
    index_0 = 1;
    k = k+1;
end


stepsLip=1;
initpointsLip=1000;
[gamma,L]= compLipConst(fun,params.U,params.R0,stepsLip,initpointsLip,options.dim_x);
eps(1)= L(1) .* gamma/2;
eps(2)= L(2) .* gamma/2;
options.Zeps = zonotope([zeros(2,1),diag(eps)]);
Zeps=options.Zeps;

% X_+ is X_1T
% X_- is X_0T
options.U_full = u(:,1:totalsamples);
options.X_0T = x_meas_vec_0(:,1:totalsamples);
options.X_1T = x_meas_vec_1(:,1:totalsamples);

% define system
sysDisc = nonlinearSysDT('nonlinearSysDiscrNonAffine',fun,params.dt,2,2);

%% Reachability Analysis ---------------------------------------------------
% compute model based reachability (R) and data driven one (R_data)
tic
% model based Reachability Analysis
R_model= reach(sysDisc,params,options);
% data-driven(LS) based Reachability Analysis
R_data = reach_LS(params,options);
tComp = toc;
disp("Computation time: " + tComp);

%%  Koopman
%initial set
options.obs = @(x) [ ...
    x(1);          % x0
    x(2);          % x1
    x(1)^2;        % x0^2
    x(1)*x(2);     % x0*x1
    x(2)^2;        % x1^2
    ];

syms x1 x2 real
phi_x_sym =[ ...
    x1;
    x2;
    x1^2;
    x1*x2;
    x2^2;
    ];

x_sym = [x1; x2];

Jacobian_phi_sym = jacobian(phi_x_sym, x_sym);

% dimension of x
options.dim_z=5;
% Full lifting (state + input)
V_x_u = @(x,u) [ ...
    u(1);                   % 0: u0
    u(2);                   % 11: u1
    x(1)*u(1);              % 12: x0*u0
    x(1)*u(2);              % 13: x0*u1
    x(2)*u(1);              % 14: x1*u0
    x(2)*u(2);              % 15: x1*u1
    u(1)^2;                 % 16: u0^2
    u(1)*u(2);              % 17: u0*u1
    u(2)^2;                 % 18: u1^2
    ];
options.V_x_u = V_x_u;
obs_edmd = @(x) [ ...
    x(1,:);          % x0
    x(2,:);          % x1
    x(1,:).^2;        % x0^2
    x(1,:).*x(2,:);     % x0*x1
    x(2,:).^2;        % x1^2
    ];

V_xu_edmd = @(x,u) [ ...
    u(1,:);                   % 0: u0
    u(2,:);                   % 11: u1
    x(1,:).*u(1,:);              % 12: x0*u0
    x(1,:).*u(2,:);              % 13: x0*u1
    x(2,:).*u(1,:);              % 14: x1*u0
    x(2,:).*u(2,:);              % 15: x1*u1
    u(1,:).^2;                 % 16: u0^2
    u(1,:).*u(2,:);              % 17: u0*u1
    u(2,:).^2;                 % 18: u1^2
    ];

Theta_curr = obs_edmd(options.X_0T);
Theta_next = obs_edmd(options.X_1T);
Upsilon    = V_xu_edmd(Theta_curr, options.U_full);

K = Theta_next * pinv([Theta_curr; Upsilon]);
A = K(:, 1:5); B = K(:, 6:end);
% [A_bilinear, B_bilinear] = Koopman_Bilinear_model();

C_proj = [1 0 0 0 0 ; 0 1 0 0 0 ];
%% ================= 3. CONSERVATIVE RESIDUAL BOUNDING =============
test_start_idx = 1;
X0 = options.X_0T(:, test_start_idx);
k=1;

for k = 1:initpoints
    X0_ =  X0_traj{k,1};
    X1_ = X1_traj{k,1};
    U_  = u(:,1:steps);

    for j = 1:length(X0_)-NN-2


        Theta_curr = obs_edmd(X0_);
        Theta_next = obs_edmd(X1_);
        Upsilon    = V_xu_edmd(Theta_curr, U_);

        Z0 = Theta_curr(:,j:j+NN) ;
        Z1 = Theta_next(:,j:j+NN) ;
        U_0 = Upsilon(:,j:j+NN);

        Z_sim = Z0(:,1);
        for m = 1:NN

            [A_lin, B_lin] = linearizeBilinear_updated(A, B, Z_sim, U_0(1:2,m));
            Z_sim = A*Z_sim + B*V_xu_edmd(Z_sim, U_0(1:2,m));
            % Propagate lifted dynamics
            Residuals_z = C_proj * (Z1(:,m) - Z_sim);

            Residuals_x{m,1}(:,k) = Residuals_z;
        end
        k = k +1;
    end
end
for i = 1:NN
    VInt = intervalMatrix(Residuals_x{i,1});
    leftLimit = infimum(VInt);
    rightLimit = supremum(VInt);
    Error_Zono{i} = zonotope(interval(min(leftLimit')', max(rightLimit')'));
end


x_star = params.R0.center;
Jacobian_phi = double(subs(Jacobian_phi_sym, [x1 x2], x_star'));
G = Jacobian_phi*params.R0.generators;
params.R0 = zonotope(options.obs(params.R0.center),G);
fun_Koopman = @(z,u) A * z + B*V_x_u(z(1:2),u);

stepsLip=1;
initpointsLip=1000;
[gamma,L]= compLipConst_Koopman(fun_Koopman,params.U,params.R0,stepsLip,initpointsLip,options.dim_z);
eps=  L.* gamma/2;
options.Zeps_Koopman = zonotope([zeros(options.dim_z,1),diag(eps)])*0;
%

%noise zonotope
w_star = options.W.center;

Jacobian_phi = double(subs(Jacobian_phi_sym, [x1 x2], w_star'));

G_W = diag(L) *sigma_v ;
options.W1 = zonotope(options.obs(options.W.center),G_W); % disturbance


%% ================= 5. GUARANTEED REACHABILITY LOOP ==============
z_nom = params.R0.center;
R = params.R0;

params.R0 = R;
fun_Koopman = @(z,u) A * z + B*V_xu_edmd(z,u);

sysDisc_koopman = nonlinearSysDT('nonlinearSysDiscrNonAffineBilinear',fun_Koopman,params.dt,options.dim_z,2);

R_Koopman0= reach_Koopman(sysDisc_koopman,params,options);

for i = 1:NN
  
     R_Koopman1{i} = C_proj * R_Koopman0.timePoint.set{i+1} + Error_Zono{i}+options.Zeps;
    t_vec(k) = k;

end
% Data driven reach set object
time_idx = 2:3:NN;
timePoint.set = R_Koopman1(time_idx);
timePoint.time = num2cell(t_vec(time_idx));
R_Koopman = reachSet(timePoint);

%% Visualization -----------------------------------------------------------
figure('Renderer', 'painters', 'Position', [10 10 700 600]);hold on; box on;


% plot initial set
% plot model based reachable set
time = 0:params.dt:params.tFinal;
idx = 2:3:length(time);
timePoint.set = R_model.timePoint.set(idx);
timePoint.time = R_model.timePoint.time(idx);
R_sel = reachSet(timePoint);

handleModel=plot(R_sel,[1 2],'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b','LineWidth',1);

% % plot data driven reachable set
%
handleData =  plot(R_data,[1 2], ...
    'EdgeColor','r', ...
    'LineWidth',1, ...
    'FaceColor','none');   % outline only, no fill

handleDataKoopman=plot(R_Koopman, [1, 2],'EdgeColor','g', ...
    'LineWidth',1, ...
    'FaceColor','none');


%get state trajectories
index=1;
% NM = 100;
% for j=1:NM
%     x = randPoint(params.R0);
%     y = [];
%     for i=1:NN
%         x = fun(x,u(:,i)) +randPoint(options.W);
%         y = [y x];
%     end
%     hold on
%     handleTraj = plot(y(1,:),y(2,:),'.k','LineWidth',1);
% end
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

% formatting
xlabel('$x_1$','Interpreter', 'latex', 'FontSize', 30);
ylabel('$x_2$','Interpreter', 'latex', 'FontSize', 20);

% skip warning for extra legend entries
handleX0=plot(RR0,[1,2],'c-','LineWidth',1);

warOrig = warning; warning('off','all');
legend([handleX0,handleTraj,handleModel,handleData,handleDataKoopman],...
    'Initial set ($\mathcal{X}_0$)','Random Trajectories','Reachable Sets ($\mathcal{R}_k$) from model (CORA)','Reachable Sets ($\tilde{\mathcal{R}}_k$) from data using LS','Reachable Sets ($\mathcal{R}^\prime_k$) from data using Koopman','Location','northwest','Interpreter','latex');
warning(warOrig);
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
ylim([-1.5 4.5]); % Example: from -2 to 5 in x_2
xlim([0.25 0.6]); % Example: from -2 to 5 in x_2


%------------- END OF CODE --------------


function [A_lin, B_lin] = linearizeBilinear_updated(A, B, z0, u0)
% Inputs:
% z0: [x0; x1; x0^2; x0*x1; x1^2] (5x1)
% u0: [u0; u1] (2x1)

u_0 = u0(1);
u_1 = u0(2);
x_0 = z0(1);
x_1 = z0(2);

%% 1. Partial Derivatives of V w.r.t z (dV/dz) -> 9x5
dVdz = zeros(9, 5);
% Row 3: x0*u0 -> d/dx0 = u0
dVdz(3, 1) = u_0;
% Row 4: x0*u1 -> d/dx0 = u1
dVdz(4, 1) = u_1;
% Row 5: x1*u0 -> d/dx1 = u0
dVdz(5, 2) = u_0;
% Row 6: x1*u1 -> d/dx1 = u1
dVdz(6, 2) = u_1;

%% 2. Partial Derivatives of V w.r.t u (dV/du) -> 9x2
dVdu = zeros(9, 2);
% Row 1: u0 -> [1, 0]
dVdu(1, 1) = 1;
% Row 2: u1 -> [0, 1]
dVdu(2, 2) = 1;
% Row 3: x0*u0 -> [x0, 0]
dVdu(3, 1) = x_0;
% Row 4: x0*u1 -> [0, x0]
dVdu(4, 2) = x_0;
% Row 5: x1*u0 -> [x1, 0]
dVdu(5, 1) = x_1;
% Row 6: x1*u1 -> [0, x1]
dVdu(6, 2) = x_1;
% Row 7: u0^2 -> [2*u0, 0]
dVdu(7, 1) = 2 * u_0;
% Row 8: u0*u1 -> [u1, u0]
dVdu(8, 1) = u_1;
dVdu(8, 2) = u_0;
% Row 9: u1^2 -> [0, 2*u1]
dVdu(9, 2) = 2 * u_1;

%% 3. Combine into Linearized Matrices
A_lin = A + B * dVdz;
B_lin = B * dVdu;
end