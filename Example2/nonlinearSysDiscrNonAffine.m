function f = nonlinearSysDiscrNonAffine(x, u, T)
% nonlinearSysDiscr - discrete-time nonlinear system using Euler method
%
%   x = [x1; x2]        state vector
%   u = [u1; u2]        input vector
%   T = sampling time dt
%
%   f = next state vector
%
%   Implements:
%       x1_next = x1 + T * (mu*x1 - x1 + x1*exp(u1))
%       x2_next = x2 + T * (lam*(x2 - x1^2) - x2 + (u1*u2 + x2*exp(u2)))
%
% Author: Converted from Python AutoKoopman example

%------------- BEGIN CODE --------------

% System parameters
mu  = -0.6;
lam = -1.0;

% State
x1 = x(1);
x2 = x(2);

% Inputs (default to zero if not provided)
if nargin < 2 || isempty(u)
    u1 = 0;
    u2 = 0;
else
    u1 = u(1);
    u2 = u(2);
end

% Continuous-time derivatives
dx1 = mu*x1 - x1 + x1*exp(u1);
dx2 = lam*(x2 - x1^2) - x2 + (u1*u2 + x2*exp(u2));

% Euler discretization
x1_next = x1 + T * dx1;
x2_next = x2 + T * dx2;

% Output
f = [x1_next; x2_next];

%------------- END OF CODE --------------

end
