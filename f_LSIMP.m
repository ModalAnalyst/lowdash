function [y,t,x,u,Info] = f_LSIMP(A,B,C,D,u,dt,x0,HO)
% Linear time-invariant system response to polynomial input
%
% [y,t,x] = f_LSIMP(A,B,C,D,u,dt,x0,HO)
%
% INPUT:
% A    [Nx,Nx] state matrix ([]: DefaultSystem)
% B    [Nx,Nu] input matrix ([]:  eye(Nx) )
% C    [Ny,Nu] output matrix ([]: eye(Nx) )
% D    [Ny,Nu] feedthrough matrix ([]: zeros(Ny,Nu) )
% u    [N,Nu] input vector ([]: randn(1000,1) );
% dt   [num] time step (1/fs) ([]: automatic selection)
% x0   [Nx,1] initial conditions ([]:zeros(Nx,1))
% HO   [0,1,3] hold order ([]:0, zero-order hold)
%
% OUTPUT:
% y    [size(u)-HO+1 Ny] outputs
% t    [size(u)-HO+1 1] time vector
% x    [size(u)-HO+1 Nx] states
%
% SYNTAX:
% [y,t,x] = f_LSIMP free decay of 1-DoF harmonic oscillator, wn = 1 rad/s, zeta = 0.01, x0 = [1 0];
% [y,t,x] = f_LSIMP(A) free decay with x0 = [1 0 ... 0];
% [y,t,x] = f_LSIMP(A,B) response to u=randn(DefalutNumberOfSamples,1) on input 1
% [y,t,x] = f_LSIMP(A,B,C) response to u=randn(DefalutNumberOfSamples,1) on input 1
% [y,t,x] = f_LSIMP(A,B,C,D) response to u=randn(DefalutNumberOfSamples,1) on input 1
% [y,t,x] = f_LSIMP(...,u) given input u
% [y,t,x] = f_LSIMP(...,dt) given time step (performs check if adequate)
% [y,t,x] = f_LSIMP(...,x0) given initial conditions
% [y,t,x] = f_LSIMP(...,HO) given hold order (0, 1, 3)
%
% SIMILAR FUNCTIONALITY: MATLAB's lsim.m (Control System Toolbox)
%
% EXAMPLE:
% Use ode45 to check first-order hold:
%
% % State-space model (https://de.mathworks.com/help/control/ug/mimo-state-space-models.html)
% A = [-0.0558 -0.9968  0.0802  0.0415
%       0.5980 -0.1150 -0.0318  0
%      -3.0500  0.3880 -0.4650  0
%       0       0.0805  1.0000  0];
% B = [ 0.0073       0
%      -0.4750  0.0077
%       0.1530  0.1430
%       0       0     ];
% C = [ 0       1       0       0;
%       0       0       0       1];
% D = zeros(2);
% 
% % Parameters
% N = 1001; % desired number of time steps
% dt = 1/20; % sample rate
% i = 1; j = 2; % desired response and excitation channels
% x0 = zeros(size(A,1),1); % initial conditions
% t0 = (0:N-1).'*dt; % time vector
% u = randn(N,numel(j)); % some excitation
% 
% % Compute reference solution
% ufcn = griddedInterpolant(t0,u,'linear'); % build interpolant
% [t0,x] = ode45(@ssfun,t0,x0,odeset('MaxStep',dt),A,B(:,j),ufcn); % integrate states
% y0 = x*C(i,:).' + u*D(i,j).'; % compute outputs
% 
% % Compute f_LSIMP solution for first-order hold
% [y2,t2] = f_LSIMP(A,B(:,j),C(i,:),D(i,j),u,dt,x0,1);
% 
% % Plot
% plot(t0,y0,'-',t2,y2,'-'), legend('ode45','f_LSIMP')
% 
% function dx = ssfun(t,x,A,B,ufcn)
% u = ufcn(t); % evaluate excitation at integrator time
% dx = A*x + B*u; % state equation
% end
%
% REFERENCE
% [1] G. Jelicic, M.Boeswald, A.Brandt 
%     "Improved computation in terms of accuracy and speed of LTI system
%     response with arbitrary input"
%     https://doi.org/10.1016/j.ymssp.2020.107252

arguments
  A(:,:) double {mustBeFinite} = [];
  B(:,:) double {mustBeFinite} = [];
  C(:,:) double {mustBeFinite} = [];
  D(:,:) double {mustBeFinite} = [];
  u(:,:) double {mustBeFinite} = [];
  dt double {mustBeScalarOrEmpty,mustBeNonnegative} = [];
  x0(:,1) double {mustBeFinite} = [];
  HO double {mustBeScalarOrEmpty,mustBeMember(HO,[0 1 3])} = [];
end

% Set defaults
DefaultStateMatrix = [0 1; -1 , -.04]; % harmonic oscillator, omega = 1 rad/s, zeta = 0.02;
DefaultNumberOfSamples = 100;
DefaultHoldOrder = 0;
DefaultTimeStep = 0; % let the function chose dt

% Check input
if isempty(A), A = DefaultStateMatrix; end
assert(size(A,1)==size(A,2),'f_LSIMP:StateMatrixSquare','The state matrix must be a square matrix')
if isempty(B), B = eye(size(A,1)); end
assert(size(A,1)==size(B,1),'f_LSIMP:InputMatrixSize','The input matrix must have as many rows as there are in the state matrix')
if isempty(C), C = eye(size(A)); end
assert(size(C,2)==size(A,1),'f_LSIMP:OutputMatrixSize','The output matrix must have as many columns as there are in the state matrix')
if isempty(D), D = zeros(size(C,1),size(B,2)); end
assert(size(D,1)==size(C,1) && size(D,2)==size(B,2),'f_LSIMP:FeedthroughMatrixSize','The feedthrough matrix must have as many rows as there are in the output matrix and as many columns as there are in the input matrix')
if isempty(u)
   u = zeros(DefaultNumberOfSamples,size(B,2));
   if any(B(:))
      u(:,1) = randn(size(u,1),1);
   elseif isempty(x0)
      x0 = zeros(size(A,1),1); 
      x0(1) = 1;
   end
end
assert(size(u,2)==size(B,2),'f_LSIMP:InputSizeMatch','The system inputs must have as many columns as there are columns of the input matrix')
if isempty(dt), dt = DefaultTimeStep; end
if isempty(x0), x0 = zeros(size(A,1),1); end
assert(size(A,1)==numel(x0),'f_LSIMP:InitialConditionsSizeMatch','The initial conditions must be as many as there are states')
if isempty(HO), HO = DefaultHoldOrder; end


%% Algorithm

% Check sample rate
[dt,HO,kf] = aux_checkSampleRate(A,HO,dt,any(u(:)));

% Compute excitation's polynomial coefficients
c = aux_getPolyCoeff(dt,u,HO);

% Discretize system matrices for desired hold order
[Ad,Bd] = aux_c2d(A,B,dt,HO); % [1] eq.(9)

% Integrate states
x = ltitr(Ad,Bd,c,x0); % [1] eq.(18) calculates x[n+1] = Ad*x[n]+Bd*c[n]

% Calculate outputs
y = x*C.'+u(1:size(x,1),:)*D.';

% Build time axis
t = (0:size(x,1)-1).'*dt;

% Collect metadata
Info = struct('TimeStep',dt,'HoldOrder',HO,'InitialConditions',x0,...
   'Algorithm',mfilename,'Duration',t(end),'OversamplingFactor',kf);


end

function [dt,HO,kf] = aux_checkSampleRate(A,HO,dt,IsThereAnyInput)
% Check sample rate, throw warning if it's inadequate
if IsThereAnyInput
   switch HO % oversampling factors from reference [1]
      case 0, kf = 2.841;
      case 1, kf = 4.016;
      case 3, kf = 1.818;
   end
   maxL = eigs(sparse(A),1,'largestabs'); % use eigenvalues to determine if sample rate is adequate
   fNy = abs(maxL)/2/pi; % Nyquist frequency
   if dt==0, dt = 1/fNy/kf/2; end
   dt0 = 1/fNy/2/kf;
else
   dt0 = 0; % the zero-input response is exact
   kf = Inf;
end
if dt>dt0*1.001 && IsThereAnyInput
   warning('The time step is too large (%g>%g) for hold order "%d": the solution''s precision may be inadequate',dt,dt0,HO),
end
end

function c = aux_getPolyCoeff(dt,u,HO)
% Write polynomial coefficients for [1] eq. 18
switch HO
   case 0 % constant (a.k.a. "previous")
      c = u; % [a0]
   case 1 % linear
      c = [u(1:end-1,:) , diff(u,1,1)./dt]; % [a0,a1]
   case 3 % spline
      p = spline((0:size(u,1)-1)*dt , u.'); % splines have C2 continuity, but pchip could be used to
      c = zeros(p.pieces,p.order*p.dim);
      for n = 1:p.pieces % extract spline coefficients to match columns of Bd in aex_c2d
         tmp = p.coefs(p.dim*(n-1)+(1:p.dim),end:-1:1);
         c(n,:) = tmp(:); % [a0,a1,a2,a3]
      end
end
end

function [Ad,Bd] = aux_c2d(A,B,dt,HO) % [1] eq.(9)
% Discretize state and input matrices for given time step and hold order
Nx = size(A,1); % number of states
Nu = size(B,2); % number of inputs

% Calculate state transition matrix, a.k.a. discretized state matrix (zero-input response)
Ad = expm(full(A)*dt); 

% Calculate the zero-state response matrices for each polynomial order
I = eye(Nx);
Hp = Ad; % seed
Bd = zeros(Nx,Nu*(HO+1));
for p = 0:HO
   F = prod(1:p);
   Hp = A\(Hp-dt^p/F*I); % impulse response's i-th integral
   Bd(:,p*Nu+(1:Nu)) = F*Hp*B; % determine response of n-th polynomial (B comes as B(:,exc))
end
end
