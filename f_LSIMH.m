function [y,t,x,u,Info] = f_LSIMH(A,B,C,D,cs,dt,x0,N)
% Linear time-invariant system response to complex exponential input
%
% [y,t,x,u] = f_LSIMH(A,B,C,D,cs,dt,x0,N)
%
% INPUT:
% A    [Nx,Nx] state matrix ([]: DefaultSystem)
% B    [Nx,Nu] input matrix ([]:  eye(Nx) )
% C    [Ny,Nu] output matrix ([]: eye(Nx) )
% D    [Ny,Nu] feedthrough matrix ([]: zeros(Ny,Nu) )
% cs   [N,Nu] excitation's FFT coefficients and frequencies ([]: cosine signal on first input);
% dt   [num] time step (1/fs) ([]: automatic selection)
% x0   [Nx,1] initial conditions ([]:zeros(Nx,1))
% N    [num] number of desired time steps ([]: DefaultStepNumber)
%
% OUTPUT:
% y    [N Ny] outputs
% t    [N 1] time vector
% x    [N,Nx] states
% u    [N,Nu] excitation
%
% EXCITATION DEFINITION: 
% LSIM-H.m assumes that the input is in the form u(t) = sum c*exp(s*t).
% *  cs(:,1:end-1) are the complex Fourier series coefficients of the input
% *  cs(:,end) are the corresponding angular frequency bins [rad/s]
% The j-th row of cs corresponds to the j-th input. 
% The Fourier coefficients cs(:,1:end-1) are obtained from the FFT of a
% time-domain excitation matrix u, e.g. cs(:,1:1:end-1) = FFT(u) or from
% known input spectra such as von Kármán continuous turbulence. With the
% excitation in this form, several frequency-domain techniques and filters
% can be applied.
%
% SYNTAX:
% [] = f_LSIMH free decay of 1-DoF harmonic oscillator, wn = 1 rad/s, zeta = 0.01, x0 = [1 0];
% [] = f_LSIMH(A) free decay with x0 = [1 0 ... 0];
% [] = f_LSIMH(A,B) response to cs=cosine on input 1
% [] = f_LSIMH(A,B,C) response to cs=cosine on input 1
% [] = f_LSIMH(A,B,C,D) response to cs=cosine on input 1
% [] = f_LSIMH(...,cs = [1/2i w ; -1/2i -w) input is sin(w*t)
% [] = f_LSIMH(...,cs = [1/2i w ; -1/2i -w ; 1/2 3*w ; 1/2 -3*w ;) input is sin(w*t)+cos(3*w*t)
% [] = f_LSIMH(...,dt) given time step
% [] = f_LSIMH(...,x0) given initial conditions
% [] = f_LSIMH(...,N) given number of time steps
%
% EXAMPLE:
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
% N = 2001; % desired number of time steps
% dt = 1/10; % sample rate for ode45.m
% kf = 10; % downsampling factor for f_LSIMH (solution accuracy independent from sample rate)
% i = 1; j = 1; % response and excitation channels
% 
% % Excitation, initial conditions
% x0 = zeros(size(A,1),1); x0(1) = .5; % initial conditions
% a = 1; % excitation amplitude
% w = .5; % excitation frequency
% % cs = [a/2i 1i*w ; -a/2i -1i*w]; % Fourier coefficients of sine input
% cs = [a/2i -0.01*w+1i*w ; -a/2i -0.01*w-1i*w]; % decaying sine
% t = (0:N-1).'*dt; % time vector
% u = exp(t*cs(:,end).')*cs(:,1:end-1); % excitation vector
% 
% % Compute
% ufcn = griddedInterpolant(t,u,'linear'); % build interpolant for ode45.m
% [t1,x] = ode45(@ssfun,t,x0,odeset('MaxStep',dt),A,B(:,j),ufcn); % integrate states
% y1 = x*C(i,:).' + u*D(i,j).'; % calculate outputs
% [y2,t2] = f_LSIMH(A,B(:,j),C(i,:),D(i,j),cs,kf*dt,x0,round(N/kf));
% 
% % Plot
% plot(t2,y2,'.',t1,y1), grid on, xlabel('Time'), ylabel("Response: "+i+", excitation:"+j), legend('LSIM-H','ode45.m')
%
% function dx = ssfun(t,x,A,B,ufcn)
% u = ufcn(t); % evaluate excitation at integrator time
% dx = A*x + B*u;
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
  cs(:,:) double {mustBeFinite} = [];
  dt(1,1) double {mustBeScalarOrEmpty,mustBeNonnegative} = [];
  x0(:,1) double {mustBeFinite} = [];
  N double {mustBeScalarOrEmpty} = [];
end

%% Set defaults
DefaultStateMatrix = [0 1; -1 , -.04]; % harmonic oscillator, omega = 1 rad/s, zeta = 0.02;
DefaultNumberOfSamples = 501;
DefaultOversampling = 1.05;

%% Check input
if isempty(A), A = DefaultStateMatrix; end
assert(size(A,1)==size(A,2),'f_LSIMH:StateMatrixSquare','The state matrix must be a square matrix')
if isempty(B), B = eye(size(A,1)); end
assert(size(A,1)==size(B,1),'f_LSIMH:InputMatrixSize','The input matrix must have as many rows as there are in the state matrix')
if isempty(C), C = eye(size(A)); end
assert(size(C,2)==size(A,1),'f_LSIMH:OutputMatrixSize','The output matrix must have as many columns as there are in the state matrix')
if isempty(D), D = zeros(size(C,1),size(B,2)); end
assert(size(D,1)==size(C,1) && size(D,2)==size(B,2),'f_LSIMH:FeedthroughMatrixSize','The feedthrough matrix must have as many rows as there are in the output matrix and as many columns as there are in the input matrix')
if isempty(dt)
   fnM = max(abs(eig(A)))/2/pi;
   dt = 1/(fnM*2*DefaultOversampling);
else
   fn = [];
end
if isempty(cs) % set default input
   if isempty(fn), fn = sort(abs(eig(A)))/2/pi; end
   fu = fn(1) + 0.05*(max(fn)-min(fn)); % excite near the first eigenfrequency
   cs = [1/2 zeros(1,size(B,2)-1) 2i*pi*fu ; 1/2 zeros(1,size(B,2)-1) -2i*pi*fu]; % cosine signal
end
assert(size(cs,2)>1,'f_LSIMH:InputMinSize','f_LSIMH requires the FFT coefficients and angular frequency')
assert(size(cs,2)==size(B,2)+1,'f_LSIMH:InputNumber','f_LSIMH requires that there are as many FFT coefficient columns as there are inputs')
if isempty(x0), x0 = zeros(size(A,1),1); end, x0 = x0(:);
assert(size(A,1)==numel(x0),'f_LSIMH:InitialConditionsSize','The initial conditions must be as many as there are states')
if isempty(N), N = DefaultNumberOfSamples; end


%% Algorithm

% Parse excitation
[c,s,t,u] = aux_getInputCoeff(dt,cs,N);

% Discretize system matrices
[Ad,Bd] = aux_c2dH(A,B,dt,s,c);

% Integrate states
x = aux_intstates(Ad,Bd,dt,s,x0,N);

% Calculate outputs
y = x*C.'+u*D.';

% Collect metadata
Info = struct('TimeStep',dt,'HoldOrder',Inf,'InitialConditions',x0,...
   'Algorithm',mfilename,'Duration',t(end),'OversamplingFactor',Inf);

end

function [c,s,t,u] = aux_getInputCoeff(dt,cs,N)
% Gets excitation coefficients, builds time axis and constructs excitation vector
c = cs(:,1:end-1); % FFT coefficients of the excitation
s = cs(:,end); % frequency bins of the excitation
t = (0:N-1).'*dt; % time axis
u = exp(t*s.')*c; % excitation in time domain
end

function [Ad,Bd] = aux_c2dH(A,B,dt,s,c) % [1] eq. (23)
% Discretizes matrices
Nw = numel(s);
Nx = size(A,1);
if ~isdiag(A) % nondiagonal state matrix - slow [1] eq. (23)
   I = eye(Nx);
   Ad = expm(A*dt); % state transition matrix
   Bd = zeros(Nx,Nw); % discretized input matrix
   for k = 1:Nw % for each excitation frequency
      Bd(:,k) = (s(k)*I-A)\(B*c(k,:).');
   end
else % the state matrix is diagonal - fast [1] eq. (49)
   Ad = diag(exp(A(1:Nx+1:Nx^2)*dt)); % state transition matrix
   Bd = zeros(Nx,Nw); % discretized input matrix
   for k = 1:Nx % for each state
      Bd(k,:) = (s-A(k,k)).\(c*B(k,:).');
   end
end
end

function [x,t] = aux_intstates(Ad,Bd,dt,s,x0,N) % [1] eq (26)
% Integrates states
Ad = Ad.'; Bd = Bd.'; % transpose for column-major state vector
Nx = size(Ad,1);

% Marching algorithm
t = (0:N-1).'*dt; % time axis
x = zeros(N,Nx); % initialize state vector
x(1,:) = x0;
xzi = x(1,:) -sum(Bd,1); % initial zero-input state: initial conditions + initial state due to forced transient
if ~isdiag(Ad)
   for n = 1:N-1 % for each time step
      xzi = xzi*Ad; % zero-input response
      xzs = exp(s.'*t(n+1))*Bd; % zero-state response
      x(n+1,:) = xzi + xzs; % state
   end
else
   x = xzi.*diag(Ad).^(0:N-1).' + exp(t*s.')*Bd;
% In some cases a loop may be faster:
%    xzs = zeros(N,Nx);
%    for k = 1:numel(s) % for each excitation frequency
%       xzs = xzs + exp(t*s(k))*Bd(k,:);
%    end
%    for k = 1:Nx % for each state
%       x(:,k) = xzi(:,k)*Ad(k,k).^(0:N-1).' + xzs(:,k);
%    end
end

% For purely imaginary cs, it is possible to use the Chebyshev recursion
% method (digital resonator).
% z1 = zeros(numel(w),1);
% z2 = ones(numel(w),1);
% z3 = exp(-1i*w*dt);
% z0 = 2*real(z3);
% for n = 1:1E6
%    z1 = z0.*z2-z3; 
%    z3 = z2;
%    z2 = z1;
% end % it is faster to use three variables
% Although the errors propagate, the relative error for N = 1E6 is
% norm(z(1,:)-exp(1i*w.'*dt*n))./norm(z(1,:)) < 1E-6
end
