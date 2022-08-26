clc

% State-space model (https://de.mathworks.com/help/control/ug/mimo-state-space-models.html)
A = [-0.0558 -0.9968  0.0802  0.0415
      0.5980 -0.1150 -0.0318  0
     -3.0500  0.3880 -0.4650  0
      0       0.0805  1.0000  0];
B = [ 0.0073       0
     -0.4750  0.0077
      0.1530  0.1430
      0       0     ];
C = [ 0       1       0       0;
      0       0       0       1];
D = zeros(2);

% Parameters
N = 1001; % desired number of time steps
dt = 1/5; % sample rate
i = 1; j = 2; % desired response and excitation channels
x0 = zeros(size(A,1),1); % initial conditions
t = (0:N-1).'*dt; % time vector
u = randn(N,numel(j)); % random excitation

% Compute LSIMP solution for first-order hold
[y0,t0] = f_LSIMP(A,B(:,j),C(i,:),D(i,j),u,dt,x0,1);

% Compute solution with ode45
ufcn = griddedInterpolant(t,u,'linear'); % build interpolant
[t1,x] = ode45(@ssfun,t,x0,odeset('MaxStep',dt),A,B(:,j),ufcn); % integrate states
y1 = x*C(i,:).' + u*D(i,j).'; % compute outputs

% Plot
plot(t1,y1,'-',t0,y0,'--')
ylabel('Response'), xlabel('Time [s]'), grid on
legend('ode45','LSIMH')

function dx = ssfun(t,x,A,B,ufcn)
u = ufcn(t); % evaluate excitation at integrator time
dx = A*x + B*u; % state equation
end
