clc
close all

% https://www.mathworks.com/matlabcentral/fileexchange/68657-ooma-toolbox
m = 10*ones(1,4);
d = 10*ones(1,4);
k = 10000*ones(1,4);

M = diag(m);

K = [k(1)+k(2) -k(2)       0          0 ;...
    -k(2)       k(2)+k(3) -k(3)       0 ;...
     0         -k(3)       k(3)+k(4) -k(4) ;
     0          0         -k(4)       k(4)];

D = [d(1)+d(2) -d(2)       0          0 ;...
    -d(2)       d(2)+d(3) -d(3)       0 ;...
     0         -d(3)       d(3)+d(4) -d(4) ;
     0          0         -d(4)       d(4)];
  
[A,B,C,D] = mdk2ss(M,D,K);

% Parameters
x0 = zeros(size(A,1),1); x0(1) = .0; % initial conditions
N = 2001; % number of time steps for LSIMP
NH = 51; % number of time steps for LSIMH (independent of time step)
T = 5; % duration
i = 1; % outputs
j = 2; % inputs

dt = T/(N-1); % time step
t = (0:N-1).'*dt; % time axis
u = 10000*(cos(12*2*pi/T*t) + sin(20*2*pi/T*t)); % some excitation

% Get frequency-domain coefficients for f_LSIMH. These can be taken as the
% FFT coefficients of u, which is a good approach if u is harmonic. If u
% contains decaying sines, then many FFT coefficients are obtained; a
% better approach, if u is known exactly, is to provide its coefficients
% directly (see help f_LSIMH)
cs = fcn_getFFTcoeffs(dt,u);

% Compute exact solution with LSIMH using a large time step
[y0,t0,~,u0] = f_LSIMH(A,B(:,j),C(i,:),D(i,j),cs,T*2/(2*NH-1),x0,NH);

% Compare solution with LSIMP (approximate)
[y1,t1] = f_LSIMP(A,B(:,j),C(i,:),D(i,j),u,dt,x0,3);

% Compare solution with ode45 (approximate)
ufcn = griddedInterpolant(t,u,'linear'); % build interpolant
[t2,x] = ode45(@ssfun,t,x0,odeset('MaxStep',dt),A,B(:,j),ufcn); % integrate states
y2 = x*C(i,:).'+u*D(i,j); % compute outputs

% Compare
figure
plot(t2,y2,'-',t1,y1,'--',t0,real(y0),'.k')
ylabel('Response'), xlabel('Time [s]'), grid on
legend('ode45','LSIMP','LSIMH')

function cw = fcn_getFFTcoeffs(dt,u)
% Example function for extracting FFT coefficients of input vector
% Exact when u is a harmonic function
TollPeriodicity = 100*eps;
TollAmplitude = 100*eps;

if norm(u(1,:)-u(end,:))/norm(u(1,:))<TollPeriodicity, u(end,:) = []; end
L = size(u,1); % length of excitation ~= N (simulation length)
c = fft(u) / L;
dtu = dt;

% Double-sided frequency vector
if mod(L,2)
   w = 0:(L+1)/2-1;
else
   w = 0:L/2;
end
w = 2*pi/dtu/(L-1)*[w,-w(end+mod(L,2)-1:-1:2)].';

% Remove frequency bins with low amplitudes
idx = all(abs(c)./max(abs(c))<TollAmplitude,2);
c(idx,:) = []; w(idx) = [];

% Output
cw = [c,1i*w];
end

function dx = ssfun(t,x,A,B,ufcn)
u = ufcn(t); % evaluate excitation at integrator time
dx = A*x + B*u; % state equation
end

function [A,B,C,D] = mdk2ss(M,D,K)
Nm = size(M,1);
A = [zeros(Nm),eye(Nm);-M\K,-M\D]; 
B = [zeros(Nm);M\eye(Nm)]; 
C = [eye(Nm),zeros(Nm)]; 
D = zeros(Nm);
end

