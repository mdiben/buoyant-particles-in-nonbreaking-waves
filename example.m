%% Example code:

clear all
close all

%% initialize parameters

% Wave parameters %%%%%%%%%%%%%%%%%%%%%%%
ka=0.15; %wave steepness
A=1; %(m) wave amplitude
g=9.81; %(m/s^2) acceleration due to gravity
k=ka/A; %(1/m) wave number
omega=sqrt(g*k); %(rad/s) wave frequency
Per=2*pi/omega; %(s) wave period

% Mixing parameters %%%%%%%%%%%%%%%%%%%%%
K=0.01; %(m^2/s) eddy viscosity
L_m=1; %(m) mixing length

% Particle parameters %%%%%%%%%%%%%%%%%%%
ws=K/L_m; %(m/s) rise velocity

% Simulation parameters %%%%%%%%%%%%%%%%%
N_particles=100; %number of particles
N_waves=500; %number of wave periods
t=[0 N_waves*Per]; %time start and end
dt=Per/50; %time step
T=[t(1):dt:t(end)];

% initial x locations of particles, randomly distributed in phase
X0s=rand(N_particles,1)*2*pi/k;

% initialize output vectors
N_steps=length(T);
XX=nan(N_steps,N_particles);
ZZ=XX;
zeta=XX;


%% loop through N_particles
for i=1:N_particles
    
    % start particles at approximate surface:
    eta=A*cos(k*X0s(i));
    X0i=[X0s(i) 0 eta];
    
    % run tracking code
    [T,X,s] = waves_turb_settling(X0i,t,dt,ka, A, K, K, ws);
    
    % put into output vectors
    XX(:,i)=X(:,1);
    ZZ(:,i)=X(:,3);
    zeta(:,i)=s;
    
    
end

% calculate relative wave phase for Stokes 3rd waves
phase=-rem((k*XX-(1+ka^2/2)*omega*T'),2*pi);


%% Make some simple output plots
% need to run many particles for a long time for statistics to converge.
% in this example, i'm using few particles but taking statistics from a
% subset of time to show example histograms with shorter computational cost

%% vertical distributions
figure(1)

subplot(2,2,1)
scatter(phase(end,:),ZZ(end,:)/L_m,'.')
lim=ylim;
xlabel('\theta')
ylabel('z/L_m')
axis([0 2*pi lim])

subplot(2,2,2)
histogram(ZZ(end-500:end,:)/L_m,'orientation','horizontal','normalization','pdf')
axis([xlim lim])
xlabel('probability density')
ylabel('z/L_m')

subplot(2,2,3)
scatter(phase(end,:), zeta(end,:)/L_m,'.')
lim=ylim;
xlabel('\theta')
ylabel('\zeta/L_m')
axis([0 2*pi lim])

subplot(2,2,4)
hold off
histogram(zeta(end-500:end,:)/L_m,'orientation','horizontal','normalization','pdf')
axis([xlim lim])
xlabel('probability density')
ylabel('\zeta/L_m')
hold all
plot(exp([-10:.1:0]/L_m),-10:.1:0,'k--') %no waves theory

%% horizontal distribution
figure(2)

histogram(phase(end-1000:end,:),linspace(0,2*pi,25),'normalization','pdf')
ylabel('probability density')
xlabel('\theta')

%% joint horizontal and vertical distribution
figure(3)

histogram2(phase(end-1000:end,:),zeta(end-1000:end,:),'normalization','pdf','displaystyle','tile','linestyle','none')
colormap(hot)
colorbar
ylabel('\zeta/L_m')
xlabel('\theta')
axis([0 2*pi -3 0])


