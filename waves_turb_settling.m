function [T,X,s] = waves_turb_settling(X0,t,dt,ka, A, K_h, K_z, ws)
%% Lagrangian particle tracking of particles in a wavy, turbulent flow with a settling/rise velocity

% This code integrates one particle trajectory in time starting at location X0, over
% the time interval t, with the time step dt, for waves with specified ka
% and A, turbulence with specified K_h K_z, and particles with
% rise/settling velocity ws.

% The particles follow the fluid flow velocity (in this case waves), plus a random
% walk based on an eddy viscosity and also a constant settling velocity.

% Written by Michelle H. DiBenedetto, last edited Feb 2020
% See the manuscript: Non-breaking wave effects on buoyant particle distributions in
% Frontiers of Marine Science, 2020, for details.

%INPUTS
%X0                     %length 3 vector: initial position of particle [x,y,z] (m)
%t                      %length 2 vector: start and end of time to integrate [t1, tend] (s)
%dt                     %time step increment, delta t (s)
%ka                     %wave steepness
%A                      %wave amplitude (m)
%K_h                    %eddy diffusivity in the horizontal (m^2/s)
%K_z                    %eddy diffusivity in the vertical (m^2/s)

%OUTPUTS
%T                      %time vector
%X                      %position matrix [x y z];
%s                      %distance from the instantaneous free surface vector

%% Set up
%make sure input is a column
if isrow(X0); X0=X0'; end

%wave parameters
kaA={ka A};

%time vector
t=t(1):dt:t(2);

%initialize output vector
NN=length(t);
Y=nan(NN,3);
Y(1,:)=X0;


%% Loop through time
for jj=1:NN-1
    tt=t(jj);
    y=Y(jj,:)';
    
    
    %% particle tracking from waves using Runge Kutta Fourth Order
    
    [ k1,~] = deepwaveflowfun_stokes( y,tt, kaA );
    k1=k1+[0; 0; 0];
    
    y2=y+k1*dt/2;
    
    tk=tt+dt/2;
    [ k2,~] = deepwaveflowfun_stokes( y2,tk, kaA );
    k2=k2+[0; 0; 0];
    
    y3=y+k2*dt/2;
    
    [ k3,~] = deepwaveflowfun_stokes( y3,tk, kaA );
    k3=k3+[0; 0; 0];
    
    y4=y+k3*dt;
    
    tk=tt+dt;
    [ k4,~] = deepwaveflowfun_stokes( y4,tk, kaA );
    k4=k4+[0; 0; 0];
    
    u=(k1+2*k2+2*k3+k4)/6;
    
    y=y+u*dt;
    
    
    %% turbulence using random walk
    RR=randn(1); %gaussian distributed
    x_turb=RR*sqrt(2*K_h*dt);
    
    RR=randn(1);
    z_turb=(RR*sqrt(2*K_z*dt));
    
    y(1:3)=y(1:3)+[x_turb;0;z_turb];
    
    [eta] = deepwaveflowfun_stokes_eta(y(1),tt+dt, kaA);
    
    %% if the particle is above the free surface, reflect below it
    if round(y(3)-eta,8)>=0
        y(3)= eta-(y(3)-eta);
    end
    
    %% next add buoyancy
    y(3)=y(3)+ws*dt;
    
    %% if the particle is above the free surface, put it at the free surface
    if round(y(3)-eta,8)>=0
        y(3)= eta;
    end
    
    Y(jj+1,:)=y;
    s(jj+1)=y(3)-eta;
    
    
end
X=Y;
T=t;

end


function [ u,eta] = deepwaveflowfun_stokes( X,t, kaA )
%this function gives the wave orbital velocities (u), free surface level, eta, using stokes 3rd order
%deep water waves, given position (X), time (t), and wave characteristics (kaA).

ka=kaA{1};
A=kaA{2};
k=ka/A;
g=9.81;
omega=sqrt(g*k);


if A~=0
    omegaprime=(1+ka^2/2)*omega;
    
    U=omega*A;
    x=X(1);
    z=X(3);
    theta=k*x-omegaprime*t;
    C=cos(theta);
    S=sin(theta);
    E=exp(k*z);
    
    
    u=[U.*E.*C; 0;U.*E.*S];
    
    C2=cos(2*theta);
    C3=cos(3*theta);
    
    eta=A*((1-ka^2/16)*C+1/2*ka*C2+3/8*ka^2*C3);
else
    u=[0;0;0];
    eta=0;
end

end




function [eta] = deepwaveflowfun_stokes_eta( X,t, kaA )
%this function gives the free surface level, eta, using stokes 3rd order
%deep water waves, given position (X), time (t), and wave characteristics
%(kaA).

ka=kaA{1};
A=kaA{2};
k=ka/A;
g=9.81;
omega=sqrt(g*k);

if A~=0
    omegaprime=(1+ka^2/2)*omega;
    
    theta=k*X-omegaprime*t;
    C=cos(theta);
    
    C2=cos(2*theta);
    C3=cos(3*theta);
    
    eta=A*((1-ka^2/16)*C+1/2*ka*C2+3/8*ka^2*C3);
else
    eta=0*X;
end

end


