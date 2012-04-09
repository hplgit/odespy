function baseball_trajectory
% Calculates the trajectory of a baseball.
% Calculates maximum range for given speed,
% with and without air friction.
% Shows shape of path at high speed.
disp(['Start time:  ' datestr(now)])
cla

% (a) ODEs are in the function rhs far below.
%     The 'event' fn that stops the integration
%     when the ball hits the ground is in 'eventfn'
%     even further below.
% (b) Coefficients for a real baseball taken
% from a google search, which finds a paper
% Sawicki et al, Am. J. Phys. 71(11), Nov 2003.
% Greg Sawicki, by the way, learned some dynamics 
% in TAM 203 from Ruina at Cornell. 

% All parameters in MKS.
m   = 0.145;    % mass of baseball, 5.1 oz
rho = 1.23;     % density of air in kg/m^3
r   = 0.0366;   % baseball radius (1.44 in)
A   = pi*r^2;   % cross sectional area of ball
C_d = 0.35;     % varies, this is typical
g   = 9.81;     % typical g on earth
b   = C_d*rho*A/2; % net coeff of v^2 in drag force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b-d)  Use typical homerun hit speed and look
% at various  angles of hit.

tspan=linspace(0,100,1001); % give plenty of time
n = 45;  % number of simulations 
angle = linspace(1,89,n);  %  launch from 1 to 89 degrees
r0=[0 0]';   % Launch x and y position.

% First case:  No air friction.
b = 0;
subplot(3,2,1)
hold off

% Try lots of launch angles, one simulation for
% each launch angle.
for i = 1:n  
inspeed = 44;  % typical homerun hit (m/s), 98 mph.

theta0 = angle(i)*pi/180; % initial angle this simulation
v0=inspeed*[cos(theta0) sin(theta0)]'; %launch velocity
z0=[r0; v0]; % initial position and velocity

options=odeset('events',@eventfn);
[t zarray]=ode45(@rhs,tspan,z0,options,g,b,m); %Solve ODE

x=zarray(:,1); y=zarray(:,2); %Unpack positions
range(i)= x(end); % x value at end,  when ball hits ground

plot(x,y); title('Jane Cho: Baseball trajectories, no air friction')
xlabel('x, meters'); ylabel('y, meters'); axis('equal')
axis([0 200 0 200])
hold on  % save plot for over-writing
end % end of for loop for no-friction trajectories

%Plot range vs angle, no friction case
subplot(3,2,2); hold off;
plot(angle,range);
title('Range vs hit angle, no air friction')
xlabel('Launch angle, in degrees')
ylabel('Hit distance, in meters')

% Pick out best angle and distance
[bestx besti] = max(range);
disp(['No friction case:'])
best_theta_deg = angle(besti)
bestx


% Second case:  WITH air friction
% Identical to code above but now b is NOT zero.
b   = C_d*rho*A/2; % net coeff of v^2 in drag force

subplot(3,2,3)
hold off % clear plot overwrites

% Try lots of launch angles
for i = 1:n  % 
inspeed = 44;  % typical homerun hit (m/s), 98 mph.

theta0 = angle(i)*pi/180; % initial angle this simulation
v0=inspeed*[cos(theta0) sin(theta0)]'; %launch velocity
z0=[r0; v0]; % initial position and velocity

options=odeset('events',@eventfn);
[t zarray]=ode45(@rhs,tspan,z0,options,g,b,m); %Solve ODE

x=zarray(:,1); y=zarray(:,2); %Unpack positions
range(i)= x(end); % x value at end,  when ball hits ground

plot(x,y); title('Baseball trajectories, with air friction')
xlabel('x, meters'); ylabel('y, meters'); axis('equal')
axis([0 120 0 120])
hold on  % save plot for over-writing
end % end of for loop for with-friction trajectories

%Plot range vs angle, no friction case
subplot(3,2,4); 
plot(angle,range);
title('Range vs hit angle, with air friction')
xlabel('Launch angle, in degrees')
ylabel('Hit distance, in meters')

%Find Max range and corresponding launch angle
[bestx besti] = max(range);
disp(['With Friction:'])
best_theta_deg = angle(besti)
bestx


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now look at trajectories at a variety of speeds
% Try lots of launch angles
subplot(3,2,6)
hold off
speeds = 10.^linspace(1,8,30); % speeds from 1 to 100 million m/s
for i = 1:30  % 
inspeed = speeds(i);  % typical homerun hit (m/s), 98 mph.

theta0 = pi/4; % initial angle is 45 degrees at all speeds
v0=inspeed*[cos(theta0) sin(theta0)]'; %launch velocity
z0=[r0; v0]; % initial position and velocity

options=odeset('events',@eventfn);
[t zarray]=ode45(@rhs,tspan,z0,options,g,b,m); %Solve ODE

x=zarray(:,1); y=zarray(:,2); %Unpack positions
range(i)= x(end); % x value at end,  when ball hits ground

plot(x,y); title('Trajectories, with air friction, various speeds ')
xlabel('x, meters'); ylabel('y, meters'); axis('equal')
axis([0 2000 0 2000])
hold on  % save plot for over-writing
end % end of for loop for range at various speeds

disp(['End time:  ' datestr(now)])
end % end of Baseball_trajectory.m




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Governing Ord Diff Eqs.
function zdot=rhs(t,z,g,b,m)
% Unpack the variables
x=z(1); y=z(2);
vx=z(3); vy=z(4);

%The ODEs
xdot=vx; ydot=vy; v = sqrt(vx^2+vy^2);
vxdot=-b*vx*v/m;
vydot=-b*vy*v/m - g;

zdot= [xdot;ydot;vxdot;vydot]; % Packed up again.
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Event' that ball hits the ground
function [value isterminal dir] = eventfn(t,z,g,b,m)
y=z(2);
value = y;      % When this is zero, integration stops
isterminal = 1; % 1 means stop.
dir= -1;        % -1 means ball is falling when it hits
end

