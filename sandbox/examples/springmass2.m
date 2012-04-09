function springmass2 ()
%Spring-mass with forcing
% mx'' + cx' + kx = F_0 sin (omega t)
m=2;  k = 5; c= 0.2;
F0 = 5;
omega = sqrt(k/m); %force at resonant frequency

tspan=[0 20*pi];
x0 = 0;  v0=0; % initial conditions
z0 = [x0 v0];


[t z_soln] =ode45(@W,tspan, z0,[],m,k,c,omega,F0);

%unpack z_soln, matrix where each row i is the
%value of z at time t(i)
x = z_soln(:,1); % first  column of z_soln
v = z_soln(:,2); % second column

figure(1)
  plot (t, x);
  xlabel('t = time')
  ylabel('x = position')
  title('Andy''s plot')
  grid on
figure(2)
 plot(x,v);
 axis('equal')
  xlabel('x = position')
  ylabel('v = velocity')
  title('Phase plot')
  grid on

end

% Here is the differential equation
function zdot=W(t,z,m,k,c,omega,F0)
%unpack z, list of values of z at time t
x = z(1);
v = z(2);

xdot = v;
vdot = -(k/m)*x - (c/m)*v + F0*sin(omega*t);

zdot = [xdot vdot]';
end