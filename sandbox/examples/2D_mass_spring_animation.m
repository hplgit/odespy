function Animation ()
%Three new commands here:
%  line, set, drawnow
%Good for animations
% This animation shows  a mass hanging from a ‘zero-rest-length’
% spring in 2D.

densityofpoints = 50;  tmax = 25;
tspan = linspace(0, tmax, densityofpoints*tmax);

%Initial conditions
 x0 = 0; vx0 =10; y0 = -1; vy0 = 0;
 z0 = [x0 vx0 y0 vy0];

%parameters
 m = 100; k = 100; g = 10;

% Solve the ODEs
[t,zall] = ode45(@rhs, tspan,z0,[],m,k,g);
    x =  zall(:,1);
    vx = zall(:,2);
    y =  zall(:,3);
    vy = zall(:,4);


figure(1)  % Snapshot of the whole trajectory
plot(x,y);
axis('equal')

%%%%%%%%%%%%
figure(2);  % the animated plot
clf;        % clear figure
axis([-10 10 -20 0]);
grid('on')

% the “line” command
dot1 = line('xdata',x(1),'ydata',y(1), 'marker','.',...
            'markersize',[30],'color','red',...
            'erase','xor'); % the moving dot
dot2 = line('xdata',x(1),'ydata',y(1), 'marker','.',...
            'markersize',[5],'color','green',...
            'erase', 'none'); % the trail
         
%Animation happens in this loop
 for i = 2:length(x)
        set(dot1, 'xdata',x(i),  'ydata',y(i));
        set(dot2, 'xdata',x(i),  'ydata',y(i));
        drawnow;
 end

fprintf('Done making plot.\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zdot = rhs(t,z,m,k,g)
%unpack variables
  x  =z(1);
  vx =z(2);
  y  = z(3);
  vy = z(4);

%governing equations
  xdot = vx;
  vxdot = -(k/m)* x;

  ydot = vy;
 vydot = -(k/m)* y - g;

%pack back in to zdot
zdot = [xdot vxdot ydot vydot]';

end