function odedemo2()
% Ball Falling in honey
% -Andy Ruina, Jan 24, 2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  VARIABLES (Assume consistent units)
%  +y  = displacement up
%   v  = dy/dt
%   z  = [y v],    z is the 'state vector'
% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CONSTANTS
  g  =  10  ;  % gravity constant
  c  =    .1  ;  % drag constant (force per unit velocity)
  m  =    .1;  % mass

% INTIAL CONDITIONS
y0 =  0 ; % initial height
v0 =  20 ; % initial velocity
z0 = [y0 v0];  % z is the list of state variables
               % z is the state vector.

tspan =[0 5]; %time interval of integration

error = 1e-4; 
% Set error tolerance and use 'event detection'
options = odeset('abstol', error, 'reltol',error,...
         'events', @stopevent) ;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask Matlab to SOLEVE odes in function 'rhs'
[t zarray] = ode45(@rhs,tspan, z0, options,m,c,g);
%  The parameters m,c and g are passed to both
%  the 'rhs' function and the 'stopevent' function
% Each row of zarray is the state z at one time.

%UNPACK the zarray (the solution) into sensible variables
y = zarray(:,1);  %  y is the first  column of z
v = zarray(:,2);  %  v is the second column of z

plot (t,y, t, v)
title('Andy''s plot of ball in honey')
xlabel('time'); ylabel('y and v')
axis([0 inf -inf inf])  %inf self scales plot

end  % end of main function odedemo2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE DIFFERENTIAL EQUATION 'The Right Hand Side'
function zdot = rhs(t,z,m,c,g)
%UNPACK state vector z into sensible variables
y = z(1);  %  y is the first  element of z
v = z(2);  %  v is the second element of z

%The equations
ydot = v;          % kinematic relation between y and v
vdot = -g -c*v/m;  %  F = m a  

% Pack the rate of change of y, v into a rate of change
%  of state:
zdot = [ydot vdot]';  % Has to be a column vector
end % end of rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value, isterminal, dir] = stopevent(t,z,m,c,g)
% Have to assign numbers to value, isterminal, dir
%UNPACK z into sensible variables
y = z(1);  %  y is the first  element of z
v = z(2);  %  v is the second element of z
value      = y - (-5);  % stop integrating when y = -5
isterminal = 1;         % 1 means stop
dir        = -1; % -1 for decreasing, +1 for increasing,
                 %  0 for any which way.
end  % end of stopevent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
