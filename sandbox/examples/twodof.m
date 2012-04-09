function twodof
%TWO MASSES IN  A LINE WITH A SPRING BETWEEN THEM
%AND A SPRING FROM EACH TO THE WALL (3 SPRINGS)
m1 = 1;   m2 = 2;  
k1 = 1; k2 = 1; k3 = 1;
M = [m1 0 ; 0 m2]; 
%NOTE SIGNS
K = [ (k1 + k2)     -k2    ;
          -k2    (k2 + k3)  ];

tspan = linspace(0,10, 1001);
xzero = [1   -1]'; % INITIAL DISPLACEMENT
 vzero = [0   0 ]';
zzero = [xzero; vzero];

[t zmat] = ode45(@rhs, tspan, zzero,[],M, K);

x = zmat(:,[1 2]);
v = zmat(:,[3 4]);

plot(t,x(:,1), t, x(:,2))
xlabel('time'); ylabel('position(s)')
title('Motion of the two masses')
disp(['It is now  '  datestr(now) ' ,   try again?'])
end


function zdot = rhs(t,z,M,K)
x=z([1 2]);
v=z([3 4]);
xdot = v;

vdot = M^(-1)* (-K) *x;

zdot = [xdot; vdot];
end