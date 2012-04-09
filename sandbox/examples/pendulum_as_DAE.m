function simp_pend
clc
close all
%%% Simple pendulum as a DAE
%%% 21 April 2009 %%%%
%%% Pranav Bhounsule %%%%

l = 1; m = 1; I = m*l^2/12; g=10;  
T = 5; steps = 300;
Tol = 1e-6; 

th0 = 1;
thdot0 = 0;

x0 =  (l/2)*cos(th0);
y0 = (l/2)*sin(th0);

xdot0 = -(l/2)*sin(th0)*thdot0;
ydot0 =  (l/2)*cos(th0)*thdot0;

z0 = [x0 xdot0 y0 ydot0 th0 thdot0]';
tspan = linspace(0,T,steps);

options=odeset('abstol',Tol,'reltol',Tol);

[t,z] = ode45(@f,[tspan],z0,options,l,m,I,g);


pendanimate(t,z,l)

function zdot = f(t,z,l,m,I,g)

x =     z(1);
xdot =  z(2);
y =     z(3);
ydot =  z(4);
th =    z(5);
thdot = z(6);

A = [m      0       0                -1               0; 
     0      m       0                 0              -1; 
     0      0       I            -(l/2)*sin(th)   (l/2)*cos(th); 
     -1     0    -(l/2)*sin(th)       0               0; 
     0     -1     (l/2)*cos(th)       0               0];
b = [ m*g; 
        0;
        0; 
    (l/2)*cos(th)*thdot^2;
    (l/2)*sin(th)*thdot^2];

X = A \ b;
xddot  = X(1);
yddot  = X(2);
thddot = X(3);

zdot = [xdot xddot ydot yddot thdot thddot]';

function pendanimate(t,z,l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

moviescaling = 1;                      % slow down factor
n=length(t); dur=t(end);
delay =floor(moviescaling*1000*dur/n); %delay per frame in .001 secs 

bar1ref = [0 1; 0 0];

r = l/2;
axis('equal')
d = 0.2*r;
axis([-2*r-d 2*r+d -2*r-d 2*r+d]);

bar1pic=line('xdata',bar1ref(1,:),'ydata',bar1ref(2,:), ...
              'linewidth', 3, 'erase','xor','color','red');

hinge=line('xdata',0,'ydata',0, 'marker','.','markersize',[30], ...
           'color','black','erase', 'xor');
			  
for i=1:n
   for j=1:100, log(1:delay*10); end %delay for graphics. 
	                                 %the number in this expression
									 %is machine dependent.
									 %The LOG is just something
									 %to keep the machine busy.

   theta = -pi/2;
   R=[cos(theta) -sin(theta)	; sin(theta) cos(theta)];
   x = z(i,1); y = z(i,3); th = z(i,5);

   bar1 = R*[x-r*cos(th) x+r*cos(th); y-r*sin(th) y+r*sin(th)];
   set(bar1pic,'xdata',bar1(1,:),'ydata',bar1(2,:));

   drawnow
   if i==1
       pause(1)
   end
end
