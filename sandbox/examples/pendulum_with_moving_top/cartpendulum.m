function cartpendulum
clc
close all

y0 = 0.0; 
ydot0 = 0;
theta0 = pi/2;
thetadot0 = 0;

Tol = 1e-6;

z0 = [y0 ydot0 theta0 thetadot0]';

tspan=linspace(0,10,300);

options=odeset('abstol',Tol,'reltol',Tol);

[t,z] = ode113(@cartpend,[tspan],z0,options);


animate_cart_pend(t,z)

%plot(t,z(:,:));

function zdot = cartpend(t,z)

m1 = 1.5; 
m2 = 1; 
l  = 1; 
g  = 10;

y        =  z(1);
ydot     =  z(2);
theta    =  z(3);
thetadot =  z(4);


MM = [m1+m2,        m2*l*cos(theta);
      cos(theta),   l              ];
  
RHS = [m2*l*thetadot^2*sin(theta);
       -g*sin(theta)             ];

X = MM \ RHS;

yddot     = X(1); 
thetaddot = X(2);

zdot = [ydot yddot thetadot thetaddot]';

function animate_cart_pend(t,z)

%close all
clf
%%%
x = z(:,1);
q = z(:,3);


a1 = 0.4;
a2 = 0.4;
b1 = 0.2;
b2 = 0.2;
r = b1/4;
c1 = a1/2;
c2 = a2/2;
% x = linspace(-0.6,1.2,200);
y = (2*r+b1)*ones(1,length(t));


%%%%%%%%%%%%%%%%%%
%%% B ---------  C
%%%  |      b2   |
%%%  |  a1  G a2 |
%%%  |      b1   |
%%% A----------- D
%%%%%% R1    R2 %%%
%x_min = min([3*min(x) -3*min(x)]); x_max = max([3*max(x) -3*max(x)]);
x_min = -1.2; x_max = 1.2;
y_min = -0.9; y_max = 1.3;

axis([x_min x_max y_min y_max]);
axis('equal');
%axis off
%Set(gcf,'Color',[1,1,1])

%bar6ref=[x_min-1 x_max+1; 0 0];
bar6ref=[-10 10; 0 0];
bar6pic=line('xdata',bar6ref(1,:),'ydata',bar6ref(2,:), ...
            'linewidth', 1,'color','black');


moviescaling = 2;                      % slow down factor
n=length(t); dur=t(end);
delay =floor(moviescaling*1000*dur/n); %delay per frame in .001 secs 


bar1ref = [0 1;0 0];
bar2ref = [0 1;0 0];
bar3ref = [0 1;0 0];
bar4ref = [0 1;0 0];
bar5ref = [0 0;0 -0.8];
bar10ref = [0 1;0 0];

bar1pic=line('xdata',bar1ref(1,:),'ydata',bar1ref(2,:), ...
             'linewidth', 2, 'erase','xor','color','red');
bar2pic=line('xdata',bar2ref(1,:),'ydata',bar2ref(2,:), ...
 'linewidth', 2, 'erase','xor','color','red');
bar3pic=line('xdata',bar3ref(1,:),'ydata',bar3ref(2,:), ...
 'linewidth', 2, 'erase','xor','color','red');
bar4pic=line('xdata',bar4ref(1,:),'ydata',bar4ref(2,:), ...
 'linewidth', 2, 'erase','xor','color','red');
bar5pic=line('xdata',bar5ref(1,:),'ydata',bar5ref(2,:), ...
 'linewidth', 2, 'erase','xor','color','blue');

bar10pic=line('xdata',bar10ref(1,:),'ydata',bar10ref(2,:), ...
 'linewidth', 1, 'erase','xor','color','white');
%text1 = text(0,1,' Thank You','FontSize',18,'erase','normal'); 

hingepic=line('xdata',0,'ydata',0, 'marker','.','markersize',[50], ...
          'erase','xor','color','blue');

theta=linspace(0,2*pi,100);

diskix=r*cos(theta); diskiy=r*sin(theta); diskiz=linspace(0,0,100);
disk = patch(diskix,diskiy,diskiz,'green');
disk2 = patch(diskix,diskiy,diskiz,'green');

for (i=1:length(x)) 
    
   for j=1:100, log(1:delay*10); end %delay for graphics. 
	                                 %the number in this expression
									 %is machine dependent.
									 %The LOG is just something
									 %to keep the machine busy.
    
    x0 = x(i); y0 = y(i);
    th = q(i);
    A = [x0-a1 y0-b1];
    B = [x0-a1 y0+b2];
    C = [x0+a2 y0+b2];
    D = [x0+a2 y0-b1];
    R1 = [x0-c1 y0-(b1+r)]; 
    R2 = [x0+c2 y0-(b1+r)];
    bar1=[B(1) C(1); B(2) C(2)];
    bar2=[C(1) D(1); C(2) D(2)];
    bar3=[D(1) A(1); D(2) A(2)];
    bar4=[A(1) B(1); A(2) B(2)];
    %bar10=[A(1) C(1); A(2) C(2)];
    
    R=[cos(th) -sin(th)	; sin(th) cos(th)]; % one rotation matrix
    bar5=R*bar5ref;
    
    
    set(bar1pic,'xdata',bar1(1,:),'ydata',bar1(2,:));
    set(bar4pic,'xdata',bar4(1,:),'ydata',bar4(2,:));
    set(bar2pic,'xdata',bar2(1,:),'ydata',bar2(2,:));
    set(bar3pic,'xdata',bar3(1,:),'ydata',bar3(2,:));
    
    set(bar5pic,'xdata',bar5(1,:)+x0,'ydata',bar5(2,:)+y0);
    
    set(hingepic,'xdata',bar5(1,2)+x0,'ydata',bar5(2,2)+y0);
    
    %set(bar10pic,'xdata',bar10(1,:),'ydata',bar10(2,:));
    %set(text1,'Position',[x0-0.2 y0]);
    
    diskx = diskix+R1(1)*ones(1,length(theta)); disky = diskiy+R1(2)*ones(1,length(theta)); diskz = diskiz; 
    set(disk,'xdata',diskx,'ydata',disky,'zdata',diskz);
    
    diskx2 = diskix+R2(1)*ones(1,length(theta)); disky2 = diskiy+R2(2)*ones(1,length(theta)); diskz2 = diskiz; 
    set(disk2,'xdata',diskx2,'ydata',disky2,'zdata',diskz2);
     
    drawnow
    
    if i==1
    pause(2);
    end
end
