% Drawing and rotation of drawing
s=linspace(0, 2*pi, 20); % 21 points from zero to 2*pi

% parametric description of a circle
circ =   [cos(s); -sin(s)]; %first row is x, second is y

% The body parts are mostly obtained by stretching the x and/or
% y components (multiplication) or displacing them (addition)
head =   [3*circ(1,:); 4*circ(2,:)];
mouth=   [circ(1,1:10); circ(2,1:10)-1.8]; %half a circle
nose =   [[0 .5 0];  [-.5 -1.5 -1.5]];  % three points, an open triangle
lefteye  = [ circ(1,:)*.7-1  ; circ(2,:)/3+.3 ];
righteye = [ circ(1,:)*.7+1  ; circ(2,:)/3+.3 ];

% The matrices above are put side to side into a big two
% row matrix called pict. The top row is x components of
% points in the picture, the bottom row is y components.
% Between each matrix is the
% column vector [nan nan]' which is the two numbers 
% "Not a Number" in a column.  When plotting, MATLAB
% skips [nan nan]'.  Thus this is a trick to keep the
% pieces of the picture from connecting to each other.
refpict = [ head     [nan nan]' ...
         mouth    [nan nan]' ...
	     nose     [nan nan]' ...
	     lefteye  [nan nan]' ...
	     righteye [nan nan]'] ;
 
refpict(2,:) = refpict(2,:) -10;  % offset

hold off;
plot(refpict(1,:),refpict(2,:))
axis('equal')
 
hold on
theta = pi/4;
R = [ cos(theta)  -sin(theta);
      sin(theta)  cos(theta)  ];
pict = R*refpict;
plot(pict(1,:),pict(2,:)); 
axis('equal')




disp(['End at:    '  datestr(now)])