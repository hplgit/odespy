function status = sphereplot(t,y,flag)
%   Example of custom output function. It has to be used in 
%   a GNI options structure as 
%
%      gniset('OutputFcn','sphereplot');
%
%   It requires to have 3*N components given. It interprets them
%   as N cartesian coordinate triplets that are supposed to lie
%   on the unit sphere.
status = 0; 
chunk = 128;

if nargin < 3 | isempty(flag)
  ud = get(gcf,'UserData');

  nt = length(t);
  chunk = max(chunk,nt);
  [cols,rows] = size(ud.x);
  
  if ud.i + 1 > rows
    ud.x = [ud.x, zeros(ud.N,chunk)];
    ud.y = [ud.y, zeros(ud.N,chunk)];
    ud.z = [ud.z, zeros(ud.N,chunk)];
  end
  ud.i = ud.i + 1;
  ud.x(1:ud.N,ud.i) = y(1:3:(3*ud.N));
  ud.y(1:ud.N,ud.i) = y(2:3:(3*ud.N));
  ud.z(1:ud.N,ud.i) = y(3:3:(3*ud.N));
  set(gcf,'UserData',ud);
  
  if ud.stop == 1
    status = 1;
    sphereplot([],[],'done');
  else
  	for i=1:ud.N
		set(ud.line(i),'Xdata',ud.x(i,1:ud.i),'Ydata',ud.y(i,1:ud.i),...
			'Zdata',ud.z(i,1:ud.i));
	end
  end
  
else
  switch(flag)
  case 'init'
    ud = [];
    ud.N = floor(size(y) / 3);
    ud.N = ud.N(1);
    ud.x = zeros(ud.N,chunk);
    ud.y = zeros(ud.N,chunk);
    ud.z = zeros(ud.N,chunk);
    ud.i = 1;
    for i=1:ud.N
   		ud.x(i,1) = y(3*i-2);
    	ud.y(i,1) = y(3*i-1);
    	ud.z(i,1) = y(3*i);
    end
    
    f = figure(gcf);
	colormap([0.8 0.8 1]);
	[X,Y,Z] = sphere(100);
	sph = surfl(X,Y,Z,'light');
	set(sph(1),'EdgeColor','none');
	hold on
	colr = ['b','r','g','m','c'];
	for i=1:ud.N
		ud.line(i) = plot3(ud.x(i,1),ud.y(i,1),ud.z(i,1),...
			strcat(colr(mod(i-1,5)+1),'o-'),'EraseMode','none');
	end
	hold off
	set(gca,'XLim',[-1 1]);
	set(gca,'YLim',[-1 1]);
	set(gca,'ZLim',[-1 1]);
    
    % The STOP button.
    h = findobj(f,'Tag','stop');
    if isempty(h)
      ud.stop = 0;
      pos = get(0,'DefaultUicontrolPosition');
      pos(1) = pos(1) - 15;
      pos(2) = pos(2) - 15;
      str = 'ud=get(gcf,''UserData''); ud.stop=1; set(gcf,''UserData'',ud);';
      uicontrol( ...
          'Style','push', ...
          'String','Stop', ...
          'Position',pos, ...
          'Callback',str, ...
          'Tag','stop');
    else
      set(h,'Visible','on');
      if ishold
        oud = get(f,'UserData');
        ud.stop = oud.stop; 
      else
        ud.stop = 0;
      end
    end
    set(f,'UserData',ud);
    
  case 'done' 
    f = gcf;
    ud = get(f,'UserData'); 
    ud.x = ud.x(1:ud.N,1:ud.i);
    ud.y = ud.y(1:ud.N,1:ud.i);
    ud.z = ud.z(1:ud.N,1:ud.i);
    set(f,'UserData',ud);
  	for i=1:ud.N
	set(ud.line(i),'Xdata',ud.x(i,1:ud.i),'Ydata',ud.y(i,1:ud.i)...
		,'Zdata',ud.z(i,1:ud.i));
	end
    set(findobj(f,'Tag','stop'),'Visible','off');
    refresh;
  end
end

drawnow;
