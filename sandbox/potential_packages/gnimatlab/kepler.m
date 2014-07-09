function [out,out2,out3] = kepler(t,q,flag,varargin)
%   Example GNI problem file for solving the 2D Kepler problem.
%   Use it as follows:
%
%       gni_meth('kepler',[],[],[],0.7);
%
%   With gni_meth replaced by either gni_irk2, gni_lmm2, or gni_comp
%   to draw the solution with eccentricity 0.7.
if (nargin < 3) | isempty(flag)
	rad=q(1,:).*q(1,:)+q(2,:).*q(2,:);
	rad=rad.*sqrt(rad);
	out(1,:)=-q(1,:)./rad;
	out(2,:)=-q(2,:)./rad;
else
	switch flag
	case 'init',
		if (nargin < 4)
			ecc = 0.5;
		else
			ecc = varargin{1};
		end
		if (ecc < 0) | (ecc >= 1)
			error('The eccentricity must lie between 0 and 1');
		end
		out = [0 2*pi];
		out2 = [1-ecc,0,0,sqrt((1+ecc)/(1-ecc))];
		out3 = gniset('NumSteps',50,'Vectorized','on','Events','off',...
			'OutputSteps',-1,'OutputFcn','phaseplot');
	end
end

