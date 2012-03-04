function [out,out2,out3] = twobodysphere(t,q,flag)
%   Example GNI problem file for solving the two-body problem
%   constrained to the sphere of radius one.
%   Use it as follows:
%
%       gni_comp('twobodysphere');
%
%   Only gni_comp works, since this problem uses a special basic
%   method and is therefore not written in the standard way. This problem
%   also uses the custom output function SPHEREPLOT.
if (nargin < 3) | isempty(flag)
	prod = q(1:3)'*q(4:6);
	out = -q([4:6,1:3])/(1-prod^2)^(3/2);
else
	switch flag
	case 'init',
		out = [0 10];
		
		phi = [1.3 -2.1];
		theta = [2.1 -1.1];
		out2([1 4]) = cos(phi).*sin(theta);
		out2([2 5]) = sin(phi).*sin(theta);
		out2([3 6]) = cos(theta);
		
		dphi = [1.2 0.1];
		dtheta = [0.1 -0.5];
		out2([7 10]) = -dphi.*sin(phi).*sin(theta) + dtheta.*cos(phi).*cos(theta);
		out2([8 11]) = dphi.*cos(phi).*sin(theta) + dtheta.*sin(phi).*cos(theta);
		out2([9 12]) = -dtheta.*sin(theta);

		out3 = gniset('StepSize',0.02,'Vectorized','off','Events','off',...
			'PartialFlow','rattwo','OutputFcn','sphereplot',...
			'OutputSteps',5,'OutputSel',[1:6]);
	end
end

