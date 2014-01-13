a = 2
R = 1E+5
A = 1

f_f77_str = """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = %.3f*u(1)*(1 - u(1)/%.1f)
      return
      end
""" % (a, R)
print f_f77_str
#import sys; sys.exit(1)

import odespy
f_f77 = odespy.compile_f77(f_f77_str)
solver = odespy.Lsode(f=None, f_f77=f_f77)

solver.set_initial_condition(A)

from numpy import linspace
T = 10  # end of simulation
N = 30  # no of time steps
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)

from matplotlib.pyplot import *

plot(t, u, 'r-')
savefig('tmppng'); savefig('tmp.pdf')
show()
