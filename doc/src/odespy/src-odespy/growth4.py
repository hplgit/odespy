"""As exgr1.py, but several consequtive calls to solve."""
c = 0.1
A = 1.5

def f(u, t):
    return c*u

import odespy, numpy
from matplotlib.pyplot import *

solver = odespy.RK4(f)

# Split time domain into subdomains and
# integrate the ODE in each subdomain
T = [0, 3, 6, 10, 20, 40]

N_tot = 30               # no of time intervals in total
dt = float(T[-1])/N_tot  # time step, kept fixed
u = []; t = []           # collectors for u and t

for i in range(len(T)-1):
    T_interval = T[i+1] - T[i]
    N = int(round(T_interval/dt))
    time_points = numpy.linspace(T[i], T[i+1], N+1)

    solver.set_initial_condition(A)  # at time_points[0]
    print 'Solving in [%s, %s] with %d intervals' % \
          (T[i], T[i+1], N)
    ui, ti = solver.solve(time_points)
    A = ui[-1]  # newest u is next initial condition

    plot(ti, ui)
    hold('on')

    u.append(ui);  t.append(ti)

# Can concatenate all the elements of u and t, if desired
u = numpy.concatenate(u);  t = numpy.concatenate(t)
#plot(t, u, 'bo')  # same curve
show()
