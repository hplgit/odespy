def f(u, t, a=1, R=1):
    return a*u*(1 - u/R)

A = 1

import odespy, numpy
from matplotlib.pyplot import plot, hold, show, axis

solver = odespy.RK4(f, f_kwargs=dict(a=2, R=1E+5))

# Split time domain into subdomains and
# integrate the ODE in each subdomain
T = [0, 1, 4, 8, 12]     # subdomain boundaries

N_tot = 30               # total no of time steps
dt = float(T[-1])/N_tot  # time step, kept fixed
u = []; t = []           # collectors for u and t in each domain

for i in range(len(T)-1):
    T_interval = T[i+1] - T[i]
    N = int(round(T_interval/dt))
    time_points = numpy.linspace(T[i], T[i+1], N+1)

    solver.set_initial_condition(A)  # at time_points[0]
    print 'Solving in [%s, %s] with %d intervals' % \
          (T[i], T[i+1], N)
    ui, ti = solver.solve(time_points)
    A = ui[-1]  # newest ui value is next initial condition

    plot(ti, ui)
    hold('on')

    u.append(ui);  t.append(ti)

axis([0, T[-1], -0.1E+5, 1.1E+5])
# Can concatenate all the elements of u and t, if desired
u = numpy.concatenate(u);  t = numpy.concatenate(t)
savefig('tmppng'); savefig('tmp.pdf')
show()
#plot(t, u, 'bo')  # same curve, but only one piece
