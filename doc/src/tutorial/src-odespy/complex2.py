"""As complex1.py, but comparison of zvode and RK4."""

def f(u, t):
    return 1j*w*u

import odespy, numpy, time

w = 2*numpy.pi
n = 600  # no of periods
r = 40   # resolution of each period
tp = numpy.linspace(0, n, n*r+1)

solvers = [odespy.Vode(f, complex_valued=True,
                       atol=1E-7, rtol=1E-6,
                       adams_or_bdf='adams'),
           odespy.RK4(f, complex_valued=True),
           odespy.RKFehlberg(f, complex_valued=True,
                             atol=1E-7, rtol=1E-6)]
cpu = []
for solver in solvers:
    solver.set_initial_condition(1+0j)
    t0 = time.clock()
    solver.solve(tp)
    t1 = time.clock()
    cpu.append(t1-t0)

# Compare solutions at the end point:
exact = numpy.exp(1j*w*tp).real[-1]
min_cpu = min(cpu); cpu = [c/min_cpu for c in cpu]  # normalize
print 'Exact: u(%g)=%g' % (tp[-1], exact)
for solver, cpu_time in zip(solvers, cpu):
    print '%-15s u(%g)=%.6f (error: %10.2E, cpu: %.1f)' % \
          (solver.__class__.__name__,
           solver.t[-1], solver.u[-1].real,
           exact - solver.u[-1].real, cpu_time)
