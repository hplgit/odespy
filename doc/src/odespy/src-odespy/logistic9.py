import odespy, numpy as np, matplotlib.pyplot as plt

def run(problem, tp, solver):
    method = solver.__class__.__name__

    solver.set_initial_condition(problem.U0)

    u, t = solver.solve(tp)
    solver.u_max = u.max()
    print '%.4f  %s' % (solver.u_max, method)

    if solver.has_u_t_all():
        plt.plot(solver.t_all, solver.u_all)
        print '%s used %d steps (%d specified)' % \
              (method, len(solver.u_all), len(tp))
    else:
        plt.plot(solver.t, solver.u)
    legend.append(method)
    plt.hold('on')

rtol = 1E-6
atol = rtol
npoints = 41
problem = odespy.problems.Logistic(a=2, R=1E+5, A=1)
L = 10
N = 30
tp = np.linspace(0, L, N+1)
min_step = 0.0001

methods = ['DormandPrince', 'BogackiShampine',
           'RKFehlberg', 'Vode', 'RKF45', 'Lsoda']
solvers = [eval('odespy.' + method)(
           problem.f, atol=atol, rtol=rtol,
           min_step=min_step)
           for method in methods]
solvers[1].set(adams_or_bdf='bdf', order=4)  # Vode

legend = []
for solver in solvers:
    run(problem, tp, solver)

plt.plot(tp, problem.u_exact(tp))
legend.append('exact')
plt.legend(legend)
plt.savefig('tmp1.png')
plt.show()

