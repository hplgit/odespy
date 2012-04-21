import odespy, numpy, matplotlib.pyplot as plt

problem = odespy.problems.Exponential(a=-1, b=0)
solvers = odespy.odelab.solvers
scheme = 'ExplicitEuler'
scheme = 'RungeKutta34'
solver = odespy.odelab(problem.f, odelab_solver=scheme, atol=1E-2, rtol=1E-2)
solver.set_initial_condition(1.0)
u, t = solver.solve(numpy.linspace(0, 3, 21))
print '# pts:', len(t)
plt.plot(t, u)
plt.show()



