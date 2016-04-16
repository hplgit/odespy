import odespy, numpy, matplotlib.pyplot as plt

# Run through all scheme classes in with a simple ODE u'=-u, u(0)=1
# and see which classes that work

problem = odespy.problems.Exponential(a=-1, b=0)

scheme = 'ExplicitEuler'
test = odespy.odelab(problem.f, odelab_solver=scheme)
schemes = test.odelab_schemes  # all scheme classes

failures = []
for scheme in schemes:
    print scheme,
    solver = odespy.odelab(problem.f, odelab_solver=scheme,
                           atol=1E-2, rtol=1E-2)
    solver.set_initial_condition(1.0)
    try:
        u, t = solver.solve(numpy.linspace(0, 3, 221))
        print 'summary:', len(u), u[-1], numpy.exp(-t[-1]) - u[-1]
    except:
        failures.append(scheme)

print '\nFailed classes:\n', ' '.join(failures)
print '\nAll solver classes:\n', ' '.join(schemes)
print '\nWorking classes:\n', ' '.join([scheme for scheme in schemes if scheme not in failures])
print '\nDocumented working classes:\n', ' '.join(odespy.odelab.solvers)

