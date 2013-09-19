import odespy, numpy as np

center_points = [2, 3, 4, 5]
s = 0.5
npoints = 41
min_step = 0.0001

results = []
for center_point in center_points:
    problem = odespy.problems.Gaussian1(c=center_point, s=s)
    tp = np.linspace(0, 2*center_point, npoints)
    for tol in (1E-4, 1E-6, 1E-10, 1E-12):
        atol = rtol = tol
        # Build list of solvers
        adaptive = ['RKFehlberg', 'Fehlberg', 'DormandPrince',
                   'CashKarp', 'BogackiShampine', 'Lsoda',] # 'AdaptiveResidual']
        solvers = [(eval('odespy.' + method)
                   (problem.f, atol=atol, rtol=rtol, min_step=min_step),
                    method)
                   for method in adaptive]
        stiff_nonstiff = ['Vode', 'Lsode']
        solvers += [(eval('odespy.' + method)
                    (problem.f, atol=atol, rtol=rtol, min_step=min_step,
                     adams_or_bdf='adams'), method+'_adams')
                    for method in stiff_nonstiff]
        solvers += [(eval('odespy.' + method)
                    (problem.f, atol=atol, rtol=rtol, min_step=min_step,
                     adams_or_bdf='bdf'), method+'_bdf')
                    for method in stiff_nonstiff]

        for solver, method in solvers:

            solver.set_initial_condition(problem.U0)

            u, t = solver.solve(tp)
            u_max = u.max()
            solver.u_max = u_max   # store key result

            print '%s: u_max=%g c/s=%g %s' % \
                  (method, solver.u_max, center_point/s, str(solver))

            results.append((method, center_point, tol, u_max))

# Make results2[method][center_point][tol] = u_max
from collections import OrderedDict
results2 = OrderedDict()
for method, center_point, tol, u_max in results:
    if method not in results2:
        results2[method] = {}
    if center_point not in results2[method]:
        results2[method][center_point] = {}
    if tol not in results2[method][center_point]:
        results2[method][center_point][tol] = {}
    results2[method][center_point][tol] = u_max

print '\n\n======================================================='
for method in results2:
    for center_point in results2[method]:
        for tol in results2[method][center_point]:
            print '%s, c=%g, tol=%E, u_max=%.6f' % \
                  (method, center_point, tol,
                   results2[method][center_point][tol])
#        print '\n'
#    print '\n\n'

