import odespy, numpy as np

atol = 1E-5
rtol = atol
center_points = [2, 3, 4, 5]
s = 0.5
npoints = 41

results = []
for center_point in center_points:
    problem = odespy.problems.Gaussian0(c=center_point, s=s)
    tp = np.linspace(0, 2*center_point, npoints)
    for min_step in (1E-2, 1E-3, 1E-4, 1E-5):
        # Build list of solvers
        adaptive = ['RKFehlberg', 'Fehlberg', 'DormandPrince',
                   'CashKarp', 'BogackiShampine', 'Lsoda',] # 'AdaptiveResidual']
        solvers = [(eval('odespy.' + method)
                   (problem.f, atol=atol, rtol=rtol, min_step=min_step),
                    method)
                   for method in adaptive]
        standard = ['RK4', 'RK2']
        solvers += [(eval('odespy.' + method)(problem.f), method)
                    for method in standard]
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

            print 'min_step: %g, default: %g' % \
                  (min_step, 0.01*(tp[1]-tp[0]))


            solver.set_initial_condition(problem.U0)

            u, t = solver.solve(tp)
            u_max = u.max()
            solver.u_max = u_max   # store key result

            print 'u_max=%g c/s=%g %s' % \
                  (solver.u_max, center_point/s, str(solver))

            results.append((method, center_point, min_step, u_max))

            # Plot the best solver
            name = 'DormandPrince'
            if solver.__class__.__name__ == name:
                import matplotlib.pyplot as plt
                plt.figure()
                if solver.has_u_t_all():
                    plt.plot(solver.t, solver.u, 'r-',
                             solver.t_all, solver.u_all, 'bo',
                             solver.t, problem.u_exact(solver.t), 'g-')
                    plt.legend([name, 'adaptive points', 'exact'])
                else:
                    plt.plot(solver.t, solver.u, 'r-',
                             solver.t, problem.u_exact(solver.t), 'g-')
                    plt.legend([name, 'exact'])
                plt.savefig('tmp_%g.png' % center_point)
plt.show()

# results2[method][center_point][min_step] = u_max
results2 = {}
for method, center_point, min_step, u_max in results:
    if method not in results2:
        results2[method] = {}
    if center_point not in results2[method]:
        results2[method][center_point] = {}
    if min_step not in results2[method][center_point]:
        results2[method][center_point][min_step] = {}
    results2[method][center_point][min_step] = u_max
