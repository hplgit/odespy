def f(u, t):
    return [u[1], 3.*(1. - u[0]*u[0])*u[1] - u[0]]

u0 = [2.0, 0.0]
import odespy, numpy, matplotlib.pyplot as plt

for method in odespy.Lsode, odespy.DormandPrince, odespy.Vode:
    solver = method(f, rtol=0.0, atol=1e-6,
                    adams_or_bdf='adams', order=10)
    solver.set_initial_condition(u0)
    t_points = numpy.linspace(0, 30, 150)
    u, t = solver.solve(t_points)

    plt.plot(t, u[:,0])
    plt.hold('on')
plt.show()
