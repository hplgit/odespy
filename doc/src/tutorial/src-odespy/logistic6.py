import numpy as np
import matplotlib.pyplot as plt
import odespy

class Logistic:
    def __init__(self, a, R, A, T):
        """
        a` is (initial growth rate), `R` the carrying capacity,
        `A` the initial amount of u, and `T` is some (very) total
        simulation time when `u` is very close to the asymptotic
        value `R`.
        """
        self.a, self.R, self.A = a, R, A
        self.tol = 0.01*R # tolerance for termination criterion

    def f(self, u, t):
        """Right-hand side of the ODE."""
        a, R = self.a, self.R  # short form
        return a*u*(1 - u/R)

    def terminate(self, u, t, step_no):
        """u[step_no] holds solution at t[step_no]."""
        return abs(u[step_no] - self.R) < self.tol

    def u_exact(self, t):
        a, R, A = self.a, self.R, self.A  # short form
        return R*A*np.exp(a*t)/(R + A*(np.exp(a*t) - 1))


class Solver:
    def __init__(self, problem, dt, method='RK4'):
        self.problem = problem
        self.dt = dt
        self.method_class = eval('odespy.' + method)
        self.N = int(round(T/dt))

    def solve(self):
        self.solver = self.method_class(self.problem.f)
        self.solver.set_initial_condition(self.problem.A)
        time_points = np.linspace(0, self.problem.T, self.N+1)
        self.u, self.t = self.solver.solve(
            time_points, self.problem.terminate)
        print 'Final u(t=%g)=%g' % (t[-1], u[-1])

    def plot(self):
        plt.plot(self.t, self.u, 'r-',
                 self.t, self.u_exact(self.t), 'bo')
        plt.legend(['numerical', 'exact'])
        plt.savefig('tmp.png'); plt.savefig('tmp.pdf')
        plt.show()

def main():
    problem = Logistic(a=2, R=1E+5, A=1, T=20)
    solver = Solver(problem, dt=0.25, method='RK4')
    solver.solve()
    solver.plot()

if __name__ == '__main__':
    main()
