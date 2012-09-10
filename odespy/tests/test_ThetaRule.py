from numpy.testing import assert_array_almost_equal, TestCase, \
    run_module_suite, assert_almost_equal

import numpy as np
import odespy

class TestThetaRule(TestCase):

    def hand_coded_version(self, I, a, T, dt, theta):
        """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
        dt = float(dt)           # avoid integer division
        N = int(round(T/dt))     # no of time intervals
        T = N*dt                 # adjust T to fit time step dt
        u = np.zeros(N+1)        # array of u[n] values
        t = np.linspace(0, T, N+1)  # time mesh

        u[0] = I                 # assign initial condition
        for n in range(0, N):    # n=0,1,...,N-1
            u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
        return u, t

    def test_odespy(self):
        """
        Compare all versions of odespy methods that correspond to
        the theta-rule with theta=0, 0.5, 1, against the
        hand_coded_version.
        """
        I = 0.8; a = 1.2; T = 4; dt = 0.5
        prms0 = dict(eps_iter=1E-6, verbose=0)

        methods = [(0, odespy.ForwardEuler, {}), # (theta, classname, odespy prms)
                   (0, odespy.ThetaRule, {'theta': 0})]
        for method in 'BackwardEuler', 'CrankNicolson', 0.5, 1:
            if method == 0.5:
                method = 'ThetaRule'
                theta = 0.5
            elif method == 1:
                method = 'ThetaRule'
                theta = 1
            elif method == 'CrankNicolson':
                theta = 0.5
            elif method == 'BackwardEuler':
                theta = 1

            for nlsolver in 'Picard', 'Newton':
                prms = prms0.copy()
                prms['theta'] = theta
                prms['nonlinear_solver'] = nlsolver
                methods.append((theta, eval('odespy.'+method), prms))
                if nlsolver == 'Newton':
                    # add a version with Jacobian
                    def myjac(u, t):
                        return -a
                    #prms['jac'] = lambda u, t: -a
                    prms2 = prms.copy()
                    prms2['jac'] = myjac
                    methods.append((theta, eval('odespy.'+method), prms2))
        for theta, odespy_solver, parameters in methods:
            u_d, t_d = self.hand_coded_version(I, a, T, dt, theta)
            solver = odespy_solver(f=lambda u, t: -a*u, **parameters)
            solver.set_initial_condition(I)
            u, t = solver.solve(np.linspace(0, T, t_d.size))
            ok = np.allclose(u_d, u,
                             rtol=10*prms0['eps_iter'],
                             atol=prms0['eps_iter'])
            # str(solver) prints all parameters *different* from the defaults
            print str(solver), 'works' if ok else 'does not work'
            if not ok:
                print 'max deviation:', abs(u_d-u).max()
            assert ok

if __name__ == "__main__":
    run_module_suite()

