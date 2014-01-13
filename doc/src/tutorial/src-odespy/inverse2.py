"""As inverse1.py, but wrapped in a reusable class."""
import numpy as np
import odespy

class Inverse:
    """
    Compute the inverse function of some given function f(x).
    Method: solve an ODE for the inverse function.
    """
    def __init__(self, f, x0, I=[0,1], h=1E-4, resolution=400):
        self.f, self.I, self.h, self.x0 = f, I, h, x0
        self.resolution = resolution

        # Check that f(x0)=I[0]
        froot = abs(f(x0)-I[0])
        if froot > 1E-3:
            raise ValueError(
                'abs(f(x0)-%s)=%s - not a root of f(x)=I[0]' % \
                (I[0], froot))

    def _rhs(self, x, y):
        dfdx = (self.f(x + self.h) - self.f(x))/self.h
        return 1./dfdx

    def discrete(self):
        """Solve dx/dxy = 1/f'(x), x(I[0])=x0."""
        self.y = np.linspace(self.I[0], self.I[1], self.resolution+1)
        solver = odespy.RungeKutta4(self._rhs)
        solver.set_initial_condition(self.x0)
        self.x, self.y = solver.solve(self.y)
        return self.x, self.y

    def continuous(self):
        if not hasattr(self, 'x') and hasattr(self, 'y'):
            self.discrete()
        from scitools.numpyutils import wrap2callable
        self.g = wrap2callable((self.y, self.x))
        return self.g

    def verify(self):
        """Check that g(f(x)) = x."""
        x = np.linspace(self.I[0], self.I[1], self.resolution+1)
        x2 = self.g(self.f(x))
        return np.abs(x - x2).max()

def _test():

    def f(x):
        return np.sqrt(x)

    inverse = Inverse(f, x0=0, I=[0, 4], resolution=10)
    x, y = inverse.discrete()
    g = inverse.continuous()
    print 'max error:', inverse.verify()

    from matplotlib.pyplot import plot, legend, savefig, show
    g_e = lambda y: y**2  # exact inverse function
    y_ = np.linspace(0, 4, 11)
    plot(y, x, 'r-',
         x, f(x), 'b-',
         y_, g_e(y_), 'go')
    legend(['computed inverse', 'original f(x)', 'exact inverse'])
    savefig('tmppng'); savefig('tmp.pdf')
    show()

if __name__ == '__main__':
    _test()




