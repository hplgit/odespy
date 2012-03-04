class Problem:
    stiff = False
    complex_ = False

    parameters = dict(
        T=dict(help='stop time', default=1, type=float)
        U0=dict(help'initial condition', default=0))

    def __init__(self):
        pass

    def f(self, u, t):
        raise NotImpelementedError

    def jac(self, u, t):
        raise NotImpelementedError

    def constraints(self, u, t):
        raise NotImpelementedError


    def define_command_line_arguments(self, parser):
        raise NotImpelementedError

    def set_parameters(self, args)
        raise NotImpelementedError

    def u_exact(self, t):
        pass

    def verify(self, u, atol=1e-6, rtol=1e-5):
        pass

    # subclasses may implement computation of derived
    # quantities, e.g., energy(u, t) etc

class LinearOscillator(Problem):
    pass
# compute potental and kinetic energy, use their sum
# to verify
# offer analytical solution with damping, also with sin/cos excitation
