from odespy import *
from problems import *
from scitools.std import plot, hold
try:
    import joblib  # need joblib for storing database nested dict
except ImportError:
    print 'This script requires joblib to be installed.'
    sys.exit(1)

newresults = {}

def test_Linear():
    # Note: MidpointIter needs max_iter=10 and eps_iter=1E-6;
    # lower tolerances do not give errors to machine precision.
    # No problem with Linear2t so it is the nonlinearity that
    # fools MidpointIter.

    # Problems 2>a>0: RKF45 gives wrong results, error 0.67
    # RungeKutta3 and BogackiShampine error approx 1E-8

    # Problems a<0: numerous, also for Linear2t which is linear

    kwargs = dict(a=1.2, b=2, c=2)
    problems = Linear1(**kwargs), Linear2(**kwargs), Linear2t(**kwargs)
    r = tester(problems, list_available_solvers(), compare_tol=1E-14)
    return r

def test_disk():
    problem = Linear2(a=1.2,b=2,c=2)
    solver = ForwardEuler(problem.f, disk_storage=True, verbose=0)
    solver.set_initial_condition(problem.U0)
    tp = np.linspace(0, 10, 1000000)
    # tp = problem.default_parameters()['time_points']
    u, t = solver.solve(tp)
    print 'made u', u.size
    u0 = u[:,0]
    print 'made u0'
    if tp.size > 500000:
        u = u[-100:]
        t = t[-100:]
        error = problem.verify(u, t)
        print error
    else:
        error = None
    print 'u memmap:', u[:10]
    u[:] = 5.5
    u.flush()
    print u[:10]
    return {'Linear2': {'t': t, 'ForwardEuler': (u, error)}}

def test_Exponential():
    problem = Exponential(a=-0.5, b=2, A=1.5)
    # test f_args, jac_kwargs etc
    methods = [Vode]
    for method in methods:
        solver = method()

def test_VanDerPol(mu=0):
    problem = VanDerPolOscillator(mu=mu, U0=[1,0])
    d = problem.default_parameters()
    tp = d['time_points']
    # test f_args, jac_kwargs etc
    methods = [Vode, RK4, RungeKutta4, ForwardEuler, BackwardEuler]
    for method in methods:
        name = method.__name__
        print name
        solver = method(problem.f, jac=problem.jac,
                        atol=d['atol'], rtol=d['rtol'])
        solver.set_initial_condition(problem.U0)
        u, t = solver.solve(tp)
        plot(t, u[:,0], legend=name,
             legend_fancybox=True, legend_loc='upper left')
        hold('on')
        e = problem.verify(u, t)
        if e is not None: print e


if __name__ == '__main__':
    resfile = 'test_problems_odespy.dat'
    if os.path.isfile(resfile):
        database = joblib.load(resfile)
    else:
        database = {}
    # Run tests
    #test_VanDerPol(float(sys.argv[1]))
    #r = test_Linear()
    r = test_disk()
    import pprint
    f = open('tmp.2', 'w')
    f.write(pprint.pformat(r))
    f.close()
    # We come here only if all tests pass
    joblib.dump(database, resfile)
