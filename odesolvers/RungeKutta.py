from ODE import Solver, Adaptive
import numpy as np

class RungeKutta(Adaptive):
    """
    Super class for explicit 1-level- Runge-Kutta methods:
    RungeKutta4, Rungekutta2,etc..
    or 2-levels adaptive Runge-Kutta methods:
    Dormand & Prince, Cash-Karp or Fehlberg, etc.

    Available subclasses are:
    RungeKutta2, RungeKutta3, ForwardEuler, DormandPrince,
    RungeKutta4, Fehlberg, CashKarp, BogackiShampine,
    MyRungeKutta (user-supplied RungeKutta methods),
    """
    _method_order = None
    # An integer for unadaptive method,
    # or a pair of integers for 2 levels in adaptive methods.

    _butcher_tableau = None
    # (n, n) array for unadaptive methods,
    # or (n+1, n) array for adaptive ones.

    _optional_parameters = Adaptive._optional_parameters + \
        ['min_step','max_step','first_step']

    def get_order(self):
        '''
        Return the order of current method, both for non-adaptive
        and adaptive methods.
        '''
        order = getattr(self, 'method_order', None)
        if order is None:       # User-supplied method in MyRungeKutta
            coefficients = self.butcher_tableau
            if coefficients.shape[0] == coefficients.shape[1]:
                # Butcher_table is a square array
                order = self.calculate_order_1_level(coefficients)
            else:
                # 2-level method
                # Seperate & extract coefficients for two levels
                table_1, table_2 = coefficients[:-1,], \
                    np.vstack((coefficients[:-2,],coefficients[-1,]))
                # Calculate order seperately
                order = self.calculate_order_1_level(table_1),\
                        self.calculate_order_1_level(table_2)
            print '''
        The order of user-defined method is not provided.
        The value of calculated order is %s''' % str(order)
        return order

    def calculate_order_1_level(self,coefficients):
        '''
        Calculate order of 1-level RungeKutta method
        with help of the known solution u = -e**t.

        "Coefficients" is a square 2d-array(butcher tableau).
        '''
        test = MyRungeKutta(lambda u,t:-u,\
                            butcher_tableau = coefficients)
        test.set_initial_condition(1.)
        u,t = test.solve([0.,.1,1.1])

        # Calculate errors with 2 time step size (dt1,dt2).
        # where error1 = O(dt1**order), error2 = O(dt2**order)
        # Then error1/error2 = O(dt1**order)/O(dt2**order)
        # Take logarithms of both side, order can be estimated.
        error1, error2 = \
            u[-1] - np.exp(-t[-1]), u[-2] - np.exp(-t[-2])
        order = int(np.log(abs(error1/error2))/np.log(t[-1]/t[-2]))
        return order

    def advance(self):
        """
        Advance the solution, either using one step or a
        series of adaptive steps, from t[n] tp t[n+1].
        """

        f, n, rtol, atol, neq = \
            self.f, self.n, self.rtol, self.atol, self.neq
        u_n, t_n, t_next = self.u[n], self.t[n], self.t[n+1]
        dt = t_next - t_n

        # default setting of step size
        min_step, max_step, first_step = \
            getattr(self, 'min_step', dt/1000.), \
            getattr(self, 'max_step', dt), \
            getattr(self, 'first_step', dt)

        def middle(x,y=.1,z=4.):    # Auxilary function
            return sorted([x,y,z])[1]

        # Extract coefficients from Butcher-tableau
        table = self._butcher_tableau
        k_len = table.shape[1] - 1   # number of internal stages

        # coefficients for internal stages
        factors_u = np.asarray(table[:k_len, 1:])
        # coefficients for t
        factors_t = table[:k_len, 0]
        # coefficients for u_new
        factors_u_new = table[k_len, 1:]

        if len(table) == k_len + 2:
            # (k+2, k+1) array -->  2-levels Adaptive methods

            # coefficients for local error between 2 levels
            factors_error = table[k_len+1, 1:] - factors_u_new

            u_intermediate = [u_n,]
            t_intermediate = [t_n,]
            u, t, h = u_n, t_n, first_step               # initial values
            k = np.zeros((k_len, self.neq), self.dtype)  # intern stages
            if not hasattr(self, 'u_all'):
                self.u_all = []  # store list of arrays for *all* steps
            if not hasattr(self, 't_all'):
                self.t_all = []  # store list of arrays for *all* steps

            # Loop until next time point is reached
            while (abs(t - t_n) < abs(t_next - t_n)):
                u, t = u_intermediate[-1], t_intermediate[-1]
                # internal steps
                k[:, :] = 0.   # initialization for next step
                for m in range(k_len):
                    k_factors = (np.dot(factors_u, k))[m]
                    #print u, u+h*k_factors, f(u+h*k_factor, 0.5), self.dtype
                    k[m] = f(u+h*k_factors, t+h*factors_t[m])
                u_new = u + h*(np.dot(factors_u_new, k))

                # local error between 2 levels
                error = h*np.abs(np.dot(factors_error, k))
                # Acceptable error tolerance
                tol = rtol*np.abs(u_new) + atol

                accurate = (error <= tol).all()

                if (accurate or h <= min_step or h == max_step):
                    # Accurate enough,
                    # or the step size exceeds valid range.
                    u_intermediate.append(u_new)
                    t_intermediate.append(t+h)

                # prevent the error of dividing absolute zero
                error = np.asarray([(1e-16 if x == 0. else x) \
                                    for x in error])

                # Normarized error rate
                rms = error/tol
                rms_norm = np.sqrt(np.sum(rms*rms)/self.neq)

                order = float(self._method_order[0])
                # factor to adjust the size of next step
                # Formula is from <Numerical Methods for Engineers,
                #  Chappra & Cannle>
                s = .8 *((1./rms_norm)**(1/order))
                # scalar should in range(0.1, 4.)
                # for better accuracy and smoothness
                s = middle(s)
                h *= s

                # step size should in range(min_step, max_step)
                h = middle(h, y=min_step, z=max_step)
                # in case for last intern step
                h = min(h, t_next - t_intermediate[-1])

            self.u_all.append(np.array(u_intermediate))
            self.t_all.append(np.array(t_intermediate))

        else:
            # Run algorithm for explicit 1-level RungeKutta method
            k = np.zeros((k_len, self.neq), float)  # intern stages
            for m in range(k_len):
                k_factors = (np.dot(factors_u, k))[m]
                k[m] = f(u_n + dt*k_factors,t_n + dt*factors_t[m])
            u_new = u_n + dt*(np.dot(factors_u_new, k))
        return u_new

class RungeKutta2(RungeKutta):
    """
    Standard Runge-Kutta method of order 2.
    Implementation in the general RungeKutta Python framework.
    """
    quick_description = "Explicit 2nd-order Runge-Kutta method"

    _butcher_tableau = np.array(\
        [[0., 0., 0.],
         [.5, .5, 0.],
         [0., 0., 1.]])
    _method_order = 2

class RungeKutta3(RungeKutta):
    """
    Standard Runge-Kutta method of order 3.
    Implementation in the general RungeKutta Python framework.
    """
    quick_description = "Explicit 3rd-order Runge-Kutta method"

    _butcher_tableau = np.array(\
        [[0., 0., 0., 0.],
         [.5, .5, 0., 0.],
         [1., -1., 2., 0.],
         [0., .16666667, .66666667, .16666667]])
    _method_order = 3

class RungeKutta1(RungeKutta):
    """
    Explicit Forward Euler method implemented
    in the general RungeKutta Python framework.
    """
    quick_description = "Explicit 1st-order Runge-Kutta method"

    _butcher_tableau = np.array(\
        [[0., 0.],
         [0., 1.]])
    _method_order = 1

class DormandPrince(RungeKutta):
    """
    Dormand&Prince Runge-Kutta method of order (5, 4).
    Implementation in the general RungeKutta Python framework.
    """
    quick_description = "Dormand & Prince RK method of order (5, 4)"

    _butcher_tableau = np.array(\
        [[0., 0., 0., 0., 0., 0., 0., 0.],
         [.2, .2, 0., 0., 0., 0., 0., 0.],
         [.3, .075, .225, 0., 0., 0., 0., 0.],
         [.8, .97777778, -3.73333333, 3.55555556, 0., 0., 0., 0.],
         [.88888889,2.95259869,-11.59579332,9.82289285,-.29080933, 0., 0., 0.],
         [1.,2.84627525,-10.75757576,8.90642272,.27840909,-.2735313, 0., 0.],
         [1.,.09114583, 0.,.4492363,.65104167,-.32237618,.13095238, 0.],
         [0.,.09114583, 0.,.4492363,.65104167,-.32237618,.13095238, 0.],
         [0.,.08991319, 0.,.45348907,.6140625,-.27151238,.08904762,.025]])
    _method_order = (5,4)

class RungeKutta4(RungeKutta):
    """
    Standard Runge-Kutta method of order 4.
    Implementation in the general RungeKutta Python framework.
    """
    quick_description = "Explicit 4th-order Runge-Kutta method"

    _butcher_tableau = np.array(\
        [[0., 0., 0., 0., 0.],
         [.5, .5, 0., 0., 0.],
         [.5, 0., .5, 0., 0.],
         [1., 0., 0., 1., 0.],
         [0., .16666667, .33333333, .33333333, .16666667]])
    _method_order = 4

class Fehlberg(RungeKutta):
    """
    Adaptive Runge-Kutta-Fehlberg method of order (4, 5).
    Implementation in the general RungeKutta Python framework.
    """
    quick_description = "Adaptive Runge-Kutta-Fehlberg (4,5) method"

    _butcher_tableau = np.array(\
        [[0., 0., 0., 0., 0., 0., 0.],
         [.25, .25, 0., 0., 0., 0., 0.],
         [.375, .09375, .28125, 0., 0., 0., 0.],
         [.92307692, .87938097, -3.27719618, 3.32089213, 0., 0., 0.],
         [1., 2.03240741,-8., 7.17348928,-.20589669, 0., 0.],
         [.5, -.2962963, 2., -1.38167641, .45297271, -.275, 0.],
         [0., .11574074, 0., .54892788, .53533138, -.2, 0.],
         [0., .11851852, 0., .51898635, .50613149, -.18, .03636364]])
    _method_order = (4,5)

class CashKarp(RungeKutta):
    """
    Adaptive Cash-Karp Runge-Kutta method of order (5, 4).
    Implementation in the general RungeKutta Python framework.
    """
    quick_description = "Adaptive Cash-Karp RK method of order (5, 4)"

    _butcher_tableau = np.array(
        [[0., 0., 0., 0., 0., 0., 0.],
         [.2, .2, 0., 0., 0., 0., 0.],
         [.3, .075, .225, 0., 0., 0., 0.],
         [.6, .3, -.9, 1.2, 0., 0., 0.],
         [1., -.2037037, 2.5, -2.59259259, 1.2962963, 0., 0.],
         [.875, .0294958, .34179688, .04159433, .40034541, .06176758, 0.],
         [0., .0978836, 0., .40257649, .21043771, 0., .2891022],
         [0., .10217737, 0., .3839079, .24459274, .01932199, .25]])
    _method_order = (5,4)

class BogackiShampine(RungeKutta):
    """
    Adaptive Bogacki-Shampine Runge-Kutta method of order (3, 2).
    Implementation in the general RungeKutta Python framework.
    """
    quick_description = "Adaptive Bogacki-Shampine RK method of order (3, 2)"

    _butcher_tableau = np.array(
        [[0., 0., 0., 0., 0.],
         [.5, .5, 0., 0., 0.],
         [.75, 0., .75, 0., 0.],
         [1., .22222222, .33333333, .44444444, 0.],
         [0., .22222222, .33333333, .44444444, 0.],
         [0., .29166667, .25, .33333333, .125]])
    _method_order = (3,2)

class MyRungeKutta(RungeKutta):
    """
    User-supplied RungeKutta method, which is defined by providing
    butcher-table in an 2d-array.
    Method order should be provided if it is known. If not, the order
    would be estimated automatically with function get_order().
    """
    _butcher_tableau = None
    _method_order = None

    _required_parameters = RungeKutta._required_parameters + \
                           ['butcher_tableau',]

    _optional_parameters =  RungeKutta._optional_parameters + \
                           ['method_order',]

    def validate_data(self):
        if not RungeKutta.validate_data(self):
            return False

        # Check for dimension of user-defined butcher table.
        array_shape = self.butcher_tableau.shape
        if len(array_shape) is not 2:
            raise ValueError,'''
        Illegal input! Your input butcher_tableau should be a 2d-array!'''
        else:
            m,n = array_shape
            if m not in (n, n + 1):
                raise ValueError, '''\
        The dimension of 2d-array <method_yours_array> should be:
        1. Either (n, n), --> For 1-level RungeKutta methods
        2. Or (n+1, n),   --> For 2-levels RungeKutta methods
        The shape of your input array is (%d, %d).''' % (m,n)
        self._butcher_tableau = self.butcher_tableau

        # Check for user-defined order,
        # which should be an integer or a pair of adjacent integers
        if hasattr(self,'method_order'):
            error_1level = '''
        method_order should be a single integer, with a square butcher
        table, which implies a single-level RungeKutta method.
        Your input is %s .''' % str(self.method_order)
            error_2level = '''
        method order should be a pair of adjacent positive integers,
        with a supplied non-square butch table, which implies a
        2-level method. Your input is %s.''' % str(self.method_order)
            if array_shape[0] == array_shape[1] + 1:
                # 2-level RungeKutta methods
                if type(self.method_order) is int:
                    raise ValueError, error_2level
                try:
                    order1, order2 = self.method_order
                    if abs(order1-order2) != 1 or \
                            order1 < 1 or order2 < 1:
                        raise ValueError, error_2level
                except:
                    raise ValueError,error_2level
            else:
                # 1-level RungeKutta methods
                if type(self.method_order) is not int or \
                        self.method_order < 1:
                    raise ValueError,error_1level
            self._method_order = self.method_order

        else:   # method_order is not specified
            if array_shape[0] == array_shape[1] + 1:
                # Calculate order for 2-level-methods
                # Method_order is required for computation
                self._method_order = self.get_order()

        # check for consistency requirement of Butcher Tableau
        for i in range(1,array_shape[1] - 1):
            if not np.allclose(self.butcher_tableau[i][0],\
                               sum(self.butcher_tableau[i][1:])):
                raise ValueError, '''
        Inconsistent data in Butcher_Tableau!
        In each lines of stage-coefficients, first number should be
        equal to the sum of other numbers.
        That is, for a butcher_table with K columns,
            a[i][0] == a[i][1] + a[i][2] + ... + a[i][K - 1]
            where 1 <= i <= K - 1
        Your input for line %d is :%s
        ''' % (i,str(self.butcher_tableau[i]))

        return True
