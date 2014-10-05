import odespy, numpy as np, scitools.std as st
k = 150.
#k = 50000.
# k instability because large k means small s, even simulation
# around equilibrium is problematic with k=50,000, there is
# round-off in s and in final force when k/m*s*ny balances g
L0 = 9.81  # unstretched length of pendulum
g = 9.81
m = 1.0
theta0 = 60.0  # degrees
theta0 = theta0*np.pi/180
# Equilibruim length
L = L0 + m*g/k
# Initial condition: theta degrees from equilibrium
x0 = L*np.sin(theta0)
# cos(theta0) = (L0-L)/y0
y0 = L0 - L*np.cos(theta0)

# Unknowns: vx, x, vy, y
# Equations: vx' = ..., x' = vx, vy' = ..., y' = vy

def f(u, t, m, k, L0):
    g = 9.81
    vx, x, vy, y = u
    # mv' = -k*s*n + m*g*j, n: normal vector along the pendulum
    L = np.sqrt(x**2 + (y-L0)**2)
    s = L - L0
    nx = x/L
    ny = (y-L0)/L
    #print 't=%g vx=%5.2f vy=%5.2f x=%5.2f, y=%5.2f, Fx=%5.2f, Fy=%5.2f' % (t, vx, vy, x, y, -k/m*s*nx, -k/m*s*ny - g)
    #print 'k=%g, s=%5.2g, ny=%5.2f, k/m*s*ny=%g, k*s*ny/m=%g' % (k,s,ny,k/m*s*ny,k*s*ny/m)
    return [-k/m*s*nx, vx,
            -k/m*s*ny - g, vy]

solver = odespy.EulerCromer(f, f_args=(m, k, L0))
solver.set_initial_condition([0, x0, 0, y0])
# First test: vertical pendulum doing harmonic vertical motion, first
# at rest, then slightly out of equilibrium
# k and m should balance

# For large k, this is theta'' + g/L0*sin(theta) = 0, so L0=g
# implies theta'' + theta = 0 equation with theta0*cos(t) as solution
P = 2*np.pi
N = 60
dt = P/N
num_periods = 6
T = num_periods*P
time_points = np.linspace(0, T, num_periods*N+1)
u, t = solver.solve(time_points)
x = u[:,1]
y = u[:,3]
theta = np.arctan(x/(L0-y))
theta_exact = lambda t: theta0*np.cos(t)
theta_e = theta_exact(t)
if abs(theta0) < 1E-14:
    st.figure()
    st.plot(t, y, title='y motion')
st.figure()
# Control perfect aspect ratio of the axis so the motion is geometrically true
st.plot(x, y, 'b-', title='xy motion', daspect=[1,1,1], daspectmode='manual',
        axis=[x.min(), x.max(), L0-L-x.max(), L0-L+x.max()])
st.figure()
st.plot(t, theta, t, theta_e, legend=('EC', 'exact'), title='theta motion')
# Difference in final value for N=60
#diff_exact = 0.0014700112828
#diff_EC = abs(x[-1] - x_e[-1])
#tol = 1E-14
#assert abs(diff_exact - diff_EC) < tol, \
#       'diff_exact=%g, diff_EulerCromer=%g' % (diff_exact, diff_EC)
raw_input()
