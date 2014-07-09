import odespy
from numpy import *
from matplotlib.pyplot import *

def f(u, t):
    omega, theta = u
    return [-c*sin(theta), omega]

c = 1
Theta0_degrees = 30

solver = odespy.EulerCromer(f)
Theta0 = Theta0_degrees*pi/180
solver.set_initial_condition([0, Theta0])
# Solve for num_periods periods using formulas for small theta
freq = sqrt(c)          # frequency of oscillations
period = 2*pi/freq      # one period
N = 40                  # intervals per period
dt = period/N           # time step
num_periods = 10
T = num_periods*period  # total simulation time

time_points = linspace(0, T, num_periods*N+1)
u, t = solver.solve(time_points)

# Extract components and plot theta
theta = u[:,1]
omega = u[:,0]
theta_linear = lambda t: Theta0*cos(sqrt(c)*t)
plot(t, theta, t, theta_linear(t))
legend(['Euler-Cromer', 'Linearized problem'], loc='lower left')
title('Error in linearized term: %8.2g' % abs(Theta0-sin(Theta0)))
savefig('tmp.pdf'); savefig('tmp.png')
show()
