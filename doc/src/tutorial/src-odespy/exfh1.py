"""The FitzHugh-Nagumo model from biology."""

def f(u, t, s=10, a=0.12, c1=0.175, c2=0.03, b=0.011, d=0.55):
    v, w = u
    fv = s*c1*v*(v - a)*(1 - v) - s*c2*w
    fw = s*b*(v - d*w)
    return [fv, fw]

def jac(u, t, s=10, a=0.12, c1=0.175, c2=0.03, b=0.011, d=0.55):
    v, w = u
    return [[(1+a)*2*v*c1 - a*c1 - 3*v**2*c1, c2],
            [b, -b*d]]


v0 = 0.1
w0 = 0.0
T = 50

import odespy, numpy as np, scitools.std as st
a = 0.12   # stable
a = -0.12  # unstable
solver = odespy.ForwardEuler(f, jac=jac, f_kwargs={'a': a})
solver.set_initial_condition([v0, w0])
time_points = np.linspace(0, T, 401)
#time_points = np.linspace(0, T, 21)
u, t = solver.solve(time_points)
v = u[:,0]
w = u[:,1]
st.plot(t, v, 'r-',
        t, w, 'b-',
        legend=('v', 'w'))


