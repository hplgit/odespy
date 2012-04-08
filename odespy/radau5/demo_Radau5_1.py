from odespy import *
import scitools.std as st

def f(u,t):
    udot = np.zeros(8, float)
    udot[0] = -1.71*u[0] + .43*u[1] + 8.32*u[2] + .0007
    udot[1] = 1.71*u[0] - 8.75*u[1]
    udot[2] = -10.03*u[2] + .43*u[3] + .035*u[4]
    udot[3] = 8.32*u[1] + 1.71*u[2] - 1.12*u[3]
    udot[4] = -1.745*u[4] + .43*u[5] + .43*u[6]
    udot[5] = -280.*u[5]*u[7] + .69*u[3] + 1.71*u[4] - \
        .43*u[5] + .69*u[6]
    udot[6] = 280.*u[5]*u[7] - 1.81*u[6]
    udot[7] = -udot[6]
    return udot

def jac(u,t):
    dfu = [[-1.71,.43,8.32,0.,0.,0.,0.,0.],
           [1.71,-8.75,0.,0.,0.,0.,0.,0.],
           [0.,0.,-10.03,.43,.035,0.,0.,0.],
           [0.,8.32,1.71,-1.12,0.,0.,0.,0.],
           [0.,0.,0.,0.,-1.745,.43,.43,0.],
           [0.,0.,0.,.69,1.71,-.43-280.*u[7],.69,-280.*u[5]],
           [0.,0.,0.,0.,0.,280.*u[7],-1.81,280.*u[5]],
           [0.,0.,0.,0.,0.,-280.*u[7],1.81,-280.*u[5]]]
    return dfu

u0 = np.zeros(8, float)
u0[0], u0[-1] = 1., .0057
time_points = np.linspace(0., 10., 20)

print "HIRES, chemical reaction, mildly stiff"
st.figure()

# Loop for all possible solvers
for solver in ['Vode', 'Radau5Explicit']:
    method = eval(solver)(f, jac=jac)
    method.set_initial_condition(u0)
    u,t = method.solve(time_points)
    st.plot(t, u[:,0], hold="on", legend=solver, axis=[0.,10.,0.,1.])
    print 'Succeed when solver is %s' % solver
