# Author: Liwei Wang

"""
Test for dlsoibt(). Demo from opkddemos.f  in ODEPACK

Simplified Galerkin Solution of Burgers Equation

This program solves a semi-discretized form of the following system
Of three PDEs (each similar to a Burgers equation):

      u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i)    (i=1,2,3),
          t                           x                xx
on the interval  -1 .le. x .le. 1, and with time t .ge. 0.
The diffusion coefficients are eta(*) = .1, .02, .01.
The boundary conditions are u(i) = 0 at x = -1 and x = 1 for all i.
The initial profile for each u(i) is a square wave:
      u(i) = 0         on 1/2 .lt. abs(x) .le. 1
      u(i) = amp(i)/2  on abs(x) = 1/2
      u(i) = amp(i)    on 0 .le. abs(x) .lt. 1/2
      where the amplitudes are amp(*) = .2, .3, .5.

A simplified Galerkin treatment of the spatial variable x is used,
with piecewise linear basis functions on a uniform mesh of 100
intervals.  The result is a system of ODEs in the discrete values
u(i,k) approximating u(i)  (i=1,2,3) at the interior points
(k = 1,...,99).  The ODEs are:
     (u(i,k-1) + 4 u(i,k) + u(i,k+1))/6  =
       -(1/6dx) (c(k-1)dul(i) + 2c(k)(dul(i)+dur(i)) + c(k+1)dur(i))
       + (eta(i)/dx**2) (dur(i) - dul(i))     (i=1,2,3,  k=1,...,99),
      where  c(j) = u(1,j)+u(2,j)+u(3,j),   dx = .02 = the interval size,
      dul(i) = u(i,k) - u(i,k-1),   dur(i) = u(i,k+1) - u(i,k).
Terms involving boundary values (subscripts 0 or 100) are dropped
from the equations for k = 1 and k = 99 above.
"""
from odespy import *
import scitools.std as st
import numpy as np

def res(y, t, s, ires):
    eodsq = [250, 50, 25]
    r = np.zeros(297, float)
    y.shape = r.shape = s.shape = (99, 3)
    cc, cr = sum(y[0]), sum(y[1])
    r[0] = -25./3.*(2.*cc*y[1] + cr*(y[1] - y[0])) + eodsq*(y[1] - 2*y[0])
    for k in range(1, 98):
        cl, cc, cr = sum(y[k-1]), sum(y[k]), sum(y[k+1])
        dli, dri = y[k] - y[k-1], y[k+1] - y[k]
        r[k] = -25./3.*(cl*dli + 2.*cc*(dli + dri) + cr*dri) + eodsq*(dri - dli)
    cl,cc = sum(y[97]), sum(y[98])
    r[98] = -25./3.*(cl*(y[98] - y[97]) - 2.*cc*y[97]) - eodsq*(2*y[98] - y[97])
    if ires!=-1:   # A depend on y
        r[0] -= (s[0]*4. + s[1])/6.
        r[1:98] -= (s[0:97] + 4*s[1:98] + s[2:99])/6.
        r[98] -= (s[97]+s[98]*4.)/6.
    r.shape = (297,)
    return r, ires

def adda(y, t, pa, pb, pc):
    for i in range(3):
        pa[i, i, :] += 4./6.
        pb[i][i][:98] += 1./6.
        pc[i][i][1:] += 1./6.
    return pa, pb, pc

def jac(y, t, s):
    y.shape = (99, 3)
    eodsq = [250, 50, 25]
    pa, pb, pc = np.zeros(891,float), np.zeros(891,float), np.zeros(891,float)
    pa.shape = pb.shape = pc.shape = (3,3,99)
    cc, cr = sum(y[0]), sum(y[1])
    terma, termb = 25.*cr/3, -25.*(2.*cc + cr)/3.
    paij, pbij = -50.*y[1]/3.,-25.*(y[1] - y[0])/3.
    for j in range(3):
        pa[:,j,0], pb[:,j,0] = paij[j], pbij[j]
        pa[j][j][0] += terma - 2.*eodsq[j]
        pb[j][j][0] += termb + eodsq[j]
    for k in range(1,98):
        cl,cc,cr = sum(y[k-1]),sum(y[k]),sum(y[k+1])
        terma,termb,termc = \
                          25.*(cr-cl)/3.,-25.*(2.*cc+cr)/3.,25*(2.*cc+cl)/3.
        dlj, drj = y[k]-y[k-1], y[k+1]-y[k]
        paij,pbij,pcij = -50.*(dlj+drj)/3., -25.*drj/3., -25.*dlj/3.
        for j in range(3):
            pa[:,j,k], pb[:,j,k], pc[:,j,k] = paij[j], pbij[j], pcij[j]
            pa[j][j][k] += terma - 2.*eodsq[j]
            pb[j][j][k] += termb + eodsq[j]
            pc[j][j][k] += termc + eodsq[j]
    cl, cc = sum(y[97]), sum(y[98])
    terma, termc = -25.*cl/3., 25.*(2.*cc + cl)/3.
    paij, pcij = 50.*y[97]/3., -25.*(y[98] - y[97])/3.
    for j in range(3):
        pa[:,j,98], pc[:,j,98] = paij[j], pcij[j]
        pa[j][j][98] += terma - 2.*eodsq[j]
        pc[j][j][98] += termc + eodsq[j]
    y.shape = (297,)     # necessary to recover the shape of array y
    return pa,pb,pc    


# (u1_exact,u2_exact,u3_exact) contains the approximate value when t = 0.4
u1_exact = np.array(
    [1.70956682e-03, 3.43398445e-03, 5.18783349e-03, 6.98515842e-03,
     8.83921016e-03, 1.07622016e-02, 1.27650806e-02, 1.48573251e-02,
     1.70467655e-02, 1.93394396e-02, 2.17394852e-02, 2.42490773e-02,
     2.68684152e-02, 2.95957660e-02, 3.24275691e-02, 3.53586054e-02,
     3.83822285e-02, 4.14906520e-02, 4.46752791e-02, 4.79270545e-02,
     5.12368132e-02, 5.45956048e-02, 5.79949684e-02, 6.14271460e-02,
     6.48852271e-02, 6.83632267e-02, 7.18561029e-02, 7.53597274e-02,
     7.88708192e-02, 8.23868545e-02, 8.59059616e-02, 8.94268082e-02,
     9.29484864e-02, 9.64703968e-02, 9.99921344e-02, 1.03513375e-01,
     1.07033760e-01, 1.10552783e-01, 1.14069668e-01, 1.17583246e-01,
     1.21091827e-01, 1.24593066e-01, 1.28083828e-01, 1.31560049e-01,
     1.35016617e-01, 1.38447256e-01, 1.41844451e-01, 1.45199401e-01,
     1.48502033e-01, 1.51741065e-01, 1.54904135e-01, 1.57977973e-01,
     1.60948623e-01, 1.63801670e-01, 1.66522463e-01, 1.69096305e-01,
     1.71508595e-01, 1.73744902e-01, 1.75790974e-01, 1.77632682e-01,
     1.79255895e-01, 1.80646319e-01, 1.81789276e-01, 1.82669470e-01,
     1.83270725e-01, 1.83575716e-01, 1.83565712e-01, 1.83220322e-01,
     1.82517279e-01, 1.81432251e-01, 1.79938706e-01, 1.78007835e-01,
     1.75608540e-01, 1.72707519e-01, 1.69269456e-01, 1.65257378e-01,
     1.60633244e-01, 1.55358941e-01, 1.49398029e-01, 1.42718981e-01,
     1.35301474e-01, 1.27148627e-01, 1.18308730e-01, 1.08905085e-01,
     9.91559295e-02, 8.93515884e-02, 7.97824293e-02, 7.06663514e-02,
     6.21244732e-02, 5.41994827e-02, 4.68848207e-02, 4.01465202e-02,
     3.39357642e-02, 2.81954415e-02, 2.28635569e-02, 1.78750916e-02,
     1.31630892e-02, 8.65933391e-03, 4.29480447e-03])
u2_exact = np.array(
    [7.17416019e-06, 1.70782645e-05, 3.31245126e-05, 6.01588363e-05,
     1.05339286e-04, 1.79174771e-04, 2.96719122e-04, 4.78862606e-04,
     7.53598916e-04, 1.15707860e-03, 1.73420412e-03, 2.53849668e-03,
     3.63099110e-03, 5.07800919e-03, 6.94782549e-03, 9.30645443e-03,
     1.22130079e-02, 1.57152366e-02, 1.98459102e-02, 2.46205841e-02,
     3.00370492e-02, 3.60764461e-02, 4.27057301e-02, 4.98809820e-02,
     5.75510102e-02, 6.56607602e-02, 7.41541974e-02, 8.29764928e-02,
     9.20754824e-02, 1.01402468e-01, 1.10912474e-01, 1.20564094e-01,
     1.30319039e-01, 1.40141489e-01, 1.49997326e-01, 1.59853293e-01,
     1.69676126e-01, 1.79431680e-01, 1.89084097e-01, 1.98595037e-01,
     2.07923034e-01, 2.17023055e-01, 2.25846345e-01, 2.34340694e-01,
     2.42451240e-01, 2.50121934e-01, 2.57297724e-01, 2.63927433e-01,
     2.69967170e-01, 2.75383917e-01, 2.80158840e-01, 2.84289739e-01,
     2.87792167e-01, 2.90698875e-01, 2.93057586e-01, 2.94927384e-01,
     2.96374262e-01, 2.97466488e-01, 2.98270390e-01, 2.98847025e-01,
     2.99249945e-01, 2.99524080e-01, 2.99705593e-01, 2.99822450e-01,
     2.99895431e-01, 2.99939301e-01, 2.99963931e-01, 2.99975129e-01,
     2.99974996e-01, 2.99961526e-01, 2.99927041e-01, 2.99854809e-01,
     2.99712769e-01, 2.99442742e-01, 2.98942676e-01, 2.98038511e-01,
     2.96441259e-01, 2.93684573e-01, 2.89040478e-01, 2.81421884e-01,
     2.69315148e-01, 2.50874185e-01, 2.24457680e-01, 1.89885662e-01,
     1.49894358e-01, 1.09927672e-01, 7.54041273e-02, 4.90259517e-02,
     3.06080023e-02, 1.85165524e-02, 1.09104125e-02, 6.27726960e-03,
     3.53002680e-03, 1.94049735e-03, 1.04218859e-03, 5.45964314e-04,
     2.77379128e-04, 1.33343739e-04, 5.32660444e-05])
u3_exact = np.array(
    [1.86765383e-10, 1.96772458e-09, 1.19111389e-08, 5.54964761e-08,
     2.18340713e-07, 7.55899524e-07, 2.35604385e-06, 6.70801745e-06,
     1.76224112e-05, 4.30351929e-05, 9.82592148e-05, 2.10736217e-04,
     4.26209304e-04, 8.15657041e-04, 1.48160943e-03, 2.56186555e-03,
     4.22851247e-03, 6.68078970e-03, 1.01317466e-02, 1.47903961e-02,
     2.08424987e-02, 2.84336008e-02, 3.76573037e-02, 4.85502549e-02,
     6.10936693e-02, 7.52198901e-02, 9.08218891e-02, 1.07763660e-01,
     1.25889931e-01, 1.45034247e-01, 1.65025016e-01, 1.85689556e-01,
     2.06856371e-01, 2.28356037e-01, 2.50021072e-01, 2.71685149e-01,
     2.93181998e-01, 3.14344301e-01, 3.35002907e-01, 3.54986687e-01,
     3.74123404e-01, 3.92241969e-01, 4.09176451e-01, 4.24772089e-01,
     4.38893320e-01, 4.51433444e-01, 4.62324969e-01, 4.71549073e-01,
     4.79142163e-01, 4.85197409e-01, 4.89859810e-01, 4.93314543e-01,
     4.95770115e-01, 4.97439231e-01, 4.98520996e-01, 4.99187563e-01,
     4.99576941e-01, 4.99791928e-01, 4.99903753e-01, 4.99958343e-01,
     4.99983239e-01, 4.99993785e-01, 4.99997902e-01, 4.99999367e-01,
     4.99999835e-01, 4.99999965e-01, 4.99999995e-01, 5.00000000e-01,
     5.00000000e-01, 4.99999997e-01, 4.99999976e-01, 4.99999863e-01,
     4.99999315e-01, 4.99996914e-01, 4.99987300e-01, 4.99951740e-01,
     4.99829328e-01, 4.99435130e-01, 4.98245007e-01, 4.94883400e-01,
     4.86081966e-01, 4.65174923e-01, 4.21856650e-01, 3.47885738e-01,
     2.49649938e-01, 1.51648615e-01, 7.80173239e-02, 3.47983164e-02,
     1.38686441e-02, 5.05765688e-03, 1.71052539e-03, 5.38966324e-04,
     1.57923694e-04, 4.27352191e-05, 1.05512005e-05, 2.33068621e-06,
     4.45404604e-07, 6.88336884e-08, 7.23875975e-09])



import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 5   # default number of time-steps    

t0, tn = 0., .4
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-6, 1e-6
mb, nb = 3, 99

u0 = np.zeros(297,float).reshape(99,3)
u0[24], u0[74] = [.1, .15, .25], [.1, .15, .25]
u0[25:74,0], u0[25:74,1], u0[25:74,2] = .2,.3,.5
u0.shape = (297,)


st.figure()
method = Lsoibt

# Test case 1: Lsoibt, with res, adda, jac, mb, nb
m = method(rtol=rtol, atol=atol, res=res, adda_lsoibt=adda, 
           mb=mb, nb=nb, jac_lsoibt=jac)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'g-', title="Lsoibt with Python functions",
        legend="with res, jac & adda", hold="on")

u_final = u[-1].reshape(99,3)
u1, u2, u3 = u_final[:, 0], u_final[:, 1], u_final[:, 2]
max_error = max(max(u1 - u1_exact), max(u2 - u2_exact), max(u3 - u3_exact))
print 'Max error with Test case 1 is %g' % max_error

# Test case 2: Lsoibt, with res, adda, mb, nb
m = method(rtol=rtol, atol=atol, res=res, adda_lsoibt=adda, 
           mb=mb, nb=nb)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r*', title="Lsoibt with Python functions",
        legend="with res & adda", hold="on")

u_final = u[-1].reshape(99,3)
u1, u2, u3 = u_final[:, 0], u_final[:, 1], u_final[:, 2]
max_error = max(max(u1 - u1_exact), max(u2 - u2_exact), max(u3 - u3_exact))
print 'Max error with Test case 2 is %g' % max_error
