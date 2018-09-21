import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
def f(r,IC,p):
    l = p[0]
    m = p[1]
    om = p[2]
#   
    R = IC[0]
    ph = IC[1]
    a = IC[2]
    al = IC[3]

    dR_dr = ph
    da_dr = a*((2*l+1)*r/2*(om**2*a**2*R**2/al**2+ph**2+l*(l+1)*a**2*R**2/r**2+m**2*a**2*R**2)-(a**2-1)/(2*r))
    dal_dr = al*(da_dr/a-l*(l+1)*(2*l+1)*a**2*R**2/r-(2*l+1)*m**2*a**2*r*R**2+(a**2-1)/r)
    dph_dr = -2*ph/r-dal_dr*ph/al+da_dr*ph/a-om**2*a**2*R/al**2+l*(l+1)*a**2*R/r**2+m**2*a**2*R
    
    return [dR_dr,dph_dr,da_dr,dal_dr]

def init(u,p,r):
    if p==0:
        return np.array([1.,r,1.,u])
    else:
        return np.array([r**p,l*r**(p-1),1,u])


l = 0.
m = 1.
ep = 0.2
n_om = 100.
omega = np.linspace(m-ep,m+ep,n_om)
#r = np.linspace(0.001, 100, 1000)
r_end = 100.
r_start = 0.001
r_step = 0.01
r_interval = np.arange(r_start,r_end,r_step)

niter = 1000
tol = 0.0001
ustep = 0.01

for j in range(len(omega)):
    print('trying with $omega =$',omega[j])
    p = (l,m,omega[j])
    u = 0.001
    for i in range(niter):
        u += ustep
        ini = init(u,p[0] ,r_start)
        o = ode(f).set_integrator('vode', method='bdf', order=15, nsteps=3000)
        o.set_initial_value(ini,r_start).set_f_params(p)
        rs = []
        ys = []
        while o.successful() and o.t < r_end:
            o.integrate(o.t+r_step)
            rs.append(o.t)
            ys.append(o.y)
        r = np.vstack(rs)
        Y = np.vstack(ys).T
        if abs(Y[2,len(r)-1]-1/Y[3,len(r)-1]) < tol:
            break
    if abs(Y[0,len(r)-1]) < tol and abs(Y[2,len(r)-1]-1/Y[3,len(r)-1]) < tol:
        print(j,'times iterations in omega')
        print(i,'times iterations in alpha')
        print("R'(inf)) = ", Y[len(Y)-1,0])        
        print("alpha(0)) = ", Y[0,3])
        print("\omega",omega[j])
        plt.subplot(2,1,1)
        plt.plot(r,Y[0,:],'r',label = '$R$')
        plt.plot(r,Y[1,:],'b',label = '$d R /dr$')
        plt.xlim([0,10])
        plt.legend()
        plt.subplot(2,1,2)
        plt.plot(r,Y[2,:],'r',label = 'a')
        plt.plot(r,Y[3,:],'b', label = '$alpha$')
        plt.xlim([0,10])
        plt.legend()
        plt.show()
        break

