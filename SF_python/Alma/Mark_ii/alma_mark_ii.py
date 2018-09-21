from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from IPython.display import Image
import math
from math import *
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import ode

#Este es el metodo de shooting para hallar las condiciones iniciales en el pasado.
def Secante_1D(A0,A1,f,m,a_i,lna_step,tol=0.0001,N=100):
    print "Trying with our first guesses for A:",A0,A1
    f_x0 = f(A0,m,a_i,lna_step)
    f_x1 = f(A1,m,a_i,lna_step)
    print "First residuos:", f_x0, f_x1
    iteration_counter = 0
    while abs(f_x1) > tol and iteration_counter < N:
        try:
            denominator = float(f_x1 - f_x0)/(A1 - A0)
            A = A1 - float(f_x1)/denominator
            print "Our new try for A:",A
        except ZeroDivisionError:
            print "Error! - denominator zero for A = ", A
            sys.exit(1)     # Abort with error
        A0 = A1
        A1 = A
        f_x0 = f_x1
        f_x1 = f(A1,m,a_i,lna_step)
        print "New residuo:",f_x1
        iteration_counter += 1
    # Here, either a solution is found, or too many iterations
    if abs(f_x1) > tol:
        iteration_counter = -1
    return A, iteration_counter

def lcdm(IC,lna):
    v,w,x,z = IC
    om_b, om_dm, om_r, om_l = [0.,0,1./3.,-1.]  

    om_t = om_b*v+om_r*x+om_l*z+om_dm*w

    return [3.0*(om_t-om_b)*v,
    		3.0*(om_t-om_dm)*w,
    		3.0*(om_t-om_r)*x,
    		3.0*(om_t-om_l)*z]


def ic_lcdm():
	#return [4.559568563844155e-06, 2.6325639011579477e-05, 0.9999691147924246, 9.276384105305625e-29]
#	return[9.803072412264933e-06, 5.4757329144085316e-05, 0.9999354395984436, 1.855276821061125e-31]
#	return [4.559568563844155e-06, 0.9999691147924246, 9.276384105305625e-29]
	return[9.803072412264933e-06, 0.9999354395984436, 1.855276821061125e-31]


#La masa m del campo escalar estara medida en terminos del factor de escala hoy en dia
def ic_lsfdm(A,m,a_i):    
    theta_i = 2.0*m*a_i**2./(5.*(1.0e-4)**(1./2.))
    y_1_i = 5.0*theta_i
    w_i = a_i*(0.22994/0.00004)*((4.*theta_i/math.pi**2.)*((9.+math.pi**2./4.)/(9.+theta_i**2.)))**(3./4.)
    v_i, x_i, z_i = ic_lcdm()
    x_i = x_i+(1.-x_i-v_i-w_i-z_i)
    #d = 1.
    #z_i = 1.-v_i-A*w_i-d*x_i
    return [y_1_i,theta_i,v_i,A*w_i,x_i,z_i]
    
def background_system_real(lna,ic):

    y_1, theta, v, w, x, z = ic
    if theta> 1000:
        om_b, om_sf, om_r, om_l = [0.,0,1./3.,-1.]  
    else:
        om_b, om_sf, om_r, om_l = [0.,-math.cos(theta),1./3.,-1.]  
    om_t = om_b*v+om_r*x+om_l*z+om_sf*w

    return [3.0*(1.+om_t)*y_1/2.,
    		-3.0*math.sin(theta)+y_1,
    		3.0*(om_t-om_b)*v,
    		3.0*(om_t-om_sf)*w,
    		3.0*(om_t-om_r)*x,
    		3.0*(om_t-om_l)*z]

def fs(A,m,a_i,lna_step):
    IC = ic_lsfdm(A,m,a_i)
    o = ode(background_system_real).set_integrator('dopri5',nsteps = 3000)
    o.set_initial_value(IC,math.log(a_i))
    rs = []
    ys = []
    while o.successful() and o.t < 0:
        o.integrate(o.t+lna_step)
        rs.append(o.t)
        ys.append(o.y)
    lna = np.vstack(rs)
    Y = np.vstack(ys).T	
    #print "residuos", Y[2,-1]-0.04, Y[3,-1]-0.22994, Y[4,-1]-1.0e-4
    return np.array([Y[3,-1]-0.22994])


a_i = 1.0e-9
lna = np.linspace(0., math.log(a_i), 500)
#IC_lcdm = [0.04,0.22994,0.0001,1-0.04-0.22994-0.0001]
#y = odeint(lcdm,IC_lcdm,lna)

#plt.figure(figsize=(10,5))
#plt.plot(lna, y[:,0],label="Baryons")
#plt.plot(lna, y[:,1],label="CDM")
#plt.plot(lna, y[:,2],label="Photons")
#plt.plot(lna, y[:,3],label="DE")
#plt.legend()
#plt.title('LCDM')
#plt.show()

#m = 10.0e19
m1 = 10.0e13
lna_step = 0.0001
A0 = 1.
A1 = 0.001
A_real, n_iter = Secante_1D(A0,A1,fs,m1,a_i,lna_step,tol=0.00001,N=200)


IC = ic_lsfdm(A_real,m1,a_i)

print "Trying with our best initial conditions:",IC
#Y = odeint(background_system_real,IC,lna,mxstep=1000000)
#o = ode(background_system_real).set_integrator('vode', method='bdf', order=15, nsteps=5000)
#o = ode(background_system_real).set_integrator('vode', method='adams', order=12, nsteps=3000)
#o = ode(background_system_real).set_integrator('lsoda',nsteps = 5000)
o = ode(background_system_real).set_integrator('dopri5',nsteps = 4000)
#o = ode(background_system_real).set_integrator('dop853',nsteps = 3000)

o.set_initial_value(IC,math.log(a_i))
rs = []
ys = []
while o.successful() and o.t < 0:
    o.integrate(o.t+lna_step)
    rs.append(o.t)
    ys.append(o.y)
lna = np.vstack(rs)
Y = np.vstack(ys).T
om_phi = lna*0.
for i in range(len(lna)): 
	om_phi[i] = -math.cos(Y[1,i])

np.savetxt("data_m_e13",Y)

#Probaremos con masa 2
m2 = 10.0e14
A0 = 1.
A1 = 0.001
A_real2, n_iter2 = Secante_1D(A0,A1,fs,m2,a_i,lna_step,tol=0.00001,N=200)
IC = ic_lsfdm(A_real2,m2,a_i)
print "Trying with our best initial conditions:",IC
o = ode(background_system_real).set_integrator('dopri5',nsteps = 4000)
o.set_initial_value(IC,math.log(a_i))
rs2 = []
ys2 = []
while o.successful() and o.t < 0:
    o.integrate(o.t+lna_step)
    rs2.append(o.t)
    ys2.append(o.y)
lna = np.vstack(rs2)
Y2 = np.vstack(ys2).T
om_phi2 = lna*0.
for i in range(len(lna)): 
	om_phi2[i] = -math.cos(Y2[1,i])

np.savetxt("data_m_e13_2",Y2)

#Probaremos con masa 3
m3 = 10.0e8 #m=10^-25
A0 = 1.
A1 = 0.001
A_real3, n_iter3 = Secante_1D(A0,A1,fs,m3,a_i,lna_step,tol=0.00001,N=200)
IC = ic_lsfdm(A_real3,m3,a_i)
print "Trying with our best initial conditions:",IC
o = ode(background_system_real).set_integrator('dopri5',nsteps = 4000)
o.set_initial_value(IC,math.log(a_i))
rs3 = []
ys3 = []
while o.successful() and o.t < 0:
    o.integrate(o.t+lna_step)
    rs3.append(o.t)
    ys3.append(o.y)
lna = np.vstack(rs3)
Y3 = np.vstack(ys3).T
om_phi3 = lna*0.
for i in range(len(lna)): 
	om_phi3[i] = -math.cos(Y3[1,i])

np.savetxt("data_m_e13_3",Y3)

#plt.figure(figsize=(10,5))
#plt.plot(lna, Y[0,:],label="$\theta$")
#plt.plot(lna, Y[1,:],label="$y_1$")
#plt.title('SFDM')
#plt.legend()
#plt.show()


#plt.figure(figsize=(10,5))
#plt.plot(lna, om_phi,label="Ec. Est.")
#plt.ylim(-1.1,1.1)
#plt.title('SFDM')
#plt.legend()
#plt.show()

#plt.figure(figsize=(10,5))
#plt.plot(lna, Y[2,:],label="Baryons")
#plt.plot(lna, Y[3,:],label="SFDM")
#plt.plot(lna, Y[4,:],label="Photons")
#plt.plot(lna, Y[5,:],label="Dark energy")
#plt.ylim(-0.1,1.1)

#plt.title('SFDM')
#plt.legend()
#plt.show()

#plt.figure(figsize=(10,5))
#plt.plot(lna,  Y[2,:]+ Y[3,:]+ Y[4,:]+ Y[5,:],label="Ec. Est.")
#plt.ylim(.95,1.05)
#plt.title('SFDM')
#plt.legend()
#plt.show()



#En esta parte se tratara de graficar con el factor de escala en lugar del logaritmo.
fact = np.exp(lna)

fig = plt.figure(figsize=(10,5))
ax3 = fig.add_subplot(111)

#Masa 1
ax3.semilogx(fact, Y[2,:], ls = 'dashed', label="$\Omega_b$")
ax3.semilogx(fact, Y[3,:], ls = 'dashed', label="$\Omega_{dm}$")
ax3.semilogx(fact, Y[4,:], ls = 'dashed', label="$\Omega_{\gamma}$")
ax3.semilogx(fact, Y[5,:], ls = 'dashed', label="$\Omega_{\Lambda}$")

#Masa 2
ax3.semilogx(fact, Y2[2,:], ls = 'dotted', label="$\Omega_b2$")
ax3.semilogx(fact, Y2[3,:], ls = 'dotted', label="$\Omega_{dm}2$")
ax3.semilogx(fact, Y2[4,:], ls = 'dotted', label="$\Omega_{\gamma}2$")
ax3.semilogx(fact, Y2[5,:], ls = 'dotted', label="$\Omega_{\Lambda}2$")

#Masa 3
ax3.semilogx(fact, Y3[2,:],label="$\Omega_b3$")
ax3.semilogx(fact, Y3[3,:],label="$\Omega_{dm}3$")
ax3.semilogx(fact, Y3[4,:],label="$\Omega_{\gamma}3$")
ax3.semilogx(fact, Y3[5,:],label="$\Omega_{\Lambda}3$")

#Leyendas
plt.xlabel('$a$', fontsize=20)
plt.ylabel('$\Omega(a)$', fontsize=20) #original
#plt.legend(('$\Omega_{dm}$', '$\Omega_{\gamma}$', '$\Omega_{\Lambda}$', '$\Omega_b$' ))
plt.xlim(1.0e-7,1)
plt.ylim(-0.1,1.1)

plt.title('SFDM')
plt.legend()
plt.show()


