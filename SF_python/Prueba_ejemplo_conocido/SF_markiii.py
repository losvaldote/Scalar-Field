import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from math import sin, pi, sqrt



class DarkM:
    def __init__(self, mquin=1, name='Alberto'):

        self.name   = name
        self.cte    = 3./2.0
	self.c = 1
	self.Theta = pi/4
	self.freq = sqrt(self.c)      # frequency of oscillations when Theta is small
	self.period = 2*pi/self.freq  # the period of the oscillations
	self.T = 10*self.period       # final time
	self.N_per_period = 20   # resolution of one period
	self.N = self.N_per_period*self.period
	self.time_points = np.linspace(0, self.T, self.N+1)


        self.NP = 100000
        self.Ni = 0
        self.Nf = 10*self.period 
        self.d  = (self.Nf- self.Ni)/self.NP
        self.t = [np.exp(self.Ni+self.d*i) for i in np.arange(self.NP)]


	
	#self.Nf = self.T

    def something(self):
        print self.name


    def rk4(self, func, y_0, t, args={}):
        # Inicia el arreglo de las aproximaciones
        y = np.zeros([4, len(y_0)])
        y[0] = y_0
        for i, t_i in enumerate(t[:3]): #Revisar el contenido del enumerate

            h = self.d #t[i+1] - t_i
            k_1 = func(t_i, y[i], args)
            k_2 = func(t_i+h/2., y[i]+h/2.*k_1, args)
            k_3 = func(t_i+h/2., y[i]+h/2.*k_2, args)
            k_4 = func(t_i+h, y[i]+h*k_3, args)

            y[i+1] = y[i] + h/6.*(k_1 + 2.*k_2 + 2.*k_3 + k_4) # RK4 step

        return y


#Adams-Bashforth 4/Moulton 4 Step Predictor/Corrector
    def ABM3(self, func, y_0, t, args={}):
	y = np.zeros([len(t), len(y_0)])
	#Se calcularan los primeros pasos con rk4
        y[0:4] = self.rk4(func,y_0, t)
        k_1 = func(t[2], y[2], args)
	k_2 = func(t[1], y[1], args)
	k_3 = func(t[0], y[0], args)
	for i in range(3,self.NP-1):
             h = self.d
	     k_4 = k_3
	     k_3 = k_2
	     k_2 = k_1
             k_1 = func(t[i], y[i], args)
	     #Adams Bashforth predictor
	     y[i+1] = y[i] + h*(55.*k_1 - 59.*k_2 + 37.*k_3 - 9.*k_4)/24.
	     k_0 = func(t[i+1],y[i+1], args)
	     #Adams Moulton corrector
	     y[i+1] = y[i] + h*(9.*k_0 + 19.*k_1 - 5.*k_2 + k_3)/24.
        return y 


    def solver(self):
        y0       = np.array([self.Theta,0])
        y_result = self.ABM3(self.RHS, y0, self.t)
        return y_result


    def RHS(self, t, y, args):
        theta, omega = y
        return [omega, -self.c*sin(theta)]


    def plot(self):
        x, z= self.solver().T
        fig = plt.figure(figsize=(9,10))
        ax3 = fig.add_subplot(111)

        i =0
        tiempo = []
        phi = []
        dphi = []
        for t, p, d  in zip(self.t, x, u):
            if i%200 ==0:
               tiempo.append(t)
               phi.append(p)
               dphi.append(d)
            i+=1
        #ax3.semilogx(tiempo, np.array(phi)*np.array(phi), 'blue')   #cinetica
        #ax3.semilogx(tiempo, np.array(dphi)*np.array(dphi), 'red') #potencial
	ax3.semilogx(tiempo, (np.array(phi)*np.array(phi) - np.array(dphi)*np.array(dphi))/(np.array(phi)*np.array(phi) + np.array(dphi)*np.array(dphi)), 'red') #(xx-uu)/(xx+uu)
        plt.xlabel('$\ln a$', fontsize=20)
        #plt.ylabel('$\Omega(a)$', fontsize=20) original
	#plt.ylabel('$x(a)$', fontsize=20) #cinetica
	plt.ylabel('$u(a)$', fontsize=20) #potencial
	#plt.ylabel('$\Omega(a)$', fontsize=20) cosa rara
        plt.savefig('Omega_dark_matter_cosa_rara.pdf')
        plt.show()








if __name__ == '__main__':
    DM = DarkM(mquin= 2.)
    print DM.plot()

