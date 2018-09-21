import numpy as np
import odespy
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class DarkM:
    def __init__(self, mquin=1, name='Alberto'):

        self.name   = name
        self.cte    = 3./2.0


        self.NP = 100000
        self.Ni = np.log(1.0e-0)
        self.Nf = np.log(1.0e-6)
        self.d  = (self.Nf- self.Ni)/self.NP
        self.t = [np.exp(self.Ni+self.d*i) for i in np.arange(self.NP)] #Es el arreglo para el tiempo


    def something(self):
        print self.name


    def RHS(self, t, y):
        x0, x1, x2, x3, x4 = y
        Pe = 2.0*x0*x0+4.0*x2*x2/3.0
        return np.array([-3.*x0-x1*x4+self.cte*Pe*x0, x0*x4+self.cte*Pe*x1, self.cte*(Pe-(4.0/3.0))*x2,
                self.cte*Pe*x3, self.cte*Pe])


    def solution(self):
        y0       = np.array([np.sqrt(0.2295), np.sqrt(0.00043), np.sqrt(0.000043), np.sqrt(0.73), 1.0E3])
 	solver = odespy.AdamsBashMoulton3(self.RHS)  
	solver.set_initial_condition(y0)     
	y_result = solver.solve(self.t)
        return y_result

    def plot(self):
        x, u, z, l, s = self.solution().T
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
	ax3.semilogx(tiempo, (np.array(phi)*np.array(phi) + np.array(dphi)*np.array(dphi)), 'red') #(xx-uu)/(xx+uu)
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

