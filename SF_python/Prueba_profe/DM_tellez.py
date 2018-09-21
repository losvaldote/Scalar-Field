import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class DarkM:
    def __init__(self, mquin=1, name='Alberto'):

        self.name   = name
        self.cte    = 3/2.0


        self.NP = 100000 #000
        self.Ni = np.log(1.0e-6)
        self.Nf = np.log(1.0e-0)
        self.d  = (self.Nf- self.Ni)/self.NP
        self.t = [i for i in np.linspace(self.Ni, self.Nf, self.NP)]


    def something(self):
        print self.name


    def rk4(self, func, y_0, t, args={}):
        # Inicia el arreglo de las aproximaciones
        #y = np.zeros([4, len(y_0)])
        y = np.zeros([len(t), len(y_0)])
        y[0] = y_0
        for i, t_i in enumerate(t[:-1]):
        #for i, t_i in enumerate(t[:3]):

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




    """ DE 
    def solver(self):
        y0       = np.array([np.sqrt(0.00), np.sqrt(0.729), np.sqrt(0.0001), 1.0E-15, np.sqrt(0.27)])
        y_result = self.rk4(self.RHS, y0, self.t)
        return y_result


    def RHS(self, t, y, args):
        gamma = 1
        x0, x1, x2, x3, x4= y
        Pe = 2.0*x0*x0 +4.0*x2*x2/3.0 + gamma*x4*x4
        return np.array([-3*x0 - x1*x3 + self.cte*Pe*x0,
                        x0*x3 + self.cte*Pe*x1,
                        self.cte*(Pe-(4.0/3.0))*x2,
                        self.cte*Pe,#*x3,
                        self.cte*(Pe-gamma)*x4])
    


    """ #DM 
    def solver(self):
        y0       = np.array([np.sqrt(0.001), np.sqrt(0.27), np.sqrt(0.0001), 1.0E3, np.sqrt(0.729)])
        y_result = self.rk4(self.RHS, y0, self.t)
        return y_result


    def RHS(self, t, y, args):
        gamma = 0
        x0, x1, x2, x3, x4= y
        Pe = 2.0*x0*x0 +4.0*x2*x2/3.0 + gamma*x4*x4
        return np.array([-3*x0 - x1*x3 + self.cte*Pe*x0,
                        x0*x3 + self.cte*Pe*x1,
                        self.cte*(Pe-(4.0/3.0))*x2,
                        self.cte*Pe,
                        self.cte*(Pe-gamma)*x4])

	






    def plot(self):
        x, u, z, s, l = self.solver().T
        fig = plt.figure(figsize=(9,10))
        ax3 = fig.add_subplot(111)

        i =0
        tiempo = []
        phi = []
        dphi = []
        lam  = []
        ese   = []
        rad   = []
        for t, p, d, r2, l2 in zip(self.t, x, u, z, l):
            if i%100 ==0:
               tiempo.append(t)
               phi.append(p**2)
               dphi.append(d**2)
               lam.append(l2**2 )
               rad.append(r2**2)
               #ese.append(c)
            i+=1
        #ax3.plot(tiempo, np.array(rad), 'blue', label ='phi')
        #ax3.plot(tiempo, np.array(lam), 'red', label ='dphi')
        #ax3.semilogx(tiempo, np.array(lam), 'green', label ='lam')
        #ax3.semilogx(tiempo, np.array(rad), 'yellow', label ='rad')
        ax3.plot(tiempo, np.array(dphi) + np.array(phi), 'black')
	#ax3.semilogx(tiempo, np.array(dphi) + np.array(phi), 'black')
        #ax3.plot(tiempo, (np.array(dphi) - np.array(phi))/(np.array(dphi) + np.array(phi)), 'black')
        plt.xlabel('$a$', fontsize=20)
        plt.ylabel('$\Omega(a)$', fontsize=20)
        plt.legend(loc = 'best')
        plt.savefig('Omega_dark_matter.pdf')
        plt.show()








if __name__ == '__main__':
    DM = DarkM(mquin= 2.)
print DM.plot()
