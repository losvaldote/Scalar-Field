import numpy as np
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
        self.t = [np.exp(self.Ni+self.d*i) for i in np.arange(self.NP)]


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
    def ABM4(self, func, y_0, t, args={}):
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


    def solver(self, masa):
        y0       = np.array([np.sqrt(0.2295), np.sqrt(0.00043), np.sqrt(0.000043), np.sqrt(0.73), masa])
        y_result = self.ABM4(self.RHS, y0, self.t)
        return y_result


    def RHS(self, t, y, args={}):
        x0, x1, x2, x3, x4 = y
        Pe = 2.0*x0*x0+4.0*x2*x2/3.0
        return np.array([-3.*x0-x1*x4+self.cte*Pe*x0, x0*x4+self.cte*Pe*x1, self.cte*(Pe-(4.0/3.0))*x2,
                self.cte*Pe*x3, self.cte*Pe])


    def plot(self,masa):
        x, u, z, l, s = self.solver(masa).T
        fig = plt.figure(figsize=(9,10))
	plt.ion()
        ax3 = fig.add_subplot(111)

        i =0
       	tiempo = []
       	phi = []
       	dphi = []
       	for t, p, d  in zip(self.t, x, u):
       	    if i%100 ==0:
       	       tiempo.append(t)
       	       phi.append(p)
       	       dphi.append(d)
       	    i+=1
       	#ax3.semilogx(tiempo, np.array(phi)*np.array(phi), 'blue')   #cinetica
       	#ax3.semilogx(tiempo, np.array(dphi)*np.array(dphi), 'red') #potencial
	ax3.semilogx(tiempo, (np.array(phi)*np.array(phi) + np.array(dphi)*np.array(dphi)), 'red') #(xx+uu)
       	plt.xlabel('$\ a$', fontsize=20)
       	#plt.ylabel('$\Omega(a)$', fontsize=20) original
	#plt.ylabel('$x(a)$', fontsize=20) #cinetica
	#plt.ylabel('$u(a)$', fontsize=20) #potencial
	plt.ylabel('$u^2(a) + x^2(a)$', fontsize=20) #cosa rara
	plt.title("Figura con masa " + str(masa))
       	plt.savefig('masa_' + str(masa) + '.pdf') #Convierte el valor de la masa a cadena para poder guardar varias graficas.
       	plt.show()
	plt.close()
	
	

   





if __name__ == '__main__':
    DM = DarkM(mquin= 2.)
    for ms in (0,1E3,1E2):
		print DM.plot(ms)

