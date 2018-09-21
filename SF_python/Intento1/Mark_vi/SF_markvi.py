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


    def solver(self,masa):
        y0       = np.array([np.sqrt(0.2295), np.sqrt(0.00043), np.sqrt(0.000043), np.sqrt(0.73), masa, np.sqrt(0.04) ])
        y_result = self.ABM4(self.RHS, y0, self.t)
	#y_result = self.rk4(self.RHS, y0, self.t)
        return y_result


    def RHS(self, t, y, args={}):
        x0, x1, x2, x3, x4, x5 = y
        Pe = 2.0*x0*x0+4.0*x2*x2/3.0+x5*x5
        return np.array([-3.*x0-x1*x4+self.cte*Pe*x0, x0*x4+self.cte*Pe*x1, self.cte*(Pe-(4.0/3.0))*x2,
                self.cte*Pe*x3, self.cte*Pe, self.cte*(Pe - 1.0)*x5])


    def plot(self):
        #Primer valor de la masa.
        m1 = 10
        x, u, z, l, s, b = self.solver(m1).T
        fig = plt.figure(figsize=(9,10))
        ax3 = fig.add_subplot(111)

        i =0
        tiempo = []
        phi = []     #velocidad
        dphi = []    #potencial
	rad = []     #radiacion
	de = []	     #dark energy
        bar = []       #bariones
        for t, p, d, aux1, aux2, aux3  in zip(self.t, x, u, z, l, b):
            if i%200 ==0:
               tiempo.append(t)
               phi.append(p)
               dphi.append(d)
	       rad.append(aux1)
               de.append(aux2)
               bar.append(aux3)
            i+=1
	#Segundo valor de la masa.
	m2 = 100
        x2, u2, z2, l2, s2, b2 = self.solver(m2).T
        #fig = plt.figure(figsize=(9,10))
        #ax3 = fig.add_subplot(111)

        i =0
        tiempo2 = []
        phi2 = []     #velocidad
        dphi2 = []    #potencial
	rad2 = []     #radiacion
	de2 = []	     #dark energy
        bar2 = []       #bariones
        for t2, p2, d2, aux21, aux22, aux23  in zip(self.t, x2, u2, z2, l2, b2):
            if i%200 ==0:
               tiempo2.append(t2)
               phi2.append(p2)
               dphi2.append(d2)
	       rad2.append(aux21)
               de2.append(aux22)
               bar2.append(aux23)
            i+=1

	#Tercer valor de la masa.
	m3 = 1000
        x3, u3, z3, l3, s3, b3 = self.solver(m3).T
        #fig = plt.figure(figsize=(9,10))
        #ax3 = fig.add_subplot(111)

        i =0
        tiempo3 = []
        phi3 = []     #velocidad
        dphi3 = []    #potencial
	rad3 = []     #radiacion
	de3 = []	     #dark energy
        bar3 = []       #bariones
        for t3, p3, d3, aux31, aux32, aux33  in zip(self.t, x3, u3, z3, l3, b3):
            if i%200 ==0:
               tiempo3.append(t3)
               phi3.append(p3)
               dphi3.append(d3)
	       rad3.append(aux31)
               de3.append(aux32)
               bar3.append(aux33)
            i+=1


	#Primer grafica todos
        #ax3.semilogx(tiempo, np.array(phi)*np.array(phi), 'blue')   #cinetica
        #ax3.semilogx(tiempo, np.array(dphi)*np.array(dphi), 'red') #potencial
        #ax3.semilogx(tiempo, (np.array(phi)*np.array(phi) + np.array(dphi)*np.array(dphi)) + np.array(rad)*np.array(rad) + np.array(de)*np.array(de) + np.array(bar)*np.array(bar), 'black') #ecuacion de estado
	#ax3.semilogx(tiempo, (np.array(phi)*np.array(phi) + np.array(dphi)*np.array(dphi)),  'black', ls = 'dashed') #(xx+uu)	
	#ax3.semilogx(tiempo, np.array(rad)*np.array(rad), 'blue', ls = 'dashed') #radiacion
	#ax3.semilogx(tiempo, np.array(de)*np.array(de), 'red', ls = 'dashed') #dark energy
	#ax3.semilogx(tiempo, np.array(bar)*np.array(bar), 'green', ls = 'dashed')#bariones
	
	#Segunda grafica todos
	#ax3.semilogx(tiempo2, (np.array(phi2)*np.array(phi2) + np.array(dphi2)*np.array(dphi2)), ls = 'dotted') #(xx+uu)	
	#ax3.semilogx(tiempo2, np.array(rad2)*np.array(rad2), ls = 'dotted') #radiacion
	#ax3.semilogx(tiempo2, np.array(de2)*np.array(de2), ls = 'dotted') #dark energy
	#ax3.semilogx(tiempo2, np.array(bar2)*np.array(bar2), ls = 'dotted')#bariones

	#Tercer grafica todos
        #ax3.semilogx(tiempo3, (np.array(phi3)*np.array(phi3) + np.array(dphi3)*np.array(dphi3))) #(xx+uu)	
	#ax3.semilogx(tiempo2, np.array(rad3)*np.array(rad3)) #radiacion
	#ax3.semilogx(tiempo2, np.array(de3)*np.array(de3)) #dark energy
	#ax3.semilogx(tiempo2, np.array(bar3)*np.array(bar3))#bariones
	
	#Leyendas todos
        #plt.xlabel('$a$', fontsize=20)
        #plt.ylabel('$\Omega(a)$', fontsize=20) #original
	#plt.ylabel('$x^2(a)$', fontsize=20) #cinetica
	#plt.ylabel('$u^2(a)$', fontsize=20) #potencial
	#plt.ylabel('$x^2(a) + u^2(a)$', fontsize=20) #Suma
	#plt.legend(('$\Omega_{dm}1$', '$\Omega_{\gamma}1$', '$\Omega_{\Lambda}1$', '$\Omega_b1$', '$\Omega_{dm}2$', '$\Omega_{\gamma}2$', '$\Omega_{\Lambda}2$', '$\Omega_b2$', '$\Omega_{dm}3$', '$\Omega_{\gamma}3$', '$\Omega_{\Lambda}3$', '$\Omega_b3$'))
        
	#Solo materia oscura
 	#ax3.semilogx(tiempo, (np.array(phi)*np.array(phi) + np.array(dphi)*np.array(dphi))) #(xx+uu)
	#ax3.semilogx(tiempo2, (np.array(phi2)*np.array(phi2) + np.array(dphi2)*np.array(dphi2))) #(xx+uu)
	#ax3.semilogx(tiempo3, (np.array(phi3)*np.array(phi3) + np.array(dphi3)*np.array(dphi3))) #(xx+uu)	

	#Ecuacion de estado
        ax3.semilogx(tiempo, (np.array(phi)*np.array(phi) - np.array(dphi)*np.array(dphi))/(np.array(phi)*np.array(phi) + np.array(dphi)*np.array(dphi)))
	ax3.semilogx(tiempo2, (np.array(phi2)*np.array(phi2) - np.array(dphi2)*np.array(dphi2))/(np.array(phi2)*np.array(phi2) + np.array(dphi2)*np.array(dphi2)))
	ax3.semilogx(tiempo3, (np.array(phi3)*np.array(phi3) - np.array(dphi3)*np.array(dphi3))/(np.array(phi3)*np.array(phi3) + np.array(dphi3)*np.array(dphi3)))

	plt.xlabel('$a$', fontsize=20)
        #plt.ylabel('$\Omega(a)_{DM}$', fontsize=20) #original
	#plt.legend(('$\Omega_{dm}1$', '$\Omega_{dm}2$', '$\Omega_{dm}3$'))
	plt.ylabel('$\omega(a)$', fontsize=20) #ecuacion de estado
	plt.legend(('$\omega 1$', '$\omega 2$', '$\omega 3$'))#ecuacion de estado

	plt.savefig('ecuacion_de_estado_mark_vi.pdf')
        plt.show()












if __name__ == '__main__':
    DM = DarkM(mquin= 2.)
    print DM.plot()

