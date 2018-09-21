from numpy import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def f(t,y):
	return (t-3.2)*y + 8*t*exp((t-3.2)**2/2)*cos(4*t**2)

def dfy(t,y):
	return t-3.2

#Runge-Kutta "Classic" Order 4 method
def RK4(t0,tn,n,y0):
	h = abs(tn-t0)/n
	t = linspace(t0,tn,n+1)
	y = zeros(n+1)
	y[0] = y0
	for i in range(0,n):
		K1 = f(t[i],y[i])
		K2 = f(t[i]+h/2,y[i]+K1*h/2)
		K3 = f(t[i]+h/2,y[i]+K2*h/2)
		K4 = f(t[i]+h,y[i]+K3*h)
		y[i+1] = y[i] + h*(K1+2*K2+2*K3+K4)/6
	return y

#Adams-Bashforth 3/Moulton 4 Step Predictor/Corrector
def PreCorr3(t0,tn,n,y0):
	h = abs(tn-t0)/n
	t = linspace(t0,tn,n+1)
	y = zeros(n+1)
	#Calculate initial steps with Runge-Kutta 4
	y[0:3] = RK4(t0,t0+2*h,2,y0)
	K1 = f(t[1],y[1])
	K2 = f(t[0],y[0])
	for i in range(2,n):
		K3 = K2
		K2 = K1
		K1 = f(t[i],y[i])
		#Adams-Bashforth Predictor
		y[i+1] = y[i] + h*(23*K1-16*K2+5*K3)/12
		K0 = f(t[i+1],y[i+1])
		#Adams-Moulton Corrector
		y[i+1] = y[i] + h*(9*K0+19*K1-5*K2+K3)/24
	return y
#Script to produce graphs
fg =1
n = 300
t0 = 0
tn = 6
y0 = .75
t = linspace(t0,tn,n+1)
ypc = PreCorr3(t0,tn,n,y0)
plt.plot(t,ypc,label='Predictorcorrector')
plt.show()

