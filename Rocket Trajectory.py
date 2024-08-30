import numpy as np
import math as mt
import matplotlib.pyplot as plt
dt = 0.1
q = np.arange(0, 100.1, dt)
t = range(100)
ux = np.zeros(101,dtype = float)
uy = np.zeros(101,dtype = float)
u = np.zeros(101,dtype = float)
un = np.zeros(101,dtype = float)
x = np.zeros(101,dtype = float)
y = np.zeros(101,dtype = float)
theta = np.zeros(101,dtype = float)
m = np.zeros(101,dtype = float)
g = np.zeros(101,dtype = float)
du = np.zeros(101,dtype = float)
dun = np.zeros(101,dtype = float)
dur = np.zeros(101,dtype = float)
dtheta = np.zeros(101,dtype = float)
dux = np.zeros(101,dtype = float)
duy = np.zeros(101,dtype = float)
dx = np.zeros(101,dtype = float)
dy = np.zeros(101,dtype = float)
D = np.zeros(101,dtype = float)
Cd = np.zeros(101,dtype = float)
rho = np.zeros(101,dtype = float)
T = np.zeros(101,dtype = float)
a = np.zeros(101,dtype = float)
M = np.zeros(101,dtype = float)
# Initial conditions
T[0] = 288.15
rho[0] = 1.225
j = mt.pi/180
theta[0] = 5*j
m[0] = 15000
g[0] = 9.81
ueq = 3000
A = 1
for i in t:
  m[i+1] = m[i] - 20
  du[i+1] = ((200*3000)/m[i] - (D[i]/m[i]) - g[i]*mt.cos(theta[i]))
  dun[i+1] = (g[i]*np.sin(theta[i]))*dt
  dur[i+1] = np.sqrt(du[i+1]**2 + dun[i+1]**2)
  dtheta[i+1] = np.arctan(dun[i+1]/du[i+1])
  theta[i+1] = theta[i] + dtheta[i+1]
  dux[i+1] = dur[i+1]*np.sin(theta[i+1])
  duy[i+1] = dur[i+1]*np.cos(theta[i+1])
  ux[i+1] = ux[i] + dux[i+1]
  uy[i+1] = uy[i] + duy[i+1]
  dx[i+1] = ux[i+1]*dt
  dy[i+1] = uy[i+1]*dt
  x[i+1] = x[i] + dx[i+1]
  y[i+1] = y[i] + dy[i+1]
  u[i+1] = np.sqrt(ux[i+1]**2 + uy[i+1]**2)
  g[i+1] = 9.81*(1 - (2*y[i+1])/6371000)
  rho[i+1] = 1.125
  if y[i+1]<=11000:
    T[i+1] = 288.15 - (0.0065*y[i+1])
  if 25000>=y[i+1]>11000:
    T[i+1] = 216.65
  if y[i+1]>25000:
    T[i+1] = 141.79 + (0.00299*y[i+1])
  a[i+1] = np.sqrt(401.8*T[i+1])
  M[i+1] = u[i+1]/a[i+1]
  if 0<=M[i+1]<0.7:
    Cd[i+1] = 0.4
  if 0.7<=M[i+1]<1:
    Cd[i+1] = 0.5
  if 1<=M[i+1]<1.5:
    Cd[i+1] = 0.7
  if 1.5<=M[i+1]<2:
    Cd[i+1] = 0.5
  if 2<=M[i+1]<2.5:
    Cd[i+1] = 0.4
  if M[i+1]>3:
    Cd[i+1] = 0.2
  D[i+1] = 0.5*Cd[i+1]*rho[i+1]*(u[i+1]**2)*A
plt.plot(x,y)
plt.show()