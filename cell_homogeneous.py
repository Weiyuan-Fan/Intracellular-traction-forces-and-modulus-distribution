#!/usr/bin/env python
# coding: utf-8

# In[1]:


# length unit is micrometer
# force unit is micronewton
# pressure unit is megapascal
# thickness is 5 micrometer
from mshr import *
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
mesh = Mesh("cell6.xml")
#plot(mesh)
#plt.gca().invert_yaxis() #flip y axis so that it looks the same as the real image

# read the positions and values of the loads from data.txt
x_p = np.loadtxt("celldata6.txt")[:, 0]
y_p = np.loadtxt("celldata6.txt")[:, 1]
f_x = np.loadtxt("celldata6.txt")[:, 2]/5
f_y = np.loadtxt("celldata6.txt")[:, 3]/5 # put a minus sign because y axis is flipped

# patchload 
class PatchLoad(UserExpression):
    def __init__(self, pt, vl,r,**kwargs):
        super().__init__(**kwargs)
        self.point = pt
        self.value = vl
        self.r     = r   
    def eval(self, values, x):
        if (x[0] - self.point[0])**2+(x[1] - self.point[1])**2 <= (self.r)**2:
            values[0] = self.value[0]/(pi*self.r**2)
            values[1] = self.value[1]/(pi*self.r**2)
        else:
            values[0] = 0
            values[1] = 0
    def value_shape(self):
        return (2,)
    
# add loads together
i = 0
f = PatchLoad(pt=(0,0), vl=(0,0), r = 0,degree = 1)
while i<len(x_p):
    z = PatchLoad(pt=(x_p[i],y_p[i]), vl=(f_x[i],f_y[i]), r = 1 ,degree = 1)
    f = f+z
    i += 1
    
# two translations + one rotation
Z = [Constant((1, 0)), Constant((0, 1)), Expression(('-x[1]', 'x[0]'), degree=1)]
V = VectorFunctionSpace(mesh, "Lagrange", 1)     # Space for displacement
M = VectorFunctionSpace(mesh, 'R', 0, 3)   # Space for all multipliers
Vs = [V, M]
W = MixedFunctionSpace(*Vs)

u, p = TrialFunctions(W)
w, q = TestFunctions(W)

def eps(w):
    return sym(grad(w))
E = Constant(1.5e-3)
nu = Constant(0.47)
model = "plane_stress"
mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)
if model == "plane_stress":
    lmbda = 2*mu*lmbda/(lmbda+2*mu)
def sigma(w):
    return lmbda*tr(eps(w))*Identity(2) + 2.0*mu*eps(w)
T = Constant((0,0))
a = inner(sigma(u), eps(w))*dx     -sum(p[i]*inner(w, Z[i])*dx for i in range(len(Z)))    -sum(q[i]*inner(u, Z[i])*dx for i in range(len(Z))) 

l = inner(f, w)*dx + inner(T, w)*ds
kk = Function(W)
solve(a == l, kk)
(u, p) = kk.split()
disp = interpolate(u,V)
s= -sigma(u)


# NOTE: if I invert y axis, y values are inverted also!!!
plt.rcParams["figure.figsize"] = (15,15)
c = plot(disp,title='deformation vector μm')
plt.colorbar(c)
cmap = plt.get_cmap('jet')
plt.set_cmap(cmap)
plt.rcParams.update({'font.size': 22})
#plt.gca().invert_yaxis()
plt.show()

plt.rcParams["figure.figsize"] = (5,5)
c = plot(u, mode="displacement",title='shape of cell at stress free state')
#plt.colorbar(c)
cmap = plt.get_cmap('jet')
plt.set_cmap(cmap)
plt.gca().invert_yaxis()
plot(mesh)
plt.rcParams.update({'font.size': 14})

I1 = s[0,0] + s[1,1]
I2 = s[0,0]*s[1,1]-s[0,1]*s[0,1]

fig, ax = plt.subplots(figsize=(20,20))
plt.subplot(1, 2, 1)
c = plot(s[0,0]*1e6+s[1,1]*1e6, title='First Stress Invariant Pa')
plt.colorbar(c,fraction=0.046, pad=0.06)
cmap = plt.get_cmap('jet')
plt.set_cmap(cmap)
plt.gca().invert_yaxis()
plt.subplots_adjust(wspace=0.4)
c = plt.subplot(1, 2, 2)
c = plot(sqrt(I1*I1 - 3*I2)*1e6, title='von Mises stress Pa')
plt.colorbar(c,fraction=0.046, pad=0.06)
cmap = plt.get_cmap('jet')
plt.set_cmap(cmap)
plt.rcParams.update({'font.size': 22})
plt.gca().invert_yaxis()

