from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

mu,sigma,n = 30,12,50
xs = np.random.normal(mu,sigma,n)
ys = np.random.normal(mu,sigma,n)
zs = np.random.normal(mu,sigma,n)
ax.scatter(xs, ys, zs, c='black', marker='o')

mu,sigma,n = 70,12,50
xs = np.random.normal(mu,sigma,n)
ys = np.random.normal(mu,sigma,n)
zs = np.random.normal(mu,sigma,n)
ax.scatter(xs, ys, zs, c='black', marker='o')

ax.set_xlabel('Trait')
ax.set_ylabel('Environment')
ax.set_zlabel('Outcome')
ax.set_xlim([0, 100])
ax.set_ylim([0, 100])
ax.set_zlim([0, 100])
# plt.show()


plt.figure(figsize=(5,4.5))
# plt.figure(figsize=(6,4))
n=100
x_arr = np.linspace(0,100,n)
# y_arr_det = 100*np.exp(-0.04*x_arr) 
y_arr_det = 100*(1.0 - 1.0/(1.0+np.exp(0.1*(x_arr - 50) )))

y_arr_sto = y_arr_det + np.random.normal(0,20,n)
# plt.scatter(x_arr, y_arr_sto, marker='o', c='black')
plt.plot(x_arr, y_arr_det, c='black')
plt.ylim(-20,120)
plt.show()
