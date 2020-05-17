import math
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)
# x = 0.3 * np.outer(np.cos(u), np.sin(v))
# y = 0.3 * np.outer(np.sin(u), np.sin(v))
# z = 0.3 * np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.5)

k = 9e+9
Q = 7e-2
q = 1e-3
m = 1e-3
V0 = 10000
dt = 5e-8

angleA = 30
angleB = 45
alpha = angleA * math.pi / 180
beta = angleB * math.pi / 180
x = -0.2
y = -0.2
z = -0.2
t = 0
Vx = V0*math.cos(alpha)*math.cos(beta)
Vy = V0*math.sin(alpha)*math.cos(beta)
Vz = V0*math.sin(beta)

for i in range(100):
    Ax = k * Q * q / x ** 2 / m
    Ay = k * Q * q / y ** 2 / m
    Az = k * Q * q / z ** 2 / m
    Vx += Ax * dt
    Vy += Ay * dt
    Vz += Az * dt
    x += (Vx * dt + Ax * dt ** 2 / 2)
    y += (Vy * dt + Ay * dt ** 2 / 2)
    z += (Vz * dt + Az * dt ** 2 / 2)
    print(Ax, Ay, Az, Vx, Vy, Vz, x, y, z)
    ax.scatter([x], [y], [z], c='red')
    if math.sqrt(x**2+y**2+z**2)<0.05:
        break

ax.scatter(0, 0, 0, c='blue')
plt.show()