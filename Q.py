import math
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

q = 1.6e-19
m = 9.1e-31
V0 = 20
dt = 1e-12

angleA = 45
angleB = 70
alpha = angleA * math.pi / 180
beta = angleB * math.pi / 180
x = 0
y = 0
z = 0
Vxy = V0*math.sin(beta)
Vz = V0*math.cos(beta)
B = 1
ax.scatter([x], [y], [z], c='blue', s=100)
w = q*B/m
t=0
for i in range(100):
    t += dt
    Vx = Vxy * math.cos(w*t)
    Vy = -Vxy * math.sin(w*t)
    x += Vx * dt
    y += Vy * dt
    z += Vz * dt

    print(round(Vx), round(Vy), round(Vz), round(x, 2), round(y, 2), round(1000000*z,5))
    if math.cos(w*t)>=0:
        ax.scatter([x], [y], [z], c='red', s=0.3)
    if math.cos(w*t)<0:
        ax.scatter([x], [y], [z], c='blue', s=0.3)
plt.show()