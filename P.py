import math
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

#ax.set_xlim(-100, 100)
#ax.set_ylim(-100, 100)
#ax.set_zlim(-100, 100)
#circle1 = plt.Circle((0, 0), 5, color='blue', fill=True)
#ax.add_patch(circle1)

k = 9e+9
Q = 7e-2
q = 1e-3
m = 1e-3
V0 = 10
dt = 5e-9

angleA = 30
angleB = 0
alpha = angleA * math.pi / 180
beta = angleB * math.pi / 180
x = -20
y = -20
z = 0
t = 0
Vx = V0*math.cos(alpha)*math.cos(beta)
Vy = V0*math.sin(alpha)*math.cos(beta)
Vz = V0*math.sin(beta)

for i in range(1000):
    r = math.sqrt(x ** 2 + y ** 2 + z ** 2) / 100
    alpha1 = math.atan2(abs(y), abs(x))
    beta1 = math.atan2(abs(z), abs(y))

    if (x<0 and y>0 and z>0):
        alpha1 = math.pi - alpha1
    if (x<0 and y<0 and z>0):
        alpha1 = math.pi + alpha1
        beta1 = math.pi - beta1
    if (x>0 and y<0 and z>0):
        alpha1 = -alpha1
        beta1 = math.pi - beta1
    if (x<0 and y>0 and z<0):
        alpha1 = math.pi - alpha1
        beta1 = -beta1
    if (x<0 and y<0 and z<0):
        alpha1 = math.pi + alpha1
        beta1 = math.pi + beta1
    if (x>0 and y>0 and z<0):
        beta1 = -beta1
    if (x>0 and y<0 and z<0):
        alpha1 = alpha1 = -alpha1
        beta1 = math.pi + beta1

    Ax = k * Q * q / r ** 2 / m * math.cos(alpha1)*math.cos(beta1)
    Ay = k * Q * q / r ** 2 / m * math.sin(alpha1)*math.cos(beta1)
    Az = k * Q * q / r ** 2 / m * math.sin(beta1)
    Vx += Ax*dt
    Vy += Ay*dt
    Vz += Az*dt
    x += (Vx * dt + Ax * dt ** 2 / 2) * 100
    y += (Vy * dt + Ay * dt ** 2 / 2) * 100
    z += (Vz * dt + Az * dt ** 2 / 2) * 100
    print(Ax, Ay, Az, Vx, Vy, Vz, x, y, z, r, alpha1*180/math.pi, beta1*180/math.pi)
    ax.scatter([x], [y], [z], c='red')
ax.scatter(0, 0, 0, c='red')
plt.show()