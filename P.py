import math
import matplotlib.pyplot as plt

ax = plt.gca()
#plt.axis('scaled')
ax.grid()
#ax.set_ylim(-1, 1)
#ax.set_xlim(-200, 200)
#circle1 = plt.Circle((0, 0), 5, color='blue', fill=True)
#ax.add_patch(circle1)

k = 9e+9
Q = 7e-2
q = 1e-3
m = 1e-3
V0 = 1
dt = 5e-12

angle=45
alpha = angle * math.pi / 180
x=-65
y=3
t=0
Vx = V0*math.cos(alpha)
Vy = V0*math.sin(alpha)

for i in range(1000):
    r = math.sqrt(x ** 2 + y ** 2) / 100
    beta = math.atan2(abs(y),abs(x))
    if (x<0 and y>0):
        beta = math.pi - beta
    if (x<0 and y<0):
        beta = math.pi + beta
    if (x>0 and y<0):
        beta = -beta
    Ax = k * Q * q / r**2 / m * math.cos(beta)
    Ay = k * Q * q / r**2 / m * math.sin(beta)
    Vx += Ax*dt
    Vy += Ay*dt
    x += (Vx*dt + Ax*dt**2/2)*100
    y += (Vy*dt + Ay*dt**2/2)*100
    #print(Ax, Ay, Vx, Vy, x, y, r, beta*180/math.pi)
    plt.scatter([x], [y], s=1, c='red')
plt.show()

# Ax = - Q * q * x / (m * math.pow((x ** 2 + y ** 2), 3 / 2))
# Ay = - Q * q * y / (m * math.pow((x ** 2 + y ** 2), 3 / 2))
# Vx = Q * q * (1 / (m * math.sqrt(x ** 2 + y ** 2) - 1 / (m * math.sqrt(x0 ** 2 + y0 ** 2)))) - Vx0
# Vy = Q * q * (1 / (m * math.sqrt(x ** 2 + y ** 2) - 1 / (m * math.sqrt(x0 ** 2 + y0 ** 2)))) - Vy0
# x1 = x
# y1 = y
# x = Q * q / m * math.log(math.fabs((math.sqrt(x1 ** 2 + y1 ** 2) + x1) / (math.sqrt(x0 ** 2 + y0 ** 2)))) - Vx0 * t - x0
# y = Q * q / m * math.log(math.fabs((math.sqrt(x1 ** 2 + y1 ** 2) + y1) / (math.sqrt(x0 ** 2 + y0 ** 2)))) - Vy0 * t - y0
# t += dt

# Ax = k * Q * q / r**2 / m * math.cos(math.pi-math.tan(y/x))
# Ay = k * Q * q / r**2 / m * math.sin(math.pi-math.tan(y/x))
# Vx += Ax * dt
# Vy += Ay * dt
# x += Vx * dt + Ax*dt**2/2
# y += Vy * dt + Ay*dt**2/2
# t += dt
# r = math.sqrt(x ** 2 + y ** 2)
# print(Ax,Ay,Vx,Vy,x,y,alpha)
# plt.scatter([x], [y], s=1, c='red')
# alpha = math.tan(y / x)
# # if x>=0 and y>=0:
# #     alpha = math.tan(y/x)
# # if x<0 and y>0:
# #     alpha = math.pi-math.tan(y/x)
# # if x>0 and y<0:
# #     alpha = math.tan(y/x)
# # if x<0 and y<0:
# #     alpha = math.pi-math.tan(y/x)
# #if (x==0 or y==0 or abs(x)>30 or abs(y)>30):
# #    break