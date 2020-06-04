import math
import matplotlib.pyplot as plt
import numpy as np

k = 9e+9
q = 1
Q = -7e-2
m = 0.0000005
dt = 1
Rz = 6370

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

r = Rz
pi = np.pi
cos = np.cos
sin = np.sin
phi1, theta1 = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
x1 = r*sin(phi1)*cos(theta1)
y1 = r*sin(phi1)*sin(theta1)
z1 = r*cos(phi1)
ax.plot_surface(x1, y1, z1,  rstride=1, cstride=1, color='blue', alpha=0.6, linewidth=0)

n = 163
UT = 10.6
npr = 5 * 1e6
na = 2.5 * 1e5
dalpha1 = 11
dalpha2 = 23.5
alpha1 = dalpha1 * math.pi / 180
alpha2 = dalpha2 * math.pi / 180

beta = math.pi/180 * (23.5*math.sin(360/365*(284+n)/180*math.pi))
phi = (15*UT-69)*math.pi/180
psi = math.asin(-math.sin(beta) * math.cos(alpha1) + math.cos(beta) * math.sin(alpha1) * math.cos(phi))
si=math.sin(psi)
co=math.cos(psi)

x =Rz * 1000 *5
y= Rz * 1000 *0
z= Rz * 1000 *1

Vx = -100000
Vy = 0
Vz = 0
ax.text2D(0.03, 0.83, str(q)+' Кл'+'\n'+str(m)+' кг'+'\n'+'Vx = '+str(Vx)+' м/с \n'+'Vy = '+str(Vy)+' м/с \n'+'Vz = '+str(Vz)+' м/с \n', transform=ax.transAxes)
ax.text2D(0.88, 0.86, 'x0 = '+ str(x/1000)+' км \n'+'y0 = '+str(y/1000)+' км \n'+'z0 = '+str(z/1000)+' км \n', transform=ax.transAxes)
ax.scatter([x/1000], [y/1000], [z/1000], c='green', s=100)

for i in range(1000):
    r = math.sqrt(x**2+y**2+z**2)
    V = math.sqrt(Vx**2+Vy**2+Vz**2)

    #r1 = 10000 * math.pow(npr+4*na,-1/6)*math.pow(400000,-1/3)*Rz*1000
    r1 = r
    B1X = (8.5*si - 39.6*si*co*x/r1 + 1.2*si*y/r1 + (21.8*co**2-17.9*si**2)*z/r1 + psi/10*180/math.pi*(2.9*co + (-3*co**2+5.5*si**2)*x/r1 + 0.2*co*y/r1 - 8.5*si*co*z/r1))*(1e-9)
    B1Y = (180/math.pi*psi/10*(-0.2*co*x/r1 -2.5 * y/r1 - 0.2*si*z/r1))*(1e-9)
    B1Z = (-8.5*co + (-21.8*si**2+17.9*co**2)*x/r1 - 1.2*co*y/r1 + 39.6*si*co*z/r1 + psi/10*180/math.pi*(2.9*si-8.5*si*co*x/r1 + 0.2*si*y/r1 + (-3*si**2+5.5*co**2)*z/r1))*(1e-9)
    B1 = math.sqrt(B1X**2+B1Y**2+B1Z**2)

    B = 4/3 * 2*math.pi/4 * 4 * 1e-7 * r * (-1.15*(1e-9)) * 2 *math.pi/24/3600

    gamma = math.acos(z/r)
    theta = math.acos(x/r)
    BX = B * math.sin(gamma) * math.cos(theta)
    BY = B * math.sin(gamma) * math.sin(theta)
    BZ = B * z/r

    bx = B1X + BX
    by = B1Y + BY
    bz = B1Z + BZ
    bo = math.sqrt(bx**2+by**2+bz**2)

    alpha1 = math.atan2(y, x)
    beta1 = math.acos(z/r)

    Ax = q * (Vy*bz - Vz*by) / m + k * Q * q / r ** 2 / m * math.cos(alpha1)*math.sin(beta1)
    Ay = q * (Vz*bx - Vx*bz) / m + k * Q * q / r ** 2 / m * math.sin(alpha1)*math.sin(beta1)
    Az = q * (Vx*by - Vy*bx) / m + k * Q * q / r ** 2 / m * math.cos(beta1)

    Vx += Ax*dt
    Vy += Ay*dt
    Vz += Az*dt
    x += Vx * dt + Ax*dt**2/2
    y += Vy * dt + Ay*dt**2/2
    z += Vz * dt + Az*dt**2/2
    print(n)
    if r <= Rz*1000:
        ax.scatter([x/1000], [y/1000], [z/1000], c='magenta', s=100)
        print(x/1000,y/1000,z/1000)
        break
    #print(n,r1,x,y,z)
    #print(n, round(B1X*10**12), round(B1Y*10**12), round(B1Z*10**12),'||', round(Ax), round(Ay), round(Az), '||',round(Vx), round(Vy), round(Vz), '||',round(x/1000), round(y/1000), round(z/1000))
    n+=1
    if n%10==0:
        ax.scatter([x / 1000], [y / 1000], [z / 1000], c='red', s=0.2)
# for angle in range(0, 360):
#    ax.view_init(30, angle)
#    plt.draw()
#    plt.pause(.0001)
plt.show()
