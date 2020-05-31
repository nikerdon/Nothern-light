import math
import matplotlib.pyplot as plt
import numpy as np

q = 1
m = 0.0000005
V0 = 200
dt = 1
Rz = 6370

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# r = Rz
# pi = np.pi
# cos = np.cos
# sin = np.sin
# phi1, theta1 = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
# x1 = r*sin(phi1)*cos(theta1)
# y1 = r*sin(phi1)*sin(theta1)
# z1 = r*cos(phi1)
# ax.plot_surface(x1, y1, z1,  rstride=1, cstride=1, color='blue', alpha=0.6, linewidth=0)

n = 1
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

x = 2*Rz *1000
y= 0.9*Rz *1000
z= Rz *1000

Vx = -100
Vy = 0
Vz = 0
ax.scatter([x/1000], [y/1000], [z/1000], c='green', s=300)

for i in range(30000):
    r = math.sqrt(x**2+y**2+z**2)
    V = math.sqrt(Vx**2+Vy**2+Vz**2)

    r1 = 10000 * math.pow(npr+4*na,-1/6)*math.pow(V,-1/3)*Rz*1000
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

    Ax = q * (Vy*bz - Vz*by) / m
    Ay = q * (Vz*bx - Vx*bz) / m
    Az = q * (Vx*by - Vy*bx) / m
    Vx += Ax*dt
    Vy += Ay*dt
    Vz += Az*dt
    x += Vx * dt + Ax*dt**2/2
    y += Vy * dt + Ay*dt**2/2
    z += Vz * dt + Az*dt**2/2
    #print(n,r1,x,y,z)
    #print(n, round(B1X*10**12), round(B1Y*10**12), round(B1Z*10**12),'||', round(Ax), round(Ay), round(Az), '||',round(Vx), round(Vy), round(Vz), '||',round(x/1000), round(y/1000), round(z/1000))
    n+=1
    if n%10==0:
        ax.scatter([x / 1000], [y / 1000], [z / 1000], c='red', s=0.3)
plt.show()

# theta = math.acos(abs(Vx*x+Vy*y+Vz*z)/(V*r))
#     gamma = math.acos(z/r)
#     delta = math.acos(x/(r*math.sin(gamma)))
#     k = (1e-7)*q/r**3
#     Bq = k*V*r*math.sin(theta)
#     B0 = math.sqrt((B1X+B*x/r+k*(Vy*z-Vz*y))**2 + (B1Y+B*y/r+k*(Vz*x-Vx*z))**2 + (B1Z+B*z/r+k*(Vx*y-Vy*x))**2)
#     Bx = B1X+B*x/r+k*(Vy*z-Vz*y)
#     By = B1Y+B*y/r+k*(Vz*x-Vx*z)
#     Bz = B1Z+B*z/r+k*(Vx*y-Vy*x)
#     #print(Bx, By, Bz, B0)
#     gamma = math.atan2(By,Bx)
#     delta = math.pi/2 - math.atan2(Bz,(math.sqrt(By**2+Bx**2)))
#     A0 = abs(q) * V * B0 / m * math.sqrt(1 - (abs(Bx*Vx+By*Vy+Bz*Vz)/V/B0)**2)
#     al = math.acos((Vx*y - Vy*x)/math.sqrt(1 - (abs(Bx*Vx+By*Vy+Bz*Vz)/V/B0)**2))
#     be = math.acos((Vy*z-Vz*y)/math.sqrt(1 - (abs(Bx*Vx+By*Vy+Bz*Vz)/V/B0)**2)/math.sin(al))
#     Ax = A0 * math.sin(al) * math.cos(be)
#     Ay = A0 * math.sin(al) * math.sin(be)
#     Az = A0 * math.cos(al)
#     # Ax = abs(q) * V * B0 * (Vy*Bz - Vz*By) / m
#     # Ay = abs(q) * V * B0 * (Vz*Bx - Vx*Bz) / m
#     # Az = abs(q) * V * B0 * (Vx*By - Vy*Bx) / m
#     print(Ax,Ay,Az)
#     t += dt
#     Vx += Ax*dt
#     Vy += Ay*dt
#     Vz += Az*dt
#     x += Vx * dt + Ax*dt**2/2
#     y += Vy * dt + Ay*dt**2/2
#     z += Vz * dt + Az*dt**2/2
#     ax.scatter([x/1000], [y/1000], [z/1000], c='red', s=0.3)
#     #print((Ax), (Ay), (Az), round(x/1000, 2), round(y/1000, 2), round(z/1000,2))
#     #if math.cos(w*t)>=0:
#         #ax.scatter([x], [y], [z], c='red', s=0.3)
#     #if math.cos(w*t)<0:
#         #ax.scatter([x], [y], [z], c='blue', s=0.3)