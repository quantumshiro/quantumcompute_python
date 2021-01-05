'''鉛直投げ上げ運動シュミレーション'''
import numpy as np
import math
from matplotlib import pyplot

'''物理パラメータ'''
#時刻の範囲
t_min=0
t_max=15.0
#質量
m=1.0
#角速度
omega=2.0*math.pi/10.0
#ばね定数
kb=m*omega*omega
#初期位置
x0=1.0
#初速度
v0=0.0
#時間間隔
dt=0.01
dtr=0.01
#データ間引き数
skip=20
#格納用変数
v=[]
x=[]
vrk=[]
xrk=[]
#データ取得回数
N=round((t_max-t_min)/dt)
Nr=round((t_max-t_min)/dtr)


'''数値計算パラメータ'''
def dvdt(x,v):
	return (-kb/m*x)
def dxdt(x,v):
	return (v)

'''ルンゲクッタ'''
vp=v0
xp=x0
for i in range(Nr):
	xrk=np.append(xrk,xp)
	vrk=np.append(vrk,vp)

	k1=dxdt(xp,vp)
	j1=dvdt(xp,vp)
	k2=dxdt(xp+0.5*k1*dtr,vp+0.5*j1*dtr)
	j2=dvdt(xp+0.5*k1*dtr,vp+0.5*j1*dtr)
	k3=dxdt(xp+0.5*k2*dtr,vp+0.5*j2*dtr)
	j3=dvdt(xp+0.5*k2*dtr,vp+0.5*j2*dtr)
	k4=dxdt(xp+k3*dtr,vp+j3*dtr)
	j4=dvdt(xp+k3*dtr,vp+j3*dtr)
	k=(k1+2*k2+2*k3+k4)/6
	j=(j1+2*j2+2*j3+j4)/6

	xp=xp+k*dtr
	vp=vp+j*dtr

'''解析解'''
def X(t):
	return x0*math.cos(omega*t)+v0/omega*math.sin(omega*t)
def V(t):
	return v0*math.cos(omega*t)-x0*omega*math.sin(omega*t);

i=0
for i in range(N):
	x=np.append(x,X(i*dt))
	v=np.append(v,V(i*dt))

'''描画'''
print(xp-X(8))
t=np.linspace(t_min,t_max,N)
trk=np.linspace(t_min,t_max,Nr)
xrk=xrk[::skip]
vrk=vrk[::skip]
trk=trk[::skip]
pyplot.title("Simple vibration")
pyplot.xlabel("Time")
pyplot.ylabel("position&verocity")
pyplot.plot(trk,xrk,label="Xrk4")
pyplot.plot(trk,vrk,label="Vrk4")
pyplot.legend()
pyplot.show()
