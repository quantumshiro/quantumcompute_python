'''鉛直投げ上げ運動シュミレーション'''
import numpy as np
from matplotlib import pyplot
import math

'''物理パラメータ'''
#時刻の範囲
t_min=0
t_max=4.0
#初期位置
x0=1.0
#時間間隔
dt=0.01
dtr=0.015625
#データ間引き数
skip=1
#格納用変数
x=[]
xrk=[]
#データ取得回数
N=round((t_max-t_min)/dt)
Nr=round((t_max-t_min)/dtr)
'''数値計算パラメータ'''
def dxdt(x):
    return x

'''ルンゲクッタ'''
xp=x0
i=0
for i in range(Nr):
    xrk=np.append(xrk,xp)
    k1=dxdt(xp)
    k2=dxdt(xp+0.5*k1*dtr)
    k3=dxdt(xp+0.5*k2*dtr)
    k4=dxdt(xp+k3*dtr)
    k=(k1+2*k2+2*k3+k4)/6
    xp=xp+k*dtr

'''解析解'''
def X(t):
    return math.e**t

i=0
for i in range(N):
    x=np.append(x,X(i*dt))

'''描画'''
t=np.linspace(t_min,t_max,N)
trk=np.linspace(t_min,t_max,Nr)
trk=trk[::skip]
xrk=xrk[::skip]
print(X(t_max)-xp)
pyplot.title("Vertical throwing")
pyplot.xlabel("Time")
pyplot.ylabel("position&verocity")
pyplot.plot(t,x,label="X")
pyplot.scatter(trk,xrk,label="Xrk4")
pyplot.legend()
pyplot.show()
