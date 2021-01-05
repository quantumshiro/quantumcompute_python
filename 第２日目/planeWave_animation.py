############################################################################
#　平面波アニメーション
############################################################################
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#図全体
fig = plt.figure(figsize=(10, 6))
#全体設定
plt.rcParams['font.family'] = 'Times New Roman' #フォント
plt.rcParams['font.size'] = 12 #フォントサイズ

######################################
#　物理定数
######################################
#プランク定数
h = 6.6260896 * 10**-34
hbar = h / (2.0 * math.pi)
#電子の質量
me = 9.10938215 * 10**-31
#電子ボルト
eV = 1.60217733 * 10**-19

######################################
#　物理系の設定
######################################
#電子のエネルギー
E1 = 0.25 * eV
E2 = 1.0 * eV
E3 = 4.0 * eV
#波数
k1 = math.sqrt(2.0 * me * E1 / (hbar * hbar))
k2 = math.sqrt(2.0 * me * E2 / (hbar * hbar))
k3 = math.sqrt(2.0 * me * E3 / (hbar * hbar))
#角振動数
omega1 = E1/hbar
omega2 = E2/hbar
omega3 = E3/hbar
#時間間隔
dt = 10**-16
#空間刻み間隔
dx = 10**-9
#空間刻み数
XN = 400
#時間刻み数
TN = 1000
#空間幅
x_min = -2.0
x_max = 2.0

#アニメーション作成用
ims=[]
	
#各時刻における計算
for tn in range(TN):
	t = dt * tn
	#データ格納用
	xl=[]
	psi1l=[]
	psi2l=[]
	psi3l=[]
	for ix in range(XN):
		x = (x_min + (x_max - x_min) * ix/XN) * dx
		psi1 = math.cos(k1 * x - omega1 * t)
		psi2 = math.cos(k2 * x - omega2 * t)
		psi3 = math.cos(k3 * x - omega3 * t)
		#描画用データ生成
		xl.append(x/dx)
		psi1l.append(psi1)
		psi2l.append(psi2)
		psi3l.append(psi3)

	#各コマを描画
	img  = plt.plot(xl, psi1l, color='red',	  linestyle='solid', linewidth = 3.0, label='E_1')
	img += plt.plot(xl, psi2l, color='green', linestyle='solid', linewidth = 3.0, label='E_2')
	img += plt.plot(xl, psi3l, color='blue',  linestyle='solid', linewidth = 3.0, label='E_3')
	ims.append( img )

#グラフタイトルの設定 
plt.title("Plane wave")
#x軸ラベル・y軸ラベルの設定 
plt.xlabel("Position")
plt.ylabel("Real part of Wave fanction")
#描画範囲を設定
plt.xlim([-2.0,2.0])
plt.ylim([-1.2,1.2])
#アニメーションの生成
ani = animation.ArtistAnimation(fig, ims, interval=10)
#アニメーションの保存
ani.save("output.html", writer=animation.HTMLWriter())
#ani.save("output.gif", writer="imagemagick")
#ani.save("output.mp4", writer="ffmpeg", dpi=300)
#グラフの表示
plt.show()
