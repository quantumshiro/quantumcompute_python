############################################################################
#　「量子井戸内の固有状態」アニメーション
############################################################################
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#図全体
fig = plt.figure(figsize=(10, 6))
#全体設定
plt.rcParams['font.family'] = 'Times New Roman' #フォント
plt.rcParams['font.size'] = 12 #フォントサイズ

#複素数
I = 0.0 + 1.0j

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
#量子井戸の幅
L = 1.0 * 10**-9
#計算区間
x_min = -L / 2.0
x_max = L / 2.0
#状態数
n_max = 5
#空間分割数
NX = 500
#空間刻み間隔
dx = 1.0 * 10**-9
#計算時間の幅
ts = -50
te = 50
#固有関数
def verphi(n, x):
	kn = math.pi * (n + 1) / L
	return math.sqrt(2.0 / L) * math.sin(kn * (x + L / 2.0))
#固有エネルギー
def Energy(n):
	kn = math.pi * (n + 1) / L
	return hbar * hbar * kn**2 / (2.0 * me)
#波動関数
def phi(n, x, t):
	kn = math.pi * (n + 1) / L
	omega = hbar * kn**2 / (2.0 * me)
	return verphi(n,x) * cmath.exp(- I * omega * t)

#基底状態の周期
T0 = 2.0 * math.pi * hbar / Energy(0)
#時間間隔
dt = T0 / (te - ts + 1)

#アニメーション作成用
ims=[]

#各時刻における波束の計算
for tn in range(ts,te):
	#実時間の取得
	t = dt * tn
	xl = []
	#２重配列の初期化
	nl = [0] * (n_max + 1)
	for j in range( n_max + 1 ):
		nl[ j ] = []

	#各空間地点ごとに計算
	for ix in range(NX):
		#空間地点の座標を取得
		x = x_min + (x_max - x_min) * ix / NX
		xl.append(x / dx)
		#エネルギー順位ごとに
		for n in range(n_max+1):
			nl[ n ].append( phi(n, x, t).real / math.sqrt(2.0 / L) * 0.5 + n )

	#各コマを描画
	img  = plt.plot(xl, nl[0], "blue")
	img += plt.plot(xl, nl[1], "green")
	img += plt.plot(xl, nl[2], "red")
	img += plt.plot(xl, nl[3], "yellow")
	img += plt.plot(xl, nl[4], "black")
	img += plt.plot(xl, nl[5], "cyan")
	#アニメーションに追加
	ims.append( img )

### グラフ描画
plt.title("QuantumWell", fontsize = 16)
plt.xlabel("Position[nm]", fontsize = 16)
plt.ylabel("Eigenfunction", fontsize = 16)
#描画範囲を設定
plt.xlim([-0.5,0.5])
#アニメーションの生成
ani = animation.ArtistAnimation(fig, ims, interval=10)
#グラフの表示
plt.show()
