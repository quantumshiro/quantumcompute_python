############################################################################
#　状態遷移シミュレーション（基底状態→励起状態）
############################################################################
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate
import numpy as np
import numpy.linalg as LA

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
#電気素量
e = 1.60217733 * 10**-19
#光速
c = 2.99792458E+8
######################################
#　物理系の設定
######################################
#量子井戸の幅
L = 1 * 10**-9
#計算区間
x_min = -L / 2.0
x_max = L / 2.0
#状態数
n_max = 3
#行列の要素数
DIM = n_max + 1
#空間分割数
NX = 500
#空間刻み間隔
dx = 1.0 * 10**-9
#時間間隔
dt = 10**-18
#計算ステップ数
Tn = 3000000
#データ間引き数
skip = 10000

#入射電磁波ベクトルポテンシャルの振幅
A0 = 1.0E-8

#固有関数
def verphi(n, x):
	kn = math.pi * (n + 1) / L
	return math.sqrt(2.0 / L) * math.sin(kn * (x + L / 2.0))

#固有エネルギー
def Energy(n):
	kn = math.pi * (n + 1) / L
	return hbar * hbar * kn * kn / (2.0 * me)

#被積分関数（Xmn）の計算
def integral_Xnm(x, n1, n2):
	return verphi(n1 ,x) * x * verphi(n2, x)

#エネルギー差
dE = Energy(1) - Energy(0)
#エネルギー差に対応する光子の角振動数
omega = dE / hbar
#電磁波の波長
_lambda =  2.0 * math.pi * c / omega

print( "エネルギー（基底状態）" + str( Energy(0) / eV ) + "[eV]" )
print( "エネルギー（励起状態）" + str( Energy(1) / eV ) + "[eV]" )
print( "エネルギー差" + str( dE /eV ) + "[eV]" )
print( "エネルギー差に対応する光子の角振動数" + str( omega ) + "[rad/s]" )
print( "エネルギー差に対応する光子の角振動数に対する周期" + str( 2.0 * math.pi / omega  ) + "[s]" )
print( "電磁波の波長" + str( _lambda / 1.0E-9  ) + "[nm]" )

######################################
#　ルンゲクッタクラス
######################################
class RK4:
	#コンストラクタ
	def __init__(self, N, dt):
		self.dt = dt
		self.N = N
		self.A0 = 0
		self.omega = 0
		self.bn  = np.array([0+0j] * N) 
		self.dbn = np.array([0+0j] * N)
		#２重配列の準備
		Xnm = [0+0j] * N 
		for i in range( N ):
			Xnm[ i ] = [0+0j] * N
		self.Xnm = np.array(Xnm)

		self.__a1 = np.array([0+0j] * N)
		self.__a2 = np.array([0+0j] * N)
		self.__a3 = np.array([0+0j] * N)
		self.__a4 = np.array([0+0j] * N)

	#１階微分を与えるメソッド
	def Db(self, t, bn, out_bn ):
		for n in range(self.N):
			#対角成分
			out_bn[n] = Energy(n) / (I * hbar) * bn[n] 
			#非対角成分
			for m in range(self.N) :
				out_bn[n] += self.A0 * e / hbar / hbar * math.cos( self.omega * t ) * (Energy(n) - Energy(m)) * self.Xnm[n][m] * bn[m]

	#時間発展を計算するメソッド
	def timeEvolution(self, t):
		#１段目
		self.Db( t, self.bn, self.__a1 )
		#２段目
		self.Db( t, self.bn + self.__a1 * 0.5 * self.dt, self.__a2 )
		#３段目
		self.Db( t, self.bn + self.__a2 * 0.5 * self.dt, self.__a3 )
		#４段目
		self.Db( t, self.bn + self.__a3 * self.dt, self.__a4 )
		#差分の計算
		self.dbn = (self.__a1 + 2.0 * self.__a2 + 2.0 * self.__a3 + self.__a4) * self.dt / 6.0

######################################
#　計算開始
######################################

#ルンゲクッタクラスのインスタンスを生成
rk4 = RK4(DIM, dt)
#電磁波の角振動数
rk4.omega = omega
#入射電磁波ベクトルポテンシャルの振幅
rk4.A0 = A0

### Xnm の計算 
for n1 in range(n_max + 1):
	for n2 in range(n_max + 1):
		#ガウス・ルジャンドル積分
		result = integrate.quad(
			integral_Xnm, #被積分関数
			x_min, x_max,		   #積分区間の下端と上端
			args=(n1, n2)			#被積分関数へ渡す引数
		)
		real = result[0]
		imag = 0
		#行列要素
		rk4.Xnm[n1][n2] = real + 1j * imag

		if( abs(real / L) < L ): real = 0
		#ターミナルへ出力
		print( "(" + str(n1) + ", " + str(n2) + ")  " + str( real / L ))

#初期状態の設定
rk4.bn = np.array(
	[ 1.0+0.0j,   #基底状態
	  0.0+0.0j,   #第１励起状態 
	  0.0+0.0j,   #第２励起状態 
	  0.0+0.0j ]  #第３励起状態
)

#グラフ描画用配列の準備
ts = []
b0s = []
b1s = []

###展開係数の時間依存性の計算
for tn in range(Tn+1):
	#実時間
	t_real = dt * tn
	#計算結果の取得
	if( tn % skip == 0 ): 
		print("t =" + str(tn / skip) + "  " + str(abs(rk4.bn[0])**2))
		#時刻の取得
		ts.append( tn / skip )
		#基底状態の存在確率
		b0s.append( abs(rk4.bn[0])**2 )
		#第１励起状態の存在確率
		b1s.append( abs(rk4.bn[1])**2 )

	#ルンゲ・クッタ法による時間発展
	rk4.timeEvolution( t_real )
	#展開係数を更新
	rk4.bn += rk4.dbn

### グラフ描画
plt.title("Expansion coefficient at time")
plt.xlabel("time[s]")
plt.ylabel("Expansion coefficient")
#描画範囲を設定
plt.xlim([0,Tn/skip])
plt.ylim([0,1])
#グラフ描画
plt.plot(ts, b0s, marker="o" , linewidth=3.0)
plt.plot(ts, b1s, marker="o",  linewidth=3.0)
#グラフ表示
plt.show()
