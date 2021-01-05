############################################################################
#　基底状態の単振動の確認
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

######################################
#　物理系の設定
######################################
#量子井戸の幅
L = 1 * 10**-9
#計算区間
x_min = -L / 2.0
x_max = L / 2.0
#状態数
n_max = 2
#行列の要素数
DIM = n_max + 1
#空間分割数
NX = 500
#空間刻み間隔
dx = 1.0 * 10**-9
#時間間隔
dt = 10**-16
#計算ステップ数
Tn = 300
#データ間引き数
skip = 1

#固有関数
def verphi(n, x):
	kn = math.pi * (n + 1) / L
	return math.sqrt(2.0 / L) * math.sin(kn * (x + L / 2.0))

#固有エネルギー
def Energy(n):
	kn = math.pi * (n + 1) / L
	return hbar * hbar * kn * kn / (2.0 * me)

#基底状態の周期
T0 = 2.0 * math.pi * hbar / Energy(0)
print( "E0 = " + str(Energy(0)/eV) + "[eV]")
print( "T0 = " + str(T0) + "[s]")

######################################
#　ルンゲクッタクラス
######################################
class RK4:
	#コンストラクタ
	def __init__(self, DIM, dt):
		self.dt = dt
		self.DIM = DIM
		self.bn  = np.array([0+0j] * DIM) 
		self.dbn = np.array([0+0j] * DIM)
		self.__a1 = np.array([0+0j] * DIM)
		self.__a2 = np.array([0+0j] * DIM)
		self.__a3 = np.array([0+0j] * DIM)
		self.__a4 = np.array([0+0j] * DIM)

	#１階微分を与えるメソッド
	def Db(self, t, bn, out_bn ):
		for n in range(DIM):
			out_bn[n] = Energy(n) / (I * hbar) * bn[n] 

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


#ルンゲクッタクラスのインスタンスを生成
rk4 = RK4(DIM, dt)
#初期状態の設定
rk4.bn = np.array(
	[ 1.0+0.0j,   #基底状態
	  0.0+0.0j,   #第１励起状態 
	  0.0+0.0j ]  #第２励起状態
)
#グラフ描画用配列の準備
ts = []
b0s = []

###展開係数の時間依存性の計算
for tn in range(Tn+1):
	#実時間
	t_real = dt * tn
	#計算結果の取得
	if( tn % skip == 0 ): 
		#時刻の取得
		ts.append( tn )
		#基底状態の展開係数の取得
		b0s.append( rk4.bn[0] )

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
plt.ylim([-1,1])
#グラフ描画
plt.plot(ts, np.real(b0s), marker="o" , linewidth=3.0)
plt.plot(ts, np.imag(b0s), marker="o",  linewidth=3.0)
#グラフ表示
plt.show()
