############################################################################
#　１量子ビット（改良版量子井戸）の Xnm
############################################################################
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate
import numpy as np
import numpy.linalg as LA

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
n_max = 20
#行列の要素数
DIM = n_max + 1
#空間分割数
NX = 500
#空間刻み間隔
dx = 1.0 * 10**-9
#壁の厚さ
W = L / 5
#壁の高さの最大値
V_max = 30.0 * eV
#電場の強さ
Ex = 1.0 * 10**8
#基底状態と励起状態
N = 2

#固有関数
def verphi(n, x):
	kn = math.pi * (n + 1) / L
	return math.sqrt(2.0 / L) * math.sin(kn * (x + L / 2.0))

#ポテンシャル項
def V(x, Ex):
	if(abs(x) <= W / 2.0): 
		return V_max + e * Ex * x
	else: 
		return e * Ex * x

#固有エネルギー
def Energy(n):
	kn = math.pi * (n + 1) / L
	return hbar * hbar * kn * kn / (2.0 * me)

#静電場が加わったときの固有状態
def Qbit_verphi(n, x, an, n_max):
	phi = 0
	#全状態の重ね合わせ
	for m  in range(n_max + 1):
		phi += an[n][m] * verphi(m, x)

	return phi

#被積分関数（行列要素計算用）
def integral_matrixElement(x, n1, n2, Ex):
	return verphi(n1 ,x) * V(x, Ex) * verphi(n2, x) / eV

#被積分関数（Xmn計算用）
def integral_Xnm(x, n1, n2, n_max, an):
    return Qbit_verphi(n1, x, an, n_max) * x * Qbit_verphi(n2, x, an, n_max)

#被積分関数（平均計算用）
def average_x(x, a):
	sum = 0
	for n in range(n_max + 1):
		sum += a[n] * verphi(n, x)
	return x * sum**2

#固有値・固有ベクトルの初期化
eigenvalues = []
vectors = []

#存在確率分布グラフ描画用の配列初期化
xs = []
phi = [0] * N
for n in range(N):
	phi[n] = [0] * (NX + 1) 

#中心の電場依存性グラフ描画用の配列初期化
averageX = [0] * N

#エルミート行列（リスト）
matrix = []

######################################
#　計算開始
######################################

###行列要素の計算
for n1 in range(n_max + 1):
	col=[]
	for n2 in range(n_max + 1):
		#ガウス・ルジャンドル積分
		result = integrate.quad(
			integral_matrixElement, #被積分関数
			x_min, x_max,		        #積分区間の下端と上端
			args=(n1, n2, Ex)			  #被積分関数へ渡す引数
		)
		real = result[0]
		imag = 0j
		#無静電場のエネルギー固有値（対角成分）
		En = Energy(n1)/eV if (n1 == n2) else 0
		#行の要素を追加
		col.append( En + real )
	#行を追加
	matrix.append( col )

#リスト → 行列
matrix = np.array( matrix )

###固有値と固有ベクトルの計算
result = LA.eig( matrix )
eig = result[0] #固有値
vec = result[1] #固有ベクトル

#小さい順に並べるためのインデックス（配列）
index = np.argsort( eig )
#固有値を小さい順に並び替え
eigenvalues = eig[ index ]

#転置行列
vec = vec.T
#固有ベクトルの並び替え
vectors = vec[ index ]

### 検算：MA-EA=0 ?
sum = 0
for i in range(DIM):
	v = matrix @ vectors[i] - eigenvalues[i] * vectors[i]
	for j in range(DIM):
		sum += abs(v[j])**2
print("|MA-EA| =" + str(sum))

###固有関数の空間分布
for nx in range(NX+1):
	x = x_min + (x_max - x_min) / NX * nx
	xs.append( x / dx )
	for n in range( len(phi) ):
		#描画用データの整形
		phi[n][nx] = abs(Qbit_verphi(n, x, vectors , n_max))**2 / (1.0 * 10**9)

for n in range(len(averageX)):
	#ガウス・ルジャンドル積分
	result = integrate.quad(
		average_x,            #被積分関数
		x_min, x_max,		  #積分区間の下端と上端
		args=(vectors[n])	  #被積分関数へ渡す引数
	)
	#計算結果の取得
	averageX[n] = result[0] * (1.0 * 10**9)

###行列要素（Xnm）の計算
for n1 in range( N ):
	for n2 in range( N ):
		#ガウス・ルジャンドル積分
		result = integrate.quad(
			integral_Xnm,                     #被積分関数
			x_min, x_max,		              #積分区間の下端と上端
			args=(n1, n2, n_max, vectors)  #被積分関数へ渡す引数
		)
		real = result[0]
		imag = 0
		if( abs(real / L) < L ): real = 0
		#ターミナルへ出力
		print( "(" + str(n1) + ", " + str(n2) + ")  " + str( real / L ))



#グラフの描画（エネルギー固有値）
fig1 = plt.figure(figsize=(10, 6))
plt.title("Energy eigenvalue")
plt.xlabel("Level")
plt.ylabel("Energy[eV]")
#描画範囲を設定
plt.xlim([0, n_max])
#x軸
exs = range( n_max + 1)
#y軸
En = []
#固有エネルギーの描画
plt.plot(exs, eigenvalues, marker="o", linewidth = 3)

#グラフの描画（基底状態と励起状態）
fig2 = plt.figure(figsize=(10, 6))
plt.title("Existence probability at Position")
plt.xlabel("Position[nm]")
plt.ylabel("|phi|^2")
#描画範囲を設定
plt.xlim([-0.5, 0.5])
plt.ylim([0, 5.0])
plt.plot(xs, phi[0] , linewidth = 3)	
plt.plot(xs, phi[1] , linewidth = 3)	

'''
#グラフの描画（期待値）
fig4 = plt.figure(figsize=(10, 6))
plt.title("Position at Electric field strength")
plt.xlabel("Electric field strength[V/m]")
plt.ylabel("Position[nm]")
#描画範囲を設定
plt.xlim([0, 10])
#x軸
exs = range( NV + 1)
plt.plot(exs, averageX[0], marker="o", linewidth = 3)
plt.plot(exs, averageX[1], marker="o", linewidth = 3)
#グラフの表示
'''
plt.show()

