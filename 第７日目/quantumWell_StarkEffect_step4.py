############################################################################
#　無限に深い量子井戸中の電子に静電場を加えたときの固有状態（ステップ４）
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
n_max = 10
#行列の要素数
DIM = n_max + 1
#空間分割数
NX = 500
#空間刻み間隔
dx = 1.0 * 10**-9
#電場の強さ
Ex_max = 1.0 * 10**10
#電場強度の分割数
NEx = 10

#固有関数
def verphi(n, x):
	kn = math.pi * (n + 1) / L
	return math.sqrt(2.0 / L) * math.sin(kn * (x + L / 2.0))

#ポテンシャル項
def V(x, Ex):
	return(e * Ex * x)
	
#固有エネルギー
def Energy(n):
	kn = math.pi * (n + 1) / L
	return hbar * hbar * kn * kn / (2.0 * me)

#被積分関数（行列要素計算用）
def integral_matrixElement(x, n1, n2, Ex):
	return verphi(n1 ,x) * V(x, Ex) * verphi(n2, x) / eV

#被積分関数（平均計算用）
def average_x(x, a):
	sum = 0
	for n in range(n_max + 1):
		sum += a[n] * verphi(n, x)
	return x * sum**2

#固有値・固有ベクトルの初期化
eigenvalues = [0] * (NEx + 1)
vectors = [0] * (NEx + 1)
for nEx in range(NEx + 1):
	eigenvalues[nEx] = []
	vectors[nEx] = []

#存在確率分布グラフ描画用の配列初期化
xs = []
phi = [0] * (NEx + 1)
for nEx in range(NEx + 1):
	phi[nEx] = [0] * 2
	for n in range( len(phi[nEx]) ):
		phi[nEx][n] = [0] * (NX + 1)

#中心の電場依存性グラフ描画用の配列初期化
averageX = [0] * 2
for n in range(len(averageX)):
	averageX[n] = [0] * (NEx + 1)


#静電場強度ごとに
for nEx in range(NEx + 1):
	print("電場強度：" + str( nEx * 100 / NEx ) + "%")
	#静電場の強度を設定
	Ex = Ex_max / NEx * nEx
	#エルミート行列（リスト）
	matrix=[]

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
	eigenvalues[nEx] = eig[ index ]

	#転置行列
	vec = vec.T
	#固有ベクトルの並び替え
	vectors[nEx] = vec[ index ]

	### 検算：MA-EA=0 ?
	sum = 0
	for i in range(DIM):
		v = matrix @ vectors[nEx][i] - eigenvalues[nEx][i] * vectors[nEx][i]
		for j in range(DIM):
			sum += abs(v[j])**2
	print("|MA-EA| =" + str(sum))

	###固有関数の空間分布
	for nx in range(NX+1):
		x = x_min + (x_max - x_min) / NX * nx
		if(nEx == 0): xs.append( x/dx )

		for n in range( len(phi[nEx]) ):
			for m in range(n_max+1):
				phi[nEx][n][nx] += vectors[nEx][n][m] * verphi(m, x)

			#描画用データの整形
			phi[nEx][n][nx] = abs(phi[nEx][n][nx])**2 / (1.0 * 10**9)


	for n in range(len(averageX)):
		#ガウス・ルジャンドル積分
		result = integrate.quad(
			average_x,            #被積分関数
			x_min, x_max,		      #積分区間の下端と上端
			args=(vectors[nEx][n])		#被積分関数へ渡す引数
		)
		#計算結果の取得
		averageX[n][nEx] = result[0] * (1.0 * 10**9)


#グラフの描画（エネルギー固有値）
fig1 = plt.figure(figsize=(10, 6))
plt.title("Energy at Electric field strength")
plt.xlabel("Electric field strength[V/m]")
plt.ylabel("Energy[eV]")
#描画範囲を設定
plt.xlim([0, 10])
#x軸
exs = range( NEx + 1)
#y軸
En_0 = []
En_1 = []
for nEx in range(NEx + 1):
	En_0.append( eigenvalues[nEx][0] )
	En_1.append( eigenvalues[nEx][1] )
#基底状態と第１励起状態のグラフを描画	
plt.plot(exs, En_0, marker="o", linewidth = 3)
plt.plot(exs, En_1, marker="o", linewidth = 3)

#グラフの描画（基底状態）
fig2 = plt.figure(figsize=(10, 6))
plt.title("Existence probability at Position (n=0)")
plt.xlabel("Position[nm]")
plt.ylabel("|phi|^2")
#描画範囲を設定
plt.xlim([-0.5, 0.5])
plt.ylim([0, 4.0])
#各
for nEx in range(NEx + 1):
	plt.plot(xs, phi[nEx][0] , linewidth = 3)	

#グラフの描画（第１励起状態）
fig3 = plt.figure(figsize=(10, 6))
plt.title("Existence probability at Position (n=1)")
plt.xlabel("Position[nm]")
plt.ylabel("|phi|^2")
#描画範囲を設定
plt.xlim([-0.5, 0.5])
plt.ylim([0, 3.0])
for nEx in range(NEx + 1):
	plt.plot(xs, phi[nEx][1] , linewidth = 3)	

#グラフの描画（期待値）
fig4 = plt.figure(figsize=(10, 6))
plt.title("Position at Electric field strength")
plt.xlabel("Electric field strength[V/m]")
plt.ylabel("Position[nm]")
#描画範囲を設定
plt.xlim([0, 10])
#x軸
exs = range( NEx + 1)
plt.plot(exs, averageX[0], marker="o", linewidth = 3)
plt.plot(exs, averageX[1], marker="o", linewidth = 3)
#グラフの表示
plt.show()

