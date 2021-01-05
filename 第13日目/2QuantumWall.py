############################################################################
#　クーロン相互作用のみを考慮したときの固有状態
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
#光速
c = 2.99792458E+8
#真空の透磁率
mu0 = 4.0 * math.pi * 1.0E-7
#真空の誘電率
epsilon0 = 1.0 / (4.0 * math.pi * c * c) * 1.0E+7

######################################
#　物理系の設定
######################################
#量子井戸の幅
L = 1.0 * 10**-9
#二つの井戸の中心距離
R_min = L + 0.1 * 10**-9
R_max = L + 2.0 * 10**-9
#距離分割数
NR = 19
#計算区間
x_min = -L / 2.0
x_max = L / 2.0
#空間分割数
NX = 500
#空間刻み間隔
dx = 1.0 * 10**-9
#状態数
n_max = 5
#行列の要素数
DIM = n_max + 1
#基底状態・第１励起状態・第２励起状態・第３励起状態
N = 4

#クーロン相互作用ポテンシャル項[eV]
def V12(x1, x2):
	return e * e / (4.0 * math.pi * epsilon0)/( abs(x2 - x1) ) / eV

#独立量子井戸の固有関数
def varphi1(n1, x, R, L):
	kn = math.pi * (n1 + 1) / L
	return math.sqrt(2.0 / L) * math.sin( kn * (x + L / 2.0 + R / 2.0 ))

#独立量子井戸の固有関数
def varphi2(n2, x, R, L):
	kn = math.pi * (n2 + 1) / L
	return math.sqrt(2.0 / L) * math.sin(kn * (x + L / 2.0 - R / 2.0))

#独立量子井戸の固有関数の積
def varphi12(n1, n2, x1, x2, R, L):
	return varphi1(n1, x1, R, L) * varphi2(n2, x2, R, L)

#独立量子井戸の固有エネルギー[eV]
def Energy0_12(n1, n2, L):
	kn1 = math.pi * (n1 + 1) / L
	kn2 = math.pi * (n2 + 1) / L
	En1 = hbar**2 * kn1**2 / (2.0 * me)
	En2 = hbar**2 * kn2**2 / (2.0 * me)
	return (En1 + En2) /eV

#2量子ビットの固有状態
def Qbit_varphi(x1, x2, an, n_max, R, L):
	phi=0
	for m1 in range(0,DIM):
		for m2 in range(0,DIM):
			m = m1 * ( n_max + 1 ) + m2
			phi += an[m] * varphi12(m1, m2, x1, x2, R, L)
	return phi.real

#確率密度分布
def rho(x, an, n_max, R, L):
	if(x <= - R / 2.0 + L / 2.0):
		x_min = R / 2 - L / 2
		x_max = R / 2 + L / 2
		#ガウス・ルジャンドル積分
		result = integrate.quad(
			integral_varphi_x2,    #被積分関数
			x_min,x_max,           #積分区間の下端と上端
			args=(n_max,x,R,L,an)  #被積分関数へ渡す引数
		)
	elif( R / 2.0 - L / 2.0 <= x):
		x_min = - R / 2 - L / 2
		x_max = - R / 2 + L / 2
		#ガウス・ルジャンドル積分
		result = integrate.quad(
			integral_varphi_x1,    #被積分関数
			x_min,x_max,           #積分区間の下端と上端
			args=(n_max,x,R,L,an)  #被積分関数へ渡す引数
		)
	real = result[0]
	return real

#被積分関数
def integral_V12(x2, x1, n1, n2, m1, m2, R, L):
	return varphi12(n1, n2, x1, x2, R, L) * V12(x1, x2) * varphi12(m1, m2, x1, x2, R, L)

#密度分布計算用の被積分関数
def integral_varphi_x1(x1, n_max, x, R, L, an):
	varphi = Qbit_varphi(x1, x, an,n_max, R, L)
	return abs(varphi)**2

#密度分布計算用の被積分関数
def integral_varphi_x2(x2, n_max, x, R, L, an):
	varphi = Qbit_varphi(x, x2, an, n_max, R, L)
	return abs(varphi)**2

#固有値・固有ベクトルの初期化
eigenvalues = [0] * (NR + 1)
vectors = [0] * (NR + 1)
for nR in range(NR + 1):
	eigenvalues[nR] = []
	vectors[nR] = []

#存在確率分布グラフ描画用の配列初期化
xs = []
phi = [0] * (NR + 1)
for nR in range(NR + 1):
	phi[nR] = [0] * N
	for n in range( len(phi[nR]) ):
		phi[nR][n] = [0] * (NX + 1)

######################################
#　計算開始
######################################
DIM1 = n_max + 1
DIM2 = DIM1 * DIM1

###Rの大きさの変更
for nR in range( NR + 1 ):

	#Rの設定
	if(NR != 0):
		print("井戸間隔：" + str( nR * 100 / NR ) + "%")
		R = R_min + (R_max - R_min) * nR / NR
	else:
		R = 2.0 * L

	#積分区間の設定
	x1_min = -R / 2.0 - L / 2.0
	x1_max = -R / 2.0 + L / 2.0
	x2_min = R / 2.0 - L / 2.0
	x2_max = R / 2.0 + L / 2.0

	#エルミート行列の設定
	matrix = [[0]*DIM2 for i in range(DIM2)]

	for m1 in range(DIM1):
		for m2 in range(DIM1):
			for n1 in range(DIM1):
				for n2 in range(DIM1):
					#ガウス・ルジャンドル２重積分
					result = integrate.dblquad(
						integral_V12,   #被積分関数
						x1_min, x1_max, #積分区間の下端と上端
						lambda x : x2_min,
						lambda x : x2_max,
						args=(n1, n2, m1, m2, R, L) #被積分関数へ渡す引数
					)
					V_real = result[0]
					V_imag = 0j
					#ターミナルへ出力
					#print(f'{m1},{m2},{n1},{n2}, {V_real}')
					#行列要素
					nn = n1 * DIM1 + n2
					mm = m1 * DIM1 + m2
					matrix[mm][nn] = V_real + V_imag
					#対角成分
					if(nn == mm): matrix[mm][nn] += Energy0_12(n1, n2, L)

	#リスト → 行列
	matrix = np.array( matrix )
	###固有値と固有ベクトルの計算
	result = LA.eig( matrix )
	eig = result[0] #固有値
	vec = result[1] #固有ベクトル

	#小さい順に並べるためのインデックス（配列）
	index = np.argsort( eig )
	#固有値を小さい順に並び替え
	eigenvalues[nR] = eig[ index ]

	#転置行列
	vec = vec.T
	#固有ベクトルの並び替え
	vectors[nR] = vec[ index ]

	### 検算：MA-EA=0 ?
	sum = 0
	for i in range(DIM):
		v = matrix @ vectors[nR][i] - eigenvalues[nR][i] * vectors[nR][i]
		for j in range(DIM):
			sum += abs(v[j])**2
	print("|MA-EA| =" + str(sum))

	#基底状態・第１励起状態・第２励起状態・第３励起状態の存在確率の空間分布
	for n in range( N ):

		x_min = -R / 2 - L / 2
		x_max = R / 2 + L / 2

		for nx in range(NX+1):
			x = x_min + (x_max - x_min) / NX * nx
			if( n == 0 ): xs.append(x/dx)
			if( x <= - R / 2.0 + L / 2.0 ):
				phi[nR][n][nx] = rho(x, vectors[nR][n], n_max, R, L) * 1E-9
			elif( R / 2.0 - L / 2.0 <= x ):
				phi[nR][n][nx] = rho(x, vectors[nR][n], n_max, R, L) * 1E-9
			else: 
				phi[nR][n][nx] = 0

if(NR > 0):
	#グラフの描画（エネルギー固有値）
	fig1 = plt.figure(figsize=(10, 6))
	plt.title("Energy at Electric field strength")
	plt.xlabel("Electric field strength[V/m]")
	plt.ylabel("Energy[eV]")
	#描画範囲を設定
	plt.xlim([1, 3])
	#x軸
	exs = []
	for nR in range(NR + 1):
		R = R_min + (R_max - R_min) * nR / NR
		exs.append( R / dx )

	for n in range(N):
		#y軸
		En = []
		for nR in range(NR + 1):
			En.append( eigenvalues[nR][n] )

		#基底状態と第１励起状態と第２励起状態と第３励起状態のグラフを描画	
		plt.plot(exs, En, marker="o", linewidth = 3)

else:
	for n in range(N):
		#グラフの描画（基底状態～第３励起状態）
		fig2 = plt.figure(figsize=(10, 6))
		plt.title("Existence probability at Position (n=" + str(n) + ")")
		plt.xlabel("Position[nm]")
		plt.ylabel("|phi|^2")
		#描画範囲を設定
		plt.xlim([-1.5, 1.5])
		plt.ylim([0, 2.2])
		nR = 0 
		plt.plot(xs, phi[nR][n], linewidth = 3)

plt.show()

