############################################################################
#　無限に深い量子井戸中の電子に静電場を加えたときの固有状態（ステップ２）
############################################################################
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate
import numpy as np
import numpy.linalg as LA

#図全体
fig1 = plt.figure(figsize=(10, 6))
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
#電場の強さ
Ex = 1.0*10**10

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

#被積分関数
def integral_matrixElement(x, n1, n2, Ex):
	return verphi(n1 ,x) * V(x, Ex) * verphi(n2, x) / eV

#エルミート行列（リスト）
matrix=[]
#行列の要素数
DIM = n_max + 1

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

#ターミナルへ出力
for i in range(DIM):
	print(f'{i}番目の固有値：{eigenvalues[i]}')
	print(f'{i}番目の固有値に対応する固有ベクトル：\n{vectors[i]}')

### 検算：MA-EA=0 ?
sum = 0
for i in range(DIM):
	v = matrix@vectors[i]-eigenvalues[i]*vectors[i]
	for j in range(DIM):
		sum += abs(v[j])**2
print("|MA-EA| =" + str(sum))
