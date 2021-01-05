############################################################################
#　「エネルギー固有関数の直交性」の確認
############################################################################
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate

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

######################################
#　物理系の設定
######################################
#量子井戸の幅
L = 1 * 10**-9
#計算区間
x_min = -L / 2.0
x_max = L / 2.0
#状態数
n_max = 5

#固有関数
def verphi(n, x):
	kn = math.pi * (n + 1) / L
	return math.sqrt(2.0 / L) * math.sin(kn * (x + L / 2.0))

#固有エネルギー
def Energy(n):
	kn = math.pi * (n + 1) / L
	return hbar * hbar * kn**2 / (2.0 * me)

#被積分関数
def integral_orthonomality(x, n1, n2):
	return verphi(n1,x) * verphi(n2,x) 

###エネルギー固有関数の直行性の確認
for n1 in range(n_max + 1):
	for n2 in range(n_max + 1):
		#ガウス・ルジャンドル積分
		result = integrate.quad(
			integral_orthonomality, #被積分関数
			x_min, x_max,		   #積分区間の下端と上端
			args=(n1,n2)			#被積分関数へ渡す引数
		)
		#ターミナルへ出力
		print( "(" + str(n1) + ", " + str(n2) + ")  " + str(result[0]))
