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
#２つの井戸の中心距離
R = 1.2E-9 #1.2nm
#計算区間
x_min = -L / 2.0
x_max = L / 2.0
#状態数
n_max = 3
#行列の要素数
DIM = n_max + 1

#クーロン相互作用ポテンシャル項[eV]
def V12(x1, x2):
	return 1
	#return e * e / (4.0 * math.pi * epsilon0)/( abs(x2 - x1) ) / eV

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

#被積分関数
def integral_V12(x2, x1, n1, n2, m1, m2, R, L):
	return varphi12(n1, n2, x1, x2, R, L) * V12(x1, x2) * varphi12(m1, m2, x1, x2, R, L)

#積分区間の設定
x1_min = -R / 2.0 - L / 2.0
x1_max = -R / 2.0 + L / 2.0
x2_min = R / 2.0 - L / 2.0
x2_max = R / 2.0 + L / 2.0

###<m1,m2|V12|n1,n2>の計算
for m1 in range(DIM):
	for m2 in range(DIM):
		for n1 in range(DIM):
			for n2 in range(DIM):
				#ガウス・ルジャンドル２重積分
				result = integrate.dblquad(
					integral_V12,                         #被積分関数
					x1_min, x1_max,                       #第１引数に対する積分区間の下端と上端
					lambda x : x2_min, lambda x : x2_max, #第２引数に対する積分区間の下端と上端
					args=(n1, n2, m1, m2, R, L)           #被積分関数へ渡す引数
				)
				V_real = result[0]
				#ターミナルへ出力
				print(f'{m1},{m2},{n1},{n2}, {V_real}')
