############################################################################
#　電子波束アニメーション
############################################################################
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#図全体
fig1 = plt.figure(figsize=(10, 6))
#全体設定
plt.rcParams['font.family'] = 'Times New Roman' #フォント
plt.rcParams['font.size'] = 12 #フォントサイズ

#複素数
I=0.0+1.0j

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
#空間分割数
NX = 500
#空間分割サイズ
dx = 1.0 * 10**-9
#計算区間
x_min = -10.0 * dx
x_max = 10.0 * dx
#重ね合わせの数
NK = 200
#パルスの幅
delta_x = 2.0 * 10**-9
sigma = 2.0 * math.sqrt(math.log(2.0)) / delta_x
#波数の間隔
dk = 20.0 / (delta_x * NK)
#波束の中心エネルギー
E0 = 10.0 * eV
#波束の中心
k0 = math.sqrt( 2.0 * me * E0 / hbar**2)
omega0 = hbar / (2.0*me) * k0**2
#計算時間の幅
ts = -50
te =  50
#時間間隔
dt = 1.0 * 10**-16

#アニメーション作成用
ims=[]

### 波束の空間分布の計算
for tn in range(ts,te):
	#実時間の取得
	t_real = dt*tn
	#データ格納用
	xl = []
	Psi_real = []
	Psi_imag = []
	Psi_abs = []
	for nx in range(NX):
		#空間地点の座標を取得
		x = x_min + (x_max - x_min) * nx / NX
		#波動関数値の初期化
		Psi = 0.0 + 0.0j
		#各波数ごとの寄与を足し合わせる
		for kn in range(NK):
			#各波数を取得
			k = k0 + dk * (kn - NK/2)
			#波数から各振動数を取得
			omega = hbar / (2.0 * me) * k**2
			#平面波を足し合わせる
			Psi += cmath.exp(I * (k * x - omega * t_real)) * cmath.exp(-((k - k0) / (2.0 * sigma))**2 )

		#波動関数の値がNK,sigmaの値に依存しないように規格化
		Psi = Psi * dk * dx /10.0

		#計算結果を表示
		#print(f"x/dx:{x/dx}\nPsi_real:{Psi.real}\nPsi_imag:{Psi.imag}\nPsi_abs={abs(Psi)}")

		#描画用データ生成
		xl.append( x / dx )
		Psi_real.append( Psi.real )
		Psi_imag.append( Psi.imag )
		Psi_abs.append( abs(Psi) )

	#各コマを描画
	img  = plt.plot(xl, Psi_real, 'red')
	img += plt.plot(xl, Psi_imag, 'green')
	img += plt.plot(xl, Psi_abs,  'blue')
	time = plt.text(0.9, 1.03, 't={:.2e}'.format(t_real),transform=plt.gca().transAxes, ha='center', va='center')
	#テキストをグラフに追加
	img.append(time)
	#アニメーションに追加
	ims.append(img)

#アニメーション作成用
kl=[]
a_kl=[]

### 波数分布の出力
for kn in range(NK):
	k = k0 + dk * (kn - NK / 2)
	a_k = math.exp(-((k-k0) / (2.0 * sigma))**2)
	kl.append( k / dk)
	a_kl.append( a_k )

### グラフ描画
plt.title("Gaussian wave packet(Spatial distribution)")
plt.xlabel("Position[nm]", fontsize = 16)
plt.ylabel("Probability amplitude", fontsize = 16)
#描画範囲を設定
plt.xlim([-10.0,10.0])
plt.ylim([-0.3,0.3])
#アニメーションの生成
ani = animation.ArtistAnimation(fig1, ims, interval=50)
#アニメーションの保存
ani.save("output.html", writer=animation.HTMLWriter())

#図全体
fig2=plt.figure(figsize=(5, 5))
### グラフ描画
plt.title("Gaussian wave packet(Wave number distribution)")
plt.xlabel("Wave number", fontsize = 16)
plt.ylabel("distribushon", fontsize = 16)
plt.plot(kl,a_kl)

#レイアウトの調整
plt.tight_layout()
#グラフの表示
plt.show()
