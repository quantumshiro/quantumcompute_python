import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import matplotlib
print(matplotlib.matplotlib_fname())

#図全体
fig = plt.figure(figsize=(10, 6))

#描画範囲
x_min = -1.0
x_max = 1.0

#アニメーション作成用
ims=[]

#描画区間数
N = 100
#アニメーション分割数
AN = 30

#アニメーションの各グラフを生成
for a in range(AN):
	phi = 2.0 * math.pi * a / AN 
	#リストの生成
	xl = []
	yl = []
	#データの生成
	for i in range(N+1):
		x = x_min + (x_max -x_min)*i/N
		y = math.sin( math.pi * x + phi)
		xl.append(x)
		yl.append(y)

	#各コマのグラフの描画
	img  = plt.plot(xl, yl, color="blue", linewidth=3.0, linestyle="solid" )
	ims.append( img )

#グラフタイトルの設定 
plt.title("sin function")
#x軸ラベル・y軸ラベルの設定 
plt.xlabel("x-axis")
plt.ylabel("y-axis")
#描画範囲を設定
plt.xlim([-1.0,1.0])
plt.ylim([-1.0,1.0])
#アニメーションの生成
ani = animation.ArtistAnimation(fig, ims, interval=50)
#アニメーションの保存
#ani.save("output.html", writer=animation.HTMLWriter())
#ani.save("output.gif", writer="imagemagick")
#ani.save("output.mp4", writer="ffmpeg", dpi=300)
#ani.save("output.mp4", writer=matplotlib.animation.AVConvWriter, dpi=300)
#グラフの表示
plt.show()
