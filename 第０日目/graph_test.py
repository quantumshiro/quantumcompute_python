import math
import matplotlib.pyplot as plt
#描画範囲
x_min = -1.0
x_max = 1.0
#描画区間数
N = 100
#リストの生成
xl = []
yl = []
#データの生成
for i in range(N+1):
	x = x_min + (x_max -x_min)*i/N
	y = math.sin( math.pi * x )
	xl.append(x)
	yl.append(y)

#描画範囲を設定
plt.xlim([-1.0,1.0])
plt.ylim([-1.0,1.0])
#グラフの描画
plt.plot(xl,yl)
#グラフの表示
plt.show()

