import math
import scipy.integrate as integrate

#積分区間
x_min = 0.0
x_max = 1.0

#被積分関数
def integrand(x):
	return math.sin(math.pi*x)
#解析値
exact = 2.0 / math.pi

#数値積分の実行
result = integrate.quad(integrand,x_min,x_max)

#計算結果を表示
print("積分結果：" + str(result[0]))
#計算誤差と推定誤差を表示
print("計算誤差：" + str(result[0]-exact) + " (推定誤差：" + str(result[1]) + ")")
