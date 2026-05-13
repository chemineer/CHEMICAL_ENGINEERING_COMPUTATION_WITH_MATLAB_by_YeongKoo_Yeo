import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. G1 정의 (Time Delay 포함)
num1 = [3]
den1 = [3, 1]
G1 = ct.tf(num1, den1)
# 매트랩의 iodelay는 파이썬에서 
G1.ioTimeDelay = 1.6

# 2. G2 정의 (다항식 곱셈은 np.polyadd나 np.convolve 사용)
# conv([0.1 1], [0.5 1]) 등은 넘파이의 convolve와 동일합니다.
p1 = [0.1, 1]
p2 = [0.5, 1]
p3 = [1, 1]
p4 = [3, 1]

den2 = np.convolve(np.convolve(p1, p2), np.convolve(p3, p4))
G2 = ct.tf([3], den2)

# 3. Step Response 계산 및 시각화
t, y1 = ct.step_response(G1)
t2, y2 = ct.step_response(G2)

plt.plot(t, y1, label='G1')
plt.plot(t2, y2, ':', label='G2') # G2는 점선(':')으로 표시

# 4. 그래프 설정 (Legend, Grid 등)
plt.legend()
plt.grid(True)
plt.title('Step Response')
plt.xlabel('Time (sec)')
plt.ylabel('Amplitude')
plt.show()