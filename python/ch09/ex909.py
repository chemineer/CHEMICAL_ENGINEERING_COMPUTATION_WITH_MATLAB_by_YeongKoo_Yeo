import control as ctrl
import numpy as np
import matplotlib.pyplot as plt

# 1. 전달함수 정의 (G(s) = 3 / (2s + 1))
num = [3]
den = [2, 1]
G = ctrl.tf(num, den)

# 2. 시간 배열 생성 (0부터 10까지 0.1 간격)
t = np.arange(0, 10.1, 0.1)

# 3. 계단 응답 계산
# t를 명시적으로 전달하여 응답을 계산합니다.
t_out, y_out = ctrl.step_response(G, T=t)

# 4. 그래프 그리기
plt.figure(figsize=(8, 5))
plt.plot(t_out, y_out)
plt.grid(True)
plt.title('Step response of a 1st-order process')
plt.xlabel('t (sec)')
plt.ylabel('Response y(t)')
plt.show()