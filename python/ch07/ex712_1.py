import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# 1. 시스템 행렬 정의
A = np.array([
    [-0.325,  0.125,  0.   ,  0.   ,  0.   ],
    [ 0.2  , -0.325,  0.125,  0.   ,  0.   ],
    [ 0.   ,  0.2  , -0.325,  0.125,  0.   ],
    [ 0.   ,  0.   ,  0.2  , -0.325,  0.125],
    [ 0.   ,  0.   ,  0.   ,  0.2  , -0.325]
])
B = np.array([[0.2, 0], [0, 0], [0, 0], [0, 0], [0, 0.25]])
C = np.array([[0, 0, 0, 0, 1], [0.5, 0, 0, 0, 0]])
D = np.array([[0, 0], [0, 0]])
Us = np.array([[0.0], [0.1]])

# 2. 정상 상태(Steady-state) 계산
# xs = -inv(A) * B * Us
xs = np.linalg.solve(A, -B @ Us)
# ys = C * xs + D * Us
ys = (C @ xs) + (D @ Us)

# 3. 상태 공간 시스템 정의 및 계단 응답 계산
# MATLAB의 step(A,B,C,D,2)은 입력 1에 대한 응답을 기본으로 함
# 여기서는 0.05 크기의 변화(0.1 Us의 차이 등 반영)를 위해 수동 계산
sys = signal.StateSpace(A, B, C, D)
t = np.linspace(0, 80, 500) # 시간 범위 설정 (MATLAB의 2는 최종 시간의 비중을 의미)

# 입력 2(Us[1]=0.1)에 대한 응답 계산
# MATLAB 코드의 0.05*y는 입력 변화량을 반영한 것으로 보임
t, y_step, x_step = signal.lsim(sys, U=np.tile([0, 0.05], (len(t), 1)), T=t)

# 4. 전체 응답 계산 (정상 상태 + 변화량)
# Y(t) = ys + y(t), X(t) = xs + x(t)
y_final = y_step + ys.flatten()
x_final = x_step + xs.flatten()

# 5. 시각화 (Subplot)
plt.figure(figsize=(12, 10))

# Subplot 1: x5 응답
plt.subplot(2, 2, 1)
plt.plot(t, y_final[:, 0])
plt.xlabel('t(min)')
plt.ylabel('$x_5$')
plt.grid(True)

# Subplot 2: y1 응답
plt.subplot(2, 2, 2)
plt.plot(t, y_final[:, 1])
plt.xlabel('t(min)')
plt.ylabel('$y_1$')
plt.grid(True)

# Subplot 3: 모든 상태 변수 x1~x5 응답
plt.subplot(2, 2, 3)
styles = ['-', '--', '-.', '*-', 'o-']
labels = ['$x_1$', '$x_2$', '$x_3$', '$x_4$', '$x_5$']
for i in range(5):
    if i < 3: # 처음 3개는 스타일 적용
        plt.plot(t, x_final[:, i], styles[i], label=labels[i])
    else: # 나머지는 마커 간격 조절하여 출력
        plt.plot(t, x_final[:, i], styles[i], label=labels[i], markevery=20)

plt.xlabel('t(min)')
plt.ylabel('x')
plt.legend()
plt.tight_layout()
plt.show()