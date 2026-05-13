import control as ctrl
import numpy as np
import matplotlib.pyplot as plt

# 매개변수 설정
tau = 0.5
t = np.arange(0, 10.1, 0.1)
num = [1]

# 감쇠비에 따른 분모 정의
# 전달함수: 1 / (tau^2 * s^2 + 2 * tau * zeta * s + 1)
def get_den(zeta):
    return [tau**2, 2 * tau * zeta, 1]

# 시스템 객체 생성
G1 = ctrl.tf(num, get_den(0.5)) # Underdamped (zeta < 1)
G2 = ctrl.tf(num, get_den(1.0)) # Critically damped (zeta = 1)
G3 = ctrl.tf(num, get_den(1.5)) # Overdamped (zeta > 1)

# 계단 응답 계산
_, y1 = ctrl.step_response(G1, T=t)
_, y2 = ctrl.step_response(G2, T=t)
_, y3 = ctrl.step_response(G3, T=t)

# 그래프 그리기
plt.figure(figsize=(8, 5))
plt.plot(t, y1, ':', label='Underdamped')
plt.plot(t, y2, '--', label='Critically damped')
plt.plot(t, y3, '-', label='Overdamped')

plt.grid(True)
plt.xlabel('Time t(sec)')
plt.ylabel('Output y(t)')
plt.title('Step responses of a 2nd-order process')
plt.legend(loc='best')
plt.show()