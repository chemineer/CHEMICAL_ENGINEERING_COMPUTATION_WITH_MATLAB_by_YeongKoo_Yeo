import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# 시스템 정의
num = [3]
den = [2, 1]
system = signal.TransferFunction(num, den)

# 시간 및 입력 정의
t = np.arange(0, 10.1, 0.1)
u = np.sin(3 * t)

# scipy.signal.lsim 사용 (이 함수는 전통적으로 lsim을 지원합니다)
t_out, y_out, x_out = signal.lsim(system, U=u, T=t)

# 4. 그래프 그리기
plt.figure(figsize=(10, 5))
plt.plot(t_out, y_out, label='Output y(t)')
plt.plot(t_out, u, ':', label='Input u(t)')
plt.axhline(0, color='black', linewidth=0.5)  # t*0 라인 (x축)
plt.legend()
plt.title('Sinusoidal response of 1st-order process')
plt.xlabel('t (sec)')
plt.ylabel('Response y(t)')
plt.grid(True)
plt.show()