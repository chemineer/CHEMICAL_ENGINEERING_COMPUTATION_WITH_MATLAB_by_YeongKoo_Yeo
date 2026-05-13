import numpy as np
import matplotlib.pyplot as plt
import control as ct

# 1. 시간축 설정 (0초부터 20초까지 0.1초 간격)
t = np.arange(0, 20.1, 0.1)

# 2. 기본 전달함수 정의: G(s) = 1 / (2s + 1)
# MATLAB의 [2 1]은 2s + 1을 의미합니다.
base_num = [1]
base_den = [2, 1]
G = ct.TransferFunction(base_num, base_den)

# 3. 차수별 시스템 계산 (MATLAB의 conv는 전달함수의 곱과 같습니다)
# n=2, n=4, n=5 시스템 생성
sys_n2 = G * G
sys_n4 = G * G * G * G
sys_n5 = G * G * G * G * G

# 4. 스텝 응답(Step Response) 계산
t1, y1 = ct.step_response(sys_n2, t)
t2, y2 = ct.step_response(sys_n4, t)
t3, y3 = ct.step_response(sys_n5, t)

# 5. 그래프 출력
plt.figure(figsize=(10, 6))
plt.plot(t1, y1, ':', label='n = 2')
plt.plot(t2, y2, '--', label='n = 4')
plt.plot(t3, y3, '-', label='n = 5')

plt.grid(True)
plt.xlabel('Time t(sec)')
plt.ylabel('Response y(t)')
plt.title('Step response of higher-order process (n=2, 4, 5)')
plt.legend(loc='best')
plt.show()