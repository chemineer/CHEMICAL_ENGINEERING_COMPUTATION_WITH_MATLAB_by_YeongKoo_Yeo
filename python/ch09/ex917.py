import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. 시스템 구성 요소 정의
s = ct.TransferFunction.s

# (1) PI Controller: Gc(s) = 2 * (1 + 1/(5s))
Gc = 2 * (1 + 1/(5*s))

# (2) Control Valve: Kv = 0.01
Kv = 0.01

# (3) Process: Gp(s) = 5 / ((s+1)(2s+1))
Gp = 5 / ((s + 1) * (2 * s + 1))

# (4) Sensor/Transducer: Km = 20
Km = 20

# 2. 피드백 제어 시스템 구성
# 전향 경로 (Forward Path): G = Gc * Kv * Gp
G_forward = Gc * Kv * Gp

# 폐루프 시스템 (Closed-loop): T = G_forward / (1 + G_forward * Km)
# ct.feedback(전향, 피드백) 함수를 사용합니다.
T = ct.feedback(G_forward, Km)

# 3. 단위 계단 응답(Unit Step Response) 계산
t = np.linspace(0, 100, 1000)  # 0초부터 100초까지 계산
time, response = ct.step_response(T, T=t)

# 4. 결과 시각화
plt.figure(figsize=(10, 6))
plt.plot(time, response, label='Step Response (Unit Step)')
plt.axhline(y=1/Km, color='r', linestyle='--', label=f'Steady State (1/Km = {1/Km})') # 최종값 확인용
plt.title('Step Response of a PI Feedback Control System (Example 9.17)')
plt.xlabel('Time (sec)')
plt.ylabel('Response y(t)')
plt.grid(True)
plt.legend()
plt.show()