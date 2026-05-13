import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. 설정값 및 시간 범위 정의
Kc_list = [5, 20, 50]
t = np.arange(0, 10.1, 0.1) # [0:0.1:10]
styles = [':', '--', '-']   # 매트랩의 선 스타일 대응
labels = ['Kc = 5', 'Kc = 20', 'Kc = 50']

# 2. 반복문을 이용한 시스템 정의 및 계단 응답 계산
for Kc, style, label in zip(Kc_list, styles, labels):
    num = [Kc]
    den = [5, 1 + Kc]
    
    # 시스템 생성
    sys = ct.tf(num, den)
    
    # 계단 응답 계산 (지정된 시간 t 사용)
    time, response = ct.step_response(sys, T=t)
    
    # 그래프 그리기
    plt.plot(time, response, style, label=label)

# 3. 그래프 꾸미기
plt.xlabel('Time t(sec)')
plt.ylabel('output y (t)')
plt.legend()
plt.grid(True)
plt.title('Step Responses for different Kc')
plt.show()