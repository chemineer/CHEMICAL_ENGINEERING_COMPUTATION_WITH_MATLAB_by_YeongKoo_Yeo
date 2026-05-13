import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. 주파수 범위 설정 (10^-2 ~ 10^1 rad/time, 300개 지점)
# 매트랩: w = logspace(-2, 1, 300);
w = np.logspace(-2, 1, 300)

# 2. 전달함수 정의
# 분자: 2
num = [2]
# 분모: (10s + 1)(2.5s + 1) -> np.convolve로 다항식 곱셈 수행
den = np.convolve([10, 1], [2.5, 1])
sys = ct.tf(num, den)

# 3. 주파수 응답 데이터 계산 (Bode 데이터 추출)
# 최신 버전에서는 plot=False(소문자)를 사용하거나 frequency_response를 권장합니다.
response = ct.frequency_response(sys, omega=w)
mag = response.magnitude
phase = response.phase  # 단위: 라디안(radian)

# 4. 결과 시각화
plt.figure(figsize=(12, 5))

# (좌측) 폴라 선도 - Polar Plot
plt.subplot(1, 2, 1, projection='polar')
# 매트랩의 polar((pi/180)*p, x)와 동일하게 위상(라디안)과 진폭을 전달합니다.
plt.plot(phase, mag)
plt.title('polar plot of a 2nd-order process')

# (우측) 나이퀴스트 선도 - Nyquist Plot
plt.subplot(1, 2, 2)
# ct.nyquist_plot은 나이퀴스트 선도를 자동으로 그려줍니다.
ct.nyquist_plot(sys, omega=w)
plt.title('Nyquist Plot')
plt.grid(True)

plt.tight_layout()
plt.show()