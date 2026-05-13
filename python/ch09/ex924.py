import control as ct
import matplotlib.pyplot as plt
import numpy as np

w = np.logspace(-2, 2, 300)
# 1. 전달함수 정의
# 분자: 0.4s + 1
num = [0.4, 1] 

# 분모: (0.3s + 1)(s + 1)(s + 1)
# np.convolve를 사용하여 다항식 곱셈 수행
den = np.convolve([0.3, 1], np.convolve([1, 1], [1, 1]))

# 2. 시간 지연(ioDelay) 설정 및 시스템 생성
# 매트랩의 'iodelay', 0.2를 ioTimeDelay 속성으로 적용
G = ct.tf(num, den)
G.ioTimeDelay = 0.2

# 3. 보드 선도 데이터 계산
# G는 전달함수, w는 주파수 범위(omega)
response = ct.frequency_response(G, omega=w)

# 데이터 추출 (이 방식은 Plot 인자 에러를 완벽히 회피합니다)
mag = response.magnitude
phase = response.phase
w = response.omega

# 4. 결과 시각화
plt.figure(figsize=(10, 8))

# Magnitude Plot (상단)
plt.subplot(2, 1, 1)
# 매트랩의 loglog와 동일한 로그 스케일 그래프
plt.loglog(w, mag)
plt.grid(True, which="both")
plt.ylabel('Amplitude')
plt.xlabel('Frequency(rad/time)')
plt.title('Bode Plot with Time Delay (Example 9.24)')

# Phase Plot (하단)
plt.subplot(2, 1, 2)
# 매트랩의 semilogx와 동일하게 주파수축만 로그 스케일
# phase는 라디안 단위이므로 도(degree) 단위로 변환하여 출력
plt.semilogx(w, np.degrees(phase))
plt.grid(True, which="both")
plt.ylabel('Phase(deg)')
plt.xlabel('Frequency(rad/time)')

plt.tight_layout()
plt.show()