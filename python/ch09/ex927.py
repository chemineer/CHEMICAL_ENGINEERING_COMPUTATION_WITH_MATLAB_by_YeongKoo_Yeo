import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. 시스템 전달함수 정의
# 분자: 2
num = [2]
# 분모: (10s + 1)(2.5s + 1) -> np.convolve로 다항식 곱셈 수행
den = np.convolve([10, 1], [2.5, 1])
sys = ct.tf(num, den)

# 2. 주파수 범위 설정 (10^-2 ~ 10^1 rad/time, 300개 지점)
# 매트랩: w = logspace(-2, 1, 300);
w = np.logspace(-2, 1, 300)

# 3. 니콜스 선도(Nichols Chart) 그리기
plt.figure(figsize=(8, 8))

# ct.nichols_plot 함수를 사용하여 선도를 생성합니다.
# grid=True는 매트랩의 ngrid 명령과 유사하게 등이득/등위상 궤적을 표시합니다.
ct.nichols_plot(sys, omega=w, grid=True)

# 4. 그래프 제목 설정
plt.title('Nichols Plot of a 2nd-order process')

# 그래프 출력
plt.show()