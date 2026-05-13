import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. 전달함수 정의: G(s) = 12.76 / (5s + 1)
num = [12.76]
den = [5, 1]
G = ct.tf(num, den)

# 2. 시간 지연(iodelay) 설정: 1초
# 파이썬 control 라이브러리에서는 ioTimeDelay 속성을 사용합니다.
G.ioTimeDelay = 1

# 3. 나이퀴스트 선도 그리기
plt.figure(figsize=(8, 8))

# ct.nyquist_plot 함수를 사용하여 나이퀴스트 선도를 생성합니다.
# arrows=True 옵션으로 경로의 방향을 표시할 수 있습니다.
ct.nyquist_plot(G)

# 4. 그래프 설정 (매트랩의 xlabel, ylabel, title 대응)
plt.xlabel('Real axis')
plt.ylabel('Imaginary axis')
plt.title('Nyquist diagram of a time delay')
plt.grid(True)

# 축 범위 자동 조정 및 출력
plt.tight_layout()
plt.show()