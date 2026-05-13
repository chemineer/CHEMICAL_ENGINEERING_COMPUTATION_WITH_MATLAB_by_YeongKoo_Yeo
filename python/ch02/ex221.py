import numpy as np
from lspolfit import lspolfit    
import matplotlib.pyplot as plt

# 데이터 정의
t = np.array([0, 2, 5, 6, 13, 20, 24, 30, 35, 41, 50])
c = np.array([0.86, 0.61, 0.47, 0.39, 0.25, 0.18, 0.15, 0.12, 0.10, 0.09, 0.08])

# 보간을 위한 시간축 설정
tint = np.linspace(t.min(), t.max(), 100)

# 설정: 선 스타일 및 다항식 차수
styles = ['--', ':', '-.']
orders = [3, 4, 5]
labels = ['3rd-order', '4th-order', '5th-order'] # MATLAB 주석 기준 (m+1 적용 시)

plt.figure(figsize=(10, 6))
plt.plot(t, c, 'o', label='Data', markersize=8)

# 각 차수별로 다항식 적합 및 시각화
for i, m in enumerate(orders):
    p = lspolfit(t, c, m)
    # numpy.polyval을 사용하여 다항식 값 계산
    cint = np.polyval(p, tint)
    plt.plot(tint, cint, styles[i], label=labels[i])

# 그래프 꾸미기
plt.xlabel('t(min)')
plt.ylabel('C_A(mol/liter)')
plt.legend(loc='best')
plt.title('Polynomial Fitting to Reaction Data')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()