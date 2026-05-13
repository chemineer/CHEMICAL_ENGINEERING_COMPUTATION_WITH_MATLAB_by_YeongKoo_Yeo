import numpy as np
import matplotlib.pyplot as plt
from ngasZ import ngasZ

# 설정
T = 60
P = np.arange(100, 5010, 10) # 100부터 5000까지 10 간격
nz = []

# 데이터 계산
for k in range(4):
    Sg = 0.5 + (k * 0.1)
    nzv = ngasZ(T, P, Sg)
    nz.append(nzv)

# nz 리스트를 (N, 4) 배열로 변환 (각 열이 각 Sg에 해당)
nz = np.array(nz).T

# 그래프 그리기
plt.figure(figsize=(10, 6))
plt.plot(P, nz[:, 0], label='Sg=0.5')
plt.plot(P, nz[:, 1], ':', label='Sg=0.6')
plt.plot(P, nz[:, 2], '.-', label='Sg=0.7')
plt.plot(P, nz[:, 3], '--', label='Sg=0.8')

plt.xlabel('P(psia)')
plt.ylabel('Compressibility factor, Z')
plt.axis([100, 5000, 0, 1.1])
plt.legend(loc='best')
plt.grid(True)
plt.title('Compressibility Factor vs Pressure')
plt.show()