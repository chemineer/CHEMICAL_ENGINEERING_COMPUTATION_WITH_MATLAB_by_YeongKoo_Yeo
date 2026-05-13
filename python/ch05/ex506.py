import numpy as np
import matplotlib.pyplot as plt

# 1. 데이터 정의
Nre = np.array([8500, 20000, 30000, 60000, 700000, 1000000, 10000000])
f = np.array([0.008, 0.0065, 0.006, 0.005, 0.003, 0.0028, 0.002])

# 2. 선형 회귀 (log-log 스케일)
x = np.log10(Nre)
y = np.log10(f)

# polyfit(x, y, 1)은 1차 다항식(y = ax + b)의 계수를 찾습니다.
# p[0]은 기울기(b), p[1]은 절편(log10(a))에 해당합니다.
p = np.polyfit(x, y, 1)

# 식 f = a * Nre^b 를 위해 a와 b 추출
a = 10**p[1]
b = p[0]

print(f"a = {a:.6f}, b = {b:.6f}")

# 3. 상관관계 계산
Nrex = np.linspace(min(Nre), max(Nre), 100)
fvx = a * (Nrex**b)

# 4. 그래프 그리기
plt.figure(figsize=(8, 5))
plt.plot(Nre, f, 'o', label='Data')
plt.plot(Nrex, fvx, label='Correlation')

# 로그 스케일 적용 (로그-로그 그래프)
plt.xscale('log')
plt.yscale('log')

plt.xlabel('N_{Re} (Reynolds number)')
plt.ylabel('f (Friction factor)')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.show()