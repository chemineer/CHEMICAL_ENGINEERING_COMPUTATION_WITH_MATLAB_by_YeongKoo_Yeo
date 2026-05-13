import numpy as np
import matplotlib.pyplot as plt
from lagrintp import lagrintp
from newtintp import newtintp
from cspintp import cspintp


# --- 2. 메인 실행 데이터 ---

t = [0, 2, 5, 6, 13, 20, 24, 30, 35, 41, 50]
c = [0.86, 0.61, 0.47, 0.39, 0.25, 0.18, 0.15, 0.12, 0.10, 0.09, 0.08]
ti = [8, 15, 25, 32]

# 각 방법으로 보간 수행
cL = lagrintp(t, c, ti)
cN = newtintp(t, c, ti)
cC = cspintp(t, c, ti)

# 결과 출력
print(f"Lagrange method: C = {' '.join(map(lambda x: f'{x:g}', cL))}")
print(f"Newton method:   C = {' '.join(map(lambda x: f'{x:g}', cN))}")
print(f"Cubic spline:    C = {' '.join(map(lambda x: f'{x:g}', cC))}")

# --- 3. 시각화 ---

plt.figure(figsize=(10, 6))
plt.plot(t, c, 'o-', label='Data')
plt.plot(ti, cL, '*', markersize=10, label='Lagrange')
plt.plot(ti, cN, '<', markersize=8, label='Newton')
plt.plot(ti, cC, 'kd', markersize=6, label='Cubic spline')

plt.xlabel('t(min)')
plt.ylabel('C_A(mol/l)')
plt.title('Comparison of Interpolation Methods')
plt.legend(loc='best')
plt.grid(True)
plt.show()