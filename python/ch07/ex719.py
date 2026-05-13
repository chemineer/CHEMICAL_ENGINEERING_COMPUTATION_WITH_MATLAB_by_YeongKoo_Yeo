import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

# 외부 모듈 임포트
from ssdist import ssdist
from dyndist import dyndist

# 1. 파라미터 설정 (dpar)
dpar = {
    'alpha': 1.5, 'n': 30, 'nf': 15, 'F': 1, 'zf': 0.5, 'q': 1,
    'R': 2.7, 'Vs': 3.2, 'D': 0.5, 'md': 5, 'mb': 5, 'mt': 0.5
}

# 2. 외란 설정 (dels)
dels = {
    'delR': 0.01 * dpar['R'], 'delRt': 10, 
    'delV': 0, 'delVt': 0, 
    'delz': 0, 'delzt': 0,
    'delF': 0, 'delFt': 0
}

# 3. 초기 정상 상태 계산 (fsolve)
x0_guess = 0.5 * np.ones(dpar['n'])
# fsolve를 사용하여 ssdist(x) = 0인 지점을 찾습니다.
x0 = fsolve(ssdist, x0_guess, args=(dpar,))

# 음수 값 방지 (MATLAB의 for 루프 로직)
x0 = np.abs(x0)

# 4. 동적 시뮬레이션 (solve_ivp)
t0, tf = 0, 400
t_eval = np.linspace(t0, tf, 500)

# solve_ivp는 (t, y) 순서의 함수를 기대하므로 lambda를 사용하거나 함수 구조를 맞춰야 함
sol = solve_ivp(dyndist, [t0, tf], x0, args=(dpar, dels), t_eval=t_eval, method='RK45')

# 5. 결과 시각화
plt.figure(figsize=(12, 5))

# Subplot 1: 시간에 따른 x_D(첫 번째 단) 및 x_B(마지막 단) 변화
plt.subplot(1, 2, 1)
plt.plot(sol.t, sol.y[0, :], label='x_D')
plt.plot(sol.t, sol.y[-1, :], '--', label='x_B')
plt.xlabel('t(min)')
plt.ylabel('x')
plt.legend()
plt.grid(True)

# Subplot 2: 최종 시간에서의 단별 농도 분포
plt.subplot(1, 2, 2)
nv = np.arange(1, dpar['n'] + 1)
plt.plot(nv, sol.y[:, -1])
plt.xlabel('n')
plt.ylabel('x_i')
plt.axis('tight')
plt.grid(True)

plt.tight_layout()
plt.show()