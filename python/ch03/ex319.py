import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 1. 데이터 설정
# T: 섭씨 온도 배열에 273.15를 더해 켈빈(K)으로 변환
T = np.array([-15, -4.3, 7.5, 20.7, 29.1, 40.7, 58.1, 78.0, 99.2]) + 273.15
# Pv: 증기압 (mmHg)
Pv = np.array([0.667, 1.333, 2.666, 5.333, 8.000, 13.33, 26.66, 53.33, 101.32])

# 2. Antoine 식 매개변수 추정
# Xa = [1, 1/T, log(Pv)/T]
Xa = np.column_stack([np.ones(len(T)), 1/T, np.log(Pv)/T])
ba = np.linalg.lstsq(Xa, np.log(Pv), rcond=None)[0]
Aa, Ca, Ba = ba[0], -ba[2], ba[0]*(-ba[2]) - ba[1]

# 3. Riedel 식 매개변수 추정
# Xr = [1, 1/T, log(T), T^6]
Xr = np.column_stack([np.ones(len(T)), 1/T, np.log(T), T**6])
br = np.linalg.lstsq(Xr, np.log(Pv), rcond=None)[0]
Ar, Br, Cr, Dr = br

# 4. Harlecher-Braun 식 매개변수 추정
# Xh = [1, 1/T, log(T), Pv/T^2]
Xh = np.column_stack([np.ones(len(T)), 1/T, np.log(T), Pv/T**2])
bh = np.linalg.lstsq(Xh, np.log(Pv), rcond=None)[0]
Ah, Bh, Ch, Dh = bh

# 5. 그래프 생성을 위한 데이터 생성
Ti = np.arange(T[0], T[-1] + 1)
Pa = np.exp(Aa - Ba / (Ti + Ca))
Pr = np.exp(Ar + Br/Ti + Cr*np.log(Ti) + Dr*Ti**6)

# Harlecher-Braun 식: 비선형 방정식 풀이
# MATLAB의 for k = 1:length(Ti) 루프를 그대로 구현
Ph = []
n_ti = len(Ti)
for k in range(n_ti):
    # MATLAB: P0 = 100/length(Ti) * k (k는 1부터 시작)
    p_guess = (100 / n_ti) * (k + 1)
    
    # 식: Ah + Bh/Ti + Ch*log(Ti) + Dh*x/Ti^2 - log(x) = 0
    t_val = Ti[k]
    func = lambda x: Ah + Bh/t_val + Ch*np.log(t_val) + Dh*x/t_val**2 - np.log(x)
    
    # 수치 해법
    root = fsolve(func, p_guess)
    Ph.append(root[0])

# 6. 시각화
plt.figure(figsize=(10, 6))
plt.plot(Ti, Pa, label='Antoine')
plt.plot(Ti, Pr, ':', label='Riedel')
plt.plot(Ti, Ph, '.-', label='Harlecher-Braun')
plt.plot(T, Pv, 'o', label='Data')

plt.xlabel('T(K)') 
plt.ylabel('Pv(mmHg)')
plt.legend(loc='best')
plt.grid(True)
plt.show()