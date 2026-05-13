import numpy as np
from scipy.optimize import fsolve

def rxnfun(x, T):
    """
    화학 평형 방정식을 정의하는 함수
    x[0:9]: 각 성분의 몰수 (CH4, C2H4, C2H2, CO2, CO, O2, H2, H2O, C2H6)
    x[9:12]: 라그랑주 승수 (lambda1, lambda2, lambda3)
    """
    R = 1.9872
    # x(1:9) -> x[0:9]. sum은 9번째 인덱스 직전까지 계산
    ns = np.sum(x[0:9])
    
    # 결과 배열 초기화 (12개의 방정식을 담을 리스트)
    f = np.zeros(12)
    
    # 1~3번: 원자 수지 보존 법칙 (Atom Balances)
    f[0] = 2*x[3] + x[4] + 2*x[5] + x[7] - 4
    f[1] = 4*x[0] + 4*x[1] + 2*x[2] + 2*x[6] + 2*x[7] + 6*x[8] - 14
    f[2] = x[0] + 2*x[1] + 2*x[2] + x[3] + x[4] + 2*x[8] - 2
    
    # 4~12번: 평형 조건 (Chemical Equilibrium via Gibbs Energy Minimization)
    # x[9], x[10], x[11]은 각각 매틀랩의 x(10), x(11), x(12)에 해당
    f[3] = x[0] - ns * np.exp(-(4.61e3/(R*T) + 1 - x[0]/ns + 4*x[10] + x[11]))
    f[4] = x[1] - ns * np.exp(-(28.249e3/(R*T) + 1 - x[1]/ns + 4*x[10] + 2*x[11]))
    f[5] = x[2] - ns * np.exp(-(40.604e3/(R*T) + 1 - x[2]/ns + 2*x[10] + 2*x[11]))
    f[6] = x[3] - ns * np.exp(-(-94.61e3/(R*T) + 1 - x[3]/ns + 2*x[9] + x[11]))
    f[7] = x[4] - ns * np.exp(-(-47.942e3/(R*T) + 1 - x[4]/ns + x[9] + x[11]))
    f[8] = x[5] - ns * np.exp(-(1 - x[5]/ns + 2*x[9]))
    f[9] = x[6] - ns * np.exp(-(1 - x[6]/ns + 2*x[10]))
    f[10] = x[7] - ns * np.exp(-(-46.03e3/(R*T) + 1 - x[7]/ns + x[9] + 2*x[10]))
    f[11] = x[8] - ns * np.exp(-(26.13e3/(R*T) + 1 - x[8]/ns + 6*x[10] + 2*x[11]))
    
    return f

# --- 실행 부분 (ex216.m) ---

# 초기 추측값 설정
x0 = np.array([0.001, 0.001, 0.001, 0.993, 1, 0.0001, 5.992, 1, 0.001, 10, 10, 10])
T_val = 1000

# fsolve 실행 (추가 인자 T_val 전달)
x_sol = fsolve(rxnfun, x0, args=(T_val,))

# 결과 출력
comp_names = ['CH4', 'C2H4', 'C2H2', 'CO2', 'CO', 'O2', 'H2', 'H2O', 'C2H6']
lambda_names = ['lambda1', 'lambda2', 'lambda3']

print(f"\n{'i':<3} {'Comp.':<10} {'Initial Val.':<15} {'Final Val.':<15}")
for i in range(len(comp_names)):
    print(f"{i+1:<3} {comp_names[i]:<10} {x0[i]:<15.9f} {x_sol[i]:<15.9f}")

print(f"\n{'i':<3} {'Lambda':<10} {'Initial Val.':<15} {'Final Val.':<15}")
for i in range(len(lambda_names)):
    idx = i + len(comp_names)
    print(f"{idx+1:<3} {lambda_names[i]:<10} {x0[idx]:<15.1f} {x_sol[idx]:<15.9f}")