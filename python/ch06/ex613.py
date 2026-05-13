import numpy as np
from scipy.optimize import fsolve

# 1. 데이터 설정
V = 1
F = 25
Caf = 10
k1 = 50
k2 = 100
k3 = 10

# 2. 비선형 방정식 시스템 정의
def equations(Cs):
    # Cs[0] = Ca, Cs[1] = Cb
    Ca, Cb = Cs
    
    # f1: 성분 A에 대한 물질 수지 (입력 - 출력 + 생성/소멸 = 0)
    # -k1*Ca: 1차 반응 소멸, -k3*Ca^2: 2차 반응 소멸
    f1 = -k1 * Ca - k3 * (Ca**2) + F * (Caf - Ca) / V
    
    # f2: 성분 B에 대한 물질 수지
    # k1*Ca: A로부터 생성, -k2*Cb: B의 소멸
    f2 = k1 * Ca - k2 * Cb - F * Cb / V
    
    return [f1, f2]

# 3. 초기 추정값 설정 및 해 구하기
Cs0 = [5, 5]
Cs_solution = fsolve(equations, Cs0)

# 4. 결과 출력
Ca_ss, Cb_ss = Cs_solution
print(f"정상 상태 농도 Ca: {Ca_ss:.4f}")
print(f"정상 상태 농도 Cb: {Cb_ss:.4f}")