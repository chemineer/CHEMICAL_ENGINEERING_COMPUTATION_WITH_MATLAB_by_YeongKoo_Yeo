import numpy as np

# 1. 주어진 물리 상수 및 변수 정의
P = 12        # 압력 (Pressure, atm)
T = 315.6     # 온도 (Temperature, K)
R = 0.08205   # 기체 상수 (Gas constant)
a = 3.592     # 반데르발스 상수 a
b = 0.04267   # 반데르발스 상수 b

# 2. 부피 V에 대한 3차 방정식 계수 설정
# 방정식 형태: P*V^3 - (b*P + R*T)*V^2 + a*V - a*b = 0
coeffs = [P, -(b * P + R * T), a, -a * b]

# 3. 다항식의 근(부피 V) 계산
volumes = np.roots(coeffs)

# 4. 결과 출력
print(f"Specific volume roots for CO2 at P={P} atm, T={T} K:")
for i, v in enumerate(volumes):
    # 실제 물리적인 의미를 갖는 근은 양의 실수입니다.
    if np.isreal(v):
        print(f"Root {i+1} (Real): {v.real:.5f}")
    else:
        print(f"Root {i+1} (Complex): {v:.5f}")