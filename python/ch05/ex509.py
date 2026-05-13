# 1. 상수 및 변수 정의
P1 = 205
P2 = 125
dP = P2 - P1
dz = -125
rho = 1000
g = 9.81
Ws = 1.2e6

# 2. 질량 유량(mdot) 계산
# MATLAB 식: mdot = -Ws / (dP/rho + g*dz)
# 분모가 0이 되는 상황을 방지하기 위해 나누기 전 확인하는 것이 좋습니다.
denominator = (dP / rho) + (g * dz)

if denominator != 0:
    mdot = -Ws / denominator
    print(f"Calculated mass flow rate (mdot): {mdot:.4f}")
else:
    print("Error: Division by zero. The denominator (dP/rho + g*dz) is zero.")