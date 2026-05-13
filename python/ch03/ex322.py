# Vapor Pressure of Cyclohexanethiol

# 1. 온도 설정
# MATLAB: T = 100
T = 100

# 2. 4차 다항식을 이용한 증기압(Pv) 계산
# MATLAB: Pv = 1.35175e-6*T^4 - 2.8e-5*T^3 - 0.0053375*T^2 + 1.1674*T - 38.62
Pv = (1.35175e-6 * T**4) - (2.8e-5 * T**3) - (0.0053375 * T**2) + (1.1674 * T) - 38.62

# 3. 결과 출력
print(f"Vapor Pressure (Pv) of Cyclohexanethiol at T={T}: {Pv}")