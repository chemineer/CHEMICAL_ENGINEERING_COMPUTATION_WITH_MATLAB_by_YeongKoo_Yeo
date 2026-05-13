from vpwagner import vpwagner
import numpy as np

# 1. 물성치 및 파라미터 설정
# MATLAB: z = [-7.670734 1.965917 -2.445437 -2.899873];
z = np.array([-7.670734, 1.965917, -2.445437, -2.899873])

# MATLAB: T = 273.15; Tc = 508.1; Pc = 4.6924;
T = 273.15       # 온도 (K)
Tc = 508.1       # 임계 온도 (K)
Pc = 4.6924      # 임계 압력 (bar)

# 2. Wagner 식을 이용한 증기압 계산 함수 호출
# MATLAB: Pv = vpwagner(T, Tc, Pc, z)
Pv = vpwagner(T, Tc, Pc, z)

# 3. 결과 출력
print(f"Vapor Pressure (Pv) of Acetone at {T}K: {Pv}")