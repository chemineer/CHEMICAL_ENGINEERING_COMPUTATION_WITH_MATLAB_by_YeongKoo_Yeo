from prVp import prVp

# 1. 온도 설정 (섭씨 85도)
# MATLAB: T = 85
T = 85

# 2. 증기압 계산 함수 호출
# MATLAB: pv = prVp(T, 'H2O')
pv = prVp(T, 'H2O')

# 3. 결과 출력
print(f"Vapor Pressure of Water at {T}°C: {pv}")