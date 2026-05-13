from visG import visG

# 1. 온도 변환 (섭씨 -> 켈빈)
# MATLAB: T = 100 + 273.15
T = 100 + 273.15

# 2. 기체 점도 계산 함수 호출
# MATLAB: mu = visG(T, 'methane')
mu = visG(T, 'methane')

# 3. 결과 출력
print(f"Viscosity of Methane at {T}K: {mu}")