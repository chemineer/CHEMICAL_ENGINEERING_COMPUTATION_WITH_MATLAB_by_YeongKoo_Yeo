from visL import visL

# 1. 온도 변환 (섭씨 -> 켈빈)
# MATLAB: T = 150 + 273.15
T = 150 + 273.15

# 2. 점도 계산 함수 호출
# MATLAB: mu = visL(T, 'C6H6O')
mu = visL(T, 'C6H6O')

# 3. 결과 출력
print(f"Viscosity of Phenol at {T}K: {mu}")