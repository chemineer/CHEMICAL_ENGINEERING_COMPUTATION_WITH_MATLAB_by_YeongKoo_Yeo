from hcapG import hcapG

# 1. 온도 설정 (섭씨 300도를 켈빈 온도로 변환)
# MATLAB: T = 300 + 273.15
T = 300 + 273.15

# 2. 기체 열용량 계산 함수 호출
# MATLAB: v = hcapG(T, 'CO2')
v = hcapG(T, 'CO2')

# 3. 결과 출력
print(f"Heat Capacity of CO2 at {T}K: {v}")