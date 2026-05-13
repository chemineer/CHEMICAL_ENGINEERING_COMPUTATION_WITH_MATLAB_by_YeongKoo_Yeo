from hcapL import hcapL

# 1. 온도 변환 (섭씨 -> 켈빈)
# MATLAB: T = 150 + 273.15
T = 150 + 273.15

# 2. 액체 열용량 계산 함수 호출
# MATLAB: hcp = hcapL(T, 'C6H6O')
hcp = hcapL(T, 'C6H6O')

# 3. 결과 출력
print(f"Heat Capacity of Phenol at {T}K: {hcp}")