from hcRB import hcRB

# 1. 물성치 및 파라미터 설정
T = 373.15      # 온도 (K)
Tc = 535.55     # 임계 온도 (K)
w = 0.323       # 편심 인자 (Acentric factor)
Cpi = 120.496   # 이상 기체 열용량 (J/mol·K 등 단위에 따라 다름)

# 2. 열용량 계산 함수 호출 (Riedel-Brazhkin 등 상관식 사용)
# MATLAB: cpL = hcRB(T, Tc, w, Cpi)
cpL = hcRB(T, Tc, w, Cpi)

# 3. 결과 출력
print(f"Heat Capacity (cpL) of MEK at {T}K: {cpL}")