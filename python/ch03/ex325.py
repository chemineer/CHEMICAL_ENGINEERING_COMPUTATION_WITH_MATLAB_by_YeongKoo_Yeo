from vhwagner import vhwagner

# 매틀랩 코드의 변수 할당
C = [-7.670734, 1.965917, -2.445437, -2.899873]
T = 273.15
Tc = 508.1
Pc = 4.6924
vL = 7.145e-5
vV = 0.2453

# 함수 호출 및 결과 계산
dHv = vhwagner(T, Tc, Pc, vL, vV, C)

print(f"계산된 증발 엔탈피 (dHv): {dHv}")