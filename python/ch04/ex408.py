from deptfun import deptfun
from delHS import delHS

# 매틀랩 코드의 변수 할당
Tc = 369.8
Pc = 42.49
w = 0.152
T1 = 378.15
P1 = 5
T2 = 463.15
P2 = 25
eos = 'pr'
state = 'v'

# 상태 1에서의 잔류 특성 계산
Z1, V1, dH1, dS1 = deptfun(state, eos, T1, P1, Tc, Pc, w)
print(f"상태 1 계산 결과: Z={Z1}, V={V1}, dH={dH1}, dS={dS1}")

# 상태 2에서의 잔류 특성 계산
Z2, V2, dH2, dS2 = deptfun(state, eos, T2, P2, Tc, Pc, w)
print(f"상태 2 계산 결과: Z={Z2}, V={V2}, dH={dH2}, dS={dS2}")

# 열용량 계수 및 파라미터 재설정
A = -4.224
B = 0.3063
C = -1.586e-4
D = 3.215e-8

# 엔탈피 및 엔트로피 변화량 계산 (delHS 함수 호출)
dH, dS = delHS(state, eos, T1, P1, T2, P2, A, B, C, D, Tc, Pc, w)

# 최종 결과 출력
print("-" * 30)
print(f"계산된 엔탈피 변화량 (dH) = {dH}")
print(f"계산된 엔트로피 변화량 (dS) = {dS}")