from deptfun import deptfun

# 매틀랩 코드의 상수 및 초기값 할당
eosset = ['VR', 'VDW', 'RK', 'SRK', 'PR']
Tc = 425.1
Pc = 37.96
w = 0.2
T = 500
P = 50
state = 'v'

# 각 상태 방정식에 대해 순차적으로 계산 및 출력
for eos in eosset:
    # deptfun 함수 호출 (Z, V, dH, dS를 각각 반환받음)
    Z, V, dH, dS = deptfun(state, eos, T, P, Tc, Pc, w)
    
    # 결과 출력 (매틀랩의 fprintf와 유사한 f-string 사용)
    print(f"The equation of state={eos}: Z={Z} H^R={dH} S^R={dS}")