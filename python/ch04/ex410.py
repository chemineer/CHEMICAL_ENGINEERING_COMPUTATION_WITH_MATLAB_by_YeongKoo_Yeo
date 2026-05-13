from phigas import phigas

# 매틀랩 코드의 변수 할당
Tc = 308.3
Pc = 61.39
w = 0.187
T = 250
P = 10
eos = 'srk'
state = 'v'

# 함수 호출 (phig와 f를 각각 반환받음)
# 매틀랩의 [phig f] = phigas(...)와 동일하게 처리
phig, f = phigas(state, eos, T, P, Tc, Pc, w)

# 결과 출력
print(f"상태 방정식({eos}) - 상태({state}) 계산 결과:")
print(f"퓨가시티 계수 (phig) = {phig}")
print(f"퓨가시티 (f) = {f}")