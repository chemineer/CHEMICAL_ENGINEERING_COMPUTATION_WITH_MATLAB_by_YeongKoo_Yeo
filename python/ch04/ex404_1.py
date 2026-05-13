from cubicEOSZ import cubicEOSZ

# 매틀랩 코드의 변수 할당
T = 350
P = 9.4573
Tc = 425.1
Pc = 37.96
w = 0.2
state = 'L'  # 액체 상태(Liquid) 지정
eos = 'rk'   # Redlich-Kwong 식 지정

# 함수 호출 (Z와 V를 각각 반환받음)
Z, V = cubicEOSZ(state, eos, T, P, Tc, Pc, w)

# 결과 출력
print(f"상태 방정식({eos}) - 상태({state}) 계산 결과:")
print(f"압축 인자 (Z) = {Z}")
print(f"몰 부피 (V) = {V}")