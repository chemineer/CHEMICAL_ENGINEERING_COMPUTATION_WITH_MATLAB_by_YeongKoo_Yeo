from virialEOS import virialEOS

# 매틀랩 코드의 변수 할당
P = 14.8
T = 323.15
Pc = 48.08
Tc = 305.3
w = 0.1

# 함수 호출 (Z와 V를 각각 반환받음)
Z, V = virialEOS(P, T, Pc, Tc, w)

# 결과 출력
print(f"압축 인자 (Z) = {Z}")
print(f"몰 부피 (V) = {V}")