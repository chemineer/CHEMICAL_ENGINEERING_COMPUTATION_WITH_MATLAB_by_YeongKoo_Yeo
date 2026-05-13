from zLK import zLK

# 매틀랩 코드의 변수 할당
T = 400
Tc = 419.6
P = 20
Pc = 40.2
w = 0.191

# 함수 호출을 통한 Z 값 계산
Z = zLK(T, Tc, P, Pc, w)

# 결과 출력
print(f"1-Butene의 압축 인자 (Z) = {Z}")