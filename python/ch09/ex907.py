import control as ctrl

# 전달함수의 분자(n)와 분모(d) 계수 정의
# G(s) = (2s + 1) / (s^3 + 3s^2 + 2s + 1)
n = [2, 1]
d = [1, 3, 2, 1]

# 전달함수 생성
G = ctrl.tf(n, d)

# 전달함수를 상태공간 모델(A, B, C, D)로 변환
# sysss는 상태공간 객체입니다.
sysss = ctrl.tf2ss(G)

# 각 행렬 추출
A = sysss.A
B = sysss.B
C = sysss.C
D = sysss.D

# 결과 출력
print("A 행렬:\n", A)
print("B 행렬:\n", B)
print("C 행렬:\n", C)
print("D 행렬:\n", D)