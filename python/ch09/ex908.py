import control as ctrl

# 전달함수 계수 정의
# G(s) = (s + 2) / (s^2 + 1.5s + 1)
n = [1, 2]
d = [1, 1.5, 1]

# 전달함수 생성 및 상태공간(State-Space) 모델로 변환
# tf2ss는 시스템을 상태공간 행렬 A, B, C, D로 반환합니다.
sys_ss = ctrl.tf2ss(n, d)

# 결과 확인
print("A 행렬:\n", sys_ss.A)
print("B 행렬:\n", sys_ss.B)
print("C 행렬:\n", sys_ss.C)
print("D 행렬:\n", sys_ss.D)