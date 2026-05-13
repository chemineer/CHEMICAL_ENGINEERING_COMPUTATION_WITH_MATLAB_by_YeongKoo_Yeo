import control as ctrl

# 1. 전방 경로(Forward path) G(s)의 분자와 분모 정의
# G(s) = (2s + 1) / (s^2 + 3s + 2)
ng = [2, 1]
dg = [1, 3, 2]
G = ctrl.tf(ng, dg)

# 2. 피드백 경로(Feedback path) H(s)의 분자와 분모 정의
# H(s) = 1 / (s + 1)
nh = [1]
dh = [1, 1]
H = ctrl.tf(nh, dh)

# 3. 피드백 결합 (Negative Feedback: sign=-1)
# MATLAB의 feedback(ng, dg, nh, dh, -1)과 동일
# 파이썬 control.feedback은 기본값이 음의 피드백(sign=-1)입니다.
Gcl = ctrl.feedback(G, H)

# 결과 출력
print("폐루프 전달함수 Gcl(s):")
print(Gcl)