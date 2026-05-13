import scipy.signal as signal

# MATLAB의 N(분자), D(분모) 다항식 계수 정의
# N(s) = 1*s^3 + 5*s^2 + 9*s + 7
# D(s) = 1*s^2 + 3*s + 2
N = [1, 5, 9, 7]
D = [1, 3, 2]

# 부분분수 전개 수행
# R: 잔차(residues), P: 극점(poles), K: 직접항(direct terms/polynomial part)
R, P, K = signal.residue(N, D)

print("잔차 (Residues):", R)
print("극점 (Poles):", P)
print("직접항 (Direct terms):", K)