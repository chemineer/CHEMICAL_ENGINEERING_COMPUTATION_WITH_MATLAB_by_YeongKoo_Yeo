from corstsg import corstsg

# 1. 물성치 및 온도 설정
# MATLAB: T = 303.15; Tc = 499; Tb = 308.15; Pc = 54.9;
T = 303.15      # 온도 (K)
Tc = 499        # 임계 온도 (K)
Tb = 308.15     # 끓는점 (K)
Pc = 54.9       # 임계 압력 (bar)

# 2. 표면 장력 계산 함수 호출 (Correlations for Surface Tension)
# MATLAB: sg = corstsg(Pc, Tc, Tb, T)
sg = corstsg(Pc, Tc, Tb, T)

# 3. 결과 출력
print(f"Surface Tension (sg) of Ethanethiol at {T}K: {sg}")