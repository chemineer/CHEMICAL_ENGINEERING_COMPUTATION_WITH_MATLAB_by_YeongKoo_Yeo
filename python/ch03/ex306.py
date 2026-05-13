from costald import costald

# 1. n-헥산의 물성치 및 온도 설정
T = 293.15      # 온도 (K)
Tc = 507.8      # 임계 온도 (K)
vc = 386.8      # 임계 부피 (cm^3/mol)
w = 0.3002      # 편심 인자 (Acentric factor)
Mw = 86.178     # 분자량 (g/mol)

# 2. COSTALD 방법론을 이용한 액체 밀도 계산
# MATLAB: rhoL = costald(T,Tc,vc,w,Mw)
rhoL = costald(T, Tc, vc, w, Mw)

# 3. 결과 출력
print(f"Liquid Density (rhoL) of n-Hexane: {rhoL} g/cm^3")