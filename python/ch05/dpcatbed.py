import math

def dpcatbed(W, Bd, Bl, Pd, Pl, ep, mu, rho):
    """
    충전층을 통과하는 흐름에 대한 압력 강하를 계산합니다.
    
    입력 매개변수:
    W: 질량 유량 (lb/h)
    Bd: Bed 직경 (ft)
    Bl: Bed 길이 (ft)
    Pd: 입자 직경 (in)
    Pl: 입자 길이 (in)
    ep: 공극률 (void fraction)
    mu: 점도 (lb/ft-h 또는 관련 단위)
    rho: 밀도 (lb/ft^3)
    
    출력:
    dPt: 압력 강하
    """
    
    # 중력 가속도 환산 계수 (lb_m-ft/lb_f-h^2)
    gc = 4.17e8 
    
    # 단면적 및 겉보기 질량 유속 (superficial mass flow rate)
    A = math.pi * Bd**2 / 4 
    G = W / A # lb/h/ft^2
    
    # 입자의 표면적(Ap, ft^2) 및 부피(Vp, ft^3) 계산
    # 144와 1728은 sq.in -> sq.ft, cu.in -> cu.ft 환산 계수
    Ap = (math.pi * Pd**2 / 2 + math.pi * Pd * Pl) / 144 
    Vp = math.pi * Pl * Pd**2 / (4 * 1728) 
    
    # 비표면적(S) 및 유효 입자 직경(ePd)
    S = Ap * (1 - ep) / Vp 
    ePd = 6 * (1 - ep) / S 
    
    # 레이놀즈 수 (Nre) 계산
    # 2.419는 점도 단위 환산과 관련된 상수로 보임
    Nre = G * ePd / (2.419 * mu * (1 - ep)) 
    
    # 마찰 계수 (fP) 결정 (Ergun 방정식의 구성 요소)
    if Nre < 1:
        fP = 150 / Nre 
    elif Nre < 1e4:
        fP = 150 / Nre + 1.75 
    else:
        fP = 1.75 
        
    # 최종 압력 강하 (dPt) 계산
    # MATLAB 원본의 복잡한 수식을 파이썬 문법으로 변환
    term1 = 150 * 2.419 * mu * (1 - ep) / (ePd * G)
    term2 = 1.75
    
    numerator = Bl * (1 - ep) * (G**2) * (term1 + term2)
    denominator = 144 * (ep**3) * ePd * gc * rho
    
    dPt = numerator / denominator 
    
    return dPt

# 사용 예시:
# dp_val = dpcatbed(1000, 2, 10, 0.5, 0.5, 0.4, 0.05, 62.4)
# print(f"Pressure Drop: {dp_val}")