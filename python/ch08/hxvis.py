import math

def hxvis(muref, Tref, T, fstate):
    """
    주어진 두 개의 기준 온도와 점도를 사용하여 온도 T에서의 점도를 계산합니다.
    
    input:
      muref: 기준 온도 벡터 Tref에서의 점도 벡터 (Ns/m^2)
      Tref: 기준 온도 벡터 (K)
      T: 점도를 계산할 온도 (K)
      fstate: 유체 상태 (1: 액체(mu = A*exp(B/T)), 2: 기체(mu = A+BT))
    output:
      mu: 점도 (Ns/m^2)
    """
    if fstate == 1:  # 액체: 지수 모델 사용
        # MATLAB: mu = muref(1)*exp(Tref(2)*(Tref(1) - T)/(T*Tref(1) - Tref(2))*log(muref(2)/ muref(1)))
        term1 = Tref[1] * (Tref[0] - T)
        term2 = T * Tref[0] - Tref[1] # 주의: 원본 식의 분모 구조를 그대로 따름
        log_ratio = math.log(muref[1] / muref[0])
        
        mu = muref[0] * math.exp((term1 / term2) * log_ratio)
        
    else:  # 기체: 선형 모델 사용
        # MATLAB: mu = (muref(2)*Tref(1) - muref(1)*Tref(2) + T*(muref(1)-muref(2)))/(Tref(1) - Tref(2))
        mu = (muref[1] * Tref[0] - muref[0] * Tref[1] + T * (muref[0] - muref[1])) / (Tref[0] - Tref[1])
        
    return mu