def hxrho(rhoref, Tref, T, fstate):
    """
    주어진 두 개의 기준 온도와 밀도를 사용하여 온도 T에서의 밀도를 계산합니다.
    
    input:
      rhoref: 기준 온도 벡터 Tref에서의 밀도 벡터 (kg/m^3)
      Tref: 기준 온도 벡터 (K)
      T: 밀도를 계산할 온도 (K)
      fstate: 유체 상태 (1: 액체(선형 보간), 2: 기체(이상기체 상태방정식 유사))
    output:
      rho: 밀도 (kg/m^3)
    """
    if fstate == 1:  # 액체: 선형 보간 사용
        # MATLAB 식: (rhoref(2)*Tref(1) - rhoref(1)*Tref(2) + T*(rhoref(1)-rhoref(2))) / (Tref(1) - Tref(2))
        # 파이썬 인덱스: rhoref[1]은 MATLAB의 rhoref(2)에 해당
        rho = (rhoref[1] * Tref[0] - rhoref[0] * Tref[1] + T * (rhoref[0] - rhoref[1])) / (Tref[0] - Tref[1])
    else:  # 기체: rho = A/T 모델
        # MATLAB 식: rhoref(1)*Tref(1)/T
        rho = rhoref[0] * Tref[0] / T
        
    return rho