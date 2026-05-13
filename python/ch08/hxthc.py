def hxthc(xkref, Tref, T):
    """
    주어진 두 개의 기준 온도와 열전도도를 사용하여 온도 T에서의 열전도도를 계산합니다.
    (선형 보간: xk = A + B*T)
    
    input:
      xkref: 기준 온도 벡터 Tref에서의 열전도도 벡터 (W/m/K)
      Tref: 기준 온도 벡터 (K)
      T: 열전도도를 계산할 온도 (K)
    output:
      xk: 열전도도 (W/m/K)
    """
    # MATLAB 식: (xkref(2)*Tref(1) - xkref(1)*Tref(2) + T*(xkref(1)-xkref(2))) / (Tref(1) - Tref(2))
    # 파이썬 인덱스 적용: xkref(1) -> xkref[0], xkref(2) -> xkref[1]
    
    xk = (xkref[1] * Tref[0] - xkref[0] * Tref[1] + T * (xkref[0] - xkref[1])) / (Tref[0] - Tref[1])
    
    return xk