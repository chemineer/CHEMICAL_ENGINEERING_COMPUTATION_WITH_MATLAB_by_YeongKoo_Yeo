import numpy as np

def hxcp(cpref, Tref, T):
    """
    온도 T에서의 열용량(Cp) 계산 (Python version)
    Cp = A + B*T (선형 보간)
    """
    # cpref, Tref는 길이가 2인 리스트 또는 numpy 배열
    # Cp = (cpref(2)*Tref(1) - cpref(1)*Tref(2) + T*(cpref(1)-cpref(2)))/(Tref(1) - Tref(2))
    
    Cp = (cpref[1] * Tref[0] - cpref[0] * Tref[1] + T * (cpref[0] - cpref[1])) / (Tref[0] - Tref[1])
    return Cp