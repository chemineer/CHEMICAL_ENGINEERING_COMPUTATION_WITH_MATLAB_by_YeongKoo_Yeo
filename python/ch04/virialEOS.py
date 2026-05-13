import numpy as np

def virialEOS(P, T, Pc, Tc, w):
    """
    Virial 상태 방정식을 사용하여 주어진 T, P에서 
    압축 인자(Z)와 몰 부피(V)를 추정합니다.
    
    입력:
    P, Pc: 압력 및 임계 압력 (atm)
    T, Tc: 온도 및 임계 온도 (K)
    w: 편심 인자 (acentric factor)
    
    출력:
    Z: 압축 인자
    V: 몰 부피 (L/gmol)
    """
    
    # 기체 상수 (atm-liter/(gmol-K))
    R = 0.08206
    
    # 환산 온도 및 환산 압력 계산
    # (입력이 리스트나 배열일 경우를 대비해 numpy 연산 사용)
    Tr = np.array(T) / np.array(Tc)
    Pr = np.array(P) / np.array(Pc)
    
    # B0 및 B1 계수 계산
    B0 = 0.083 - 0.422 / (Tr**1.6)
    B1 = 0.139 - 0.172 / (Tr**4.2)
    
    # 제2 비리알 계수 B 및 압축 인자 Z 계산
    # MATLAB 코드의 로직을 그대로 따름: B = (R*Tc/Pc) * (B0 + w*B1)
    B = R * Tc * (B0 + w * B1) / Pc
    Z = 1 + B * P / (R * T)
    
    # 몰 부피 V 계산
    V = R * Z * T / P
    
    return Z, V

# 사용 예시:
# z_val, v_val = virialEOS(1.0, 300.0, 45.4, 190.6, 0.011)
# print(f"Z: {z_val}, V: {v_val}")

