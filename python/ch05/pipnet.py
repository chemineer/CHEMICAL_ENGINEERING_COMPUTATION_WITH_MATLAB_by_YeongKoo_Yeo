import numpy as np

def pipnet(m, rho, mu, D, L, rf):
    # 수치적 안정성을 위해 작은 값(eps)을 더함
    eps = 1e-10
    m_safe = np.abs(m) + eps
    
    # 기초 유체 역학 계산
    v = 4 * m_safe / (np.pi * D**2 * rho)
    Nre = D * v * rho / mu
    
    # 마찰 계수 f 계산 (Chen 방정식)
    C = (7 / Nre)**0.9 + 0.27 * rf / D
    C = np.maximum(C, 1e-10) # 로그 함수 오류 방지
    
    A = (2.457 * np.log(C))**16
    B = (37530 / Nre)**16
    f = 2 * ((8 / Nre)**12 + 1 / ((A + B)**1.5))**(1/12)
    
    # 도함수 성분 계산
    dC = -(0.9 / Nre) * (7 / Nre)**0.9
    dA = (39.312 / C) * dC * (2.457 * np.log(C))**15
    dB = -16 * B / Nre
    
    # 파라미터 b, c, a, phi 계산
    b = (Nre / 8) / ((f / 2)**12) * ((8 / Nre)**13 + (dA + dB) / ((A + B)**2.5))
    c = 32 * L / (rho * np.pi**2 * D**5)
    
    # m_safe를 사용하여 파워 연산 안정화
    a = f * (m_safe**(b - 1))
    phi = a * c
    
    # 잔차(fm) 벡터 초기화
    fm = np.zeros(12)
    
    # 질량 보존 방정식
    fm[0] = m[0] + m[3] - 334
    fm[1] = m[0] - m[1] - m[4]
    fm[2] = m[2] - m[3] + m[7]
    fm[3] = m[1] + m[2] - m[6] - m[8]
    fm[4] = m[4] - m[5] - 42
    fm[5] = m[5] + m[6] - m[11] - 108
    fm[6] = m[10] + m[11] - 9
    fm[7] = m[8] + m[9] - m[10] - 88
    fm[8] = m[7] - m[9] - 87
    
    # 에너지 손실 방정식 (m_safe 사용)
    fm[9] = phi[0]*m_safe[0]**(2-b[0]) + phi[1]*m_safe[1]**(2-b[1]) - phi[2]*m_safe[2]**(2-b[2]) - phi[3]*m_safe[3]**(2-b[3])
    fm[10] = -phi[1]*m_safe[1]**(2-b[1]) + phi[4]*m_safe[4]**(2-b[4]) + phi[5]*m_safe[5]**(2-b[5]) - phi[6]*m_safe[6]**(2-b[6])
    fm[11] = phi[2]*m_safe[2]**(2-b[2]) - phi[7]*m_safe[7]**(2-b[7]) + phi[8]*m_safe[8]**(2-b[8]) - phi[9]*m_safe[9]**(2-b[9])
    
    return fm