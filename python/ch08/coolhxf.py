import numpy as np
from coolhxdat import *

def coolhxf(b, data):
    # 데이터 압축 해제 (coolhxdat의 결과물 사용)
    Q, Ui, F12, dTlm = data['Q'], data['Ui'], data['F12'], data['dTlm']
    At, Nt, Di, Do, L = data['At'], data['Nt'], data['Di'], data['Do'], data['L']
    Pt, Ds, cl, mc, kw, dx, Dlm = data['Pt'], data['Ds'], data['cl'], data['mc'], data['kw'], data['dx'], data['Dlm']
    muc, muwc, kc, Cpc = data['muc'], data['muwc'], data['kc'], data['Cpc']
    rhoh, ui, muh, muwh, kh, Cph = data['rhoh'], data['ui'], data['muh'], data['muwh'], data['kh'], data['Cph']
    Rfi, Rfo = data['Rfi'], data['Rfo']

    # 계산 로직
    Ai = Q / (Ui * F12 * dTlm)
    Np = np.ceil(Ai / (At * Nt))
    Ai = Np * Nt * np.pi * Di * L
    
    De = (4 / (np.pi * Do)) * (Pt**2 - np.pi * Do**2 / 4)
    Acf = Ds * cl * b / Pt
    Go = mc / Acf
    
    # 쉘 측(shell-side) 계수 계산
    Nreo = De * Go / (muc * 3600 / 1488)
    Npro = Cpc * (muc * 3600 / 1488) / kc
    Nuo = 0.36 * Nreo**0.55 * Npro**(1/3) * (muc / muwc)**0.14
    ho = Nuo * kc / De
    
    # 튜브 측(tube-side) 계수 계산
    Nrei = Di * rhoh * ui / (muh / 1488)
    Npri = Cph * (muh * 3600 / 1488) / kh
    Nui = 0.027 * Nrei**0.8 * Npri**(1/3) * (muh / muwh)**0.14
    hi = Nui * kh / Di
    
    # 총괄 열전달 계수 재계산 및 잔차 반환
    Ui_new = 1 / ((Di/Do)/ho + (Di*dx)/(Dlm*kw) + 1/hi + Rfi + (Di/Do)*Rfo)
    Qd = Ai * Ui_new * F12 * dTlm
    
    return Q - Qd