import numpy as np
from heatH2Odat import heatH2Odat # 데이터 로드 함수 임포트

def heatH2O(b):
    """
    heatH2O.m: H2O 가열기 설계를 위한 잔차(Residual) 계산 함수
    """
    # 데이터 가져오기
    data = heatH2Odat()
    
    # 변수 언패킹 (수식 가독성을 위해)
    Qr, Ui_init, F12, dTlm = data['Qr'], data['Ui_init'], data['F12'], data['dTlm']
    At, Nt, Di, Do, L = data['At'], data['Nt'], data['Di'], data['Do'], data['L']
    Ds, cl, Pt, mh = data['Ds'], data['cl'], data['Pt'], data['mh']
    De, muh, muwh, kh, Cph = data['De'], data['muh'], data['muwh'], data['kh'], data['Cph']
    rhoc, ui, muc, muwc, kc, Cpc = data['rhoc'], data['ui'], data['muc'], data['muwc'], data['kc'], data['Cpc']
    dx, Dlm, kw, Rfi, Rfo = data['dx'], data['Dlm'], data['kw'], data['Rfi'], data['Rfo']

    # 1. 초기 설계 및 패스 계산
    Ai = Qr / (Ui_init * F12 * dTlm)
    Np = np.ceil(Ai / (At * Nt))
    
    # 2. Ai 및 Ui 재조정
    Ai = Np * Nt * np.pi * Di * L
    Ui_adj = Qr / (Ai * F12 * dTlm)
    
    # 3. 쉘 측(Shell-side) 계산
    Acf = Ds * cl * b / Pt
    Go = mh / Acf
    Nreo = De * Go / (muh * 3600 / 1488)
    Npro = Cph * (muh * 3600 / 1488) / kh
    Nuo = 0.36 * Nreo**0.55 * Npro**(1/3) * (muh / muwh)**0.14
    ho = Nuo * kh / De
    
    # 4. 튜브 측(Tube-side) 계산
    Nrei = Di * rhoc * ui / (muc / 1488)
    Npri = Cpc * (muc * 3600 / 1488) / kc
    Nui = 0.027 * Nrei**0.8 * Npri**(1/3) * (muc / muwc)**0.14
    hi = Nui * kc / Di
    
    # 5. 새로운 총괄 열전달 계수 및 잔차 계산
    Ui_new = 1 / ((Di / Do) / ho + (Di * dx) / (Dlm * kw) + 1 / hi + Rfi + (Di / Do) * Rfo)
    Qd = Ai * Ui_new * F12 * dTlm
    
    return Qr - Qd