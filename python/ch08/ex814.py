import numpy as np
from scipy.optimize import fsolve
import math
from heatH2Odat import heatH2Odat  # 데이터 정의 함수 임포트
from heatH2O import heatH2O        # 잔차 계산 함수 임포트

# 1. 초기 데이터 로드
data = heatH2Odat()

# 2. fsolve를 이용한 배플 간격(b) 결정
# MATLAB의 b0 = 0.0 대신, 계산 안정성을 위해 0.5(ft)를 초기값으로 설정합니다.
b0 = 0.5 
b_sol = fsolve(heatH2O, b0)[0]

# 3. 최종 결과 산출을 위한 재계산 (ex814.m 로직 반영)
b = b_sol
Qr, Ui_init, F12, dTlm = data['Qr'], data['Ui_init'], data['F12'], data['dTlm']
At, Nt, Di, Do, L, Pt, Ds, cl = data['At'], data['Nt'], data['Di'], data['Do'], data['L'], data['Pt'], data['Ds'], data['cl']
mh, rhoc, rhoh, ui, gc = data['mh'], data['rhoc'], data['rhoh'], data['ui'], data['gc']
muh, muwh, kh, Cph = data['muh'], data['muwh'], data['kh'], data['Cph']
muc, muwc, kc, Cpc = data['muc'], data['muwc'], data['kc'], data['Cpc']
dx, Dlm, kw, Rfi, Rfo, De = data['dx'], data['Dlm'], data['kw'], data['Rfi'], data['Rfo'], data['De']

# 패스 수 및 면적 확정
Ai_temp = Qr / (Ui_init * F12 * dTlm)
Np = math.ceil(Ai_temp / (At * Nt))
Acf = Ds * cl * b / Pt
Go = mh / Acf

# 쉘 측(Shell-side) 계수 계산
Nreo = De * Go / (muh * 3600 / 1488)
Npro = Cph * (muh * 3600 / 1488) / kh
ho = 0.36 * Nreo**0.55 * Npro**(1/3) * (muh / muwh)**0.14 * kh / De

# 튜브 측(Tube-side) 계수 계산
Nrei = Di * rhoc * ui / (muc / 1488)
Npri = Cpc * (muc * 3600 / 1488) / kc
hi = 0.027 * Nrei**0.8 * Npri**(1/3) * (muc / muwc)**0.14 * kc / Di

# 총괄 열전달 계수 및 최종 면적
Ui = 1 / ((Di / Do) / ho + (Di * dx) / (Dlm * kw) + 1 / hi + Rfi + (Di / Do) * Rfo)
Uo = Ui * Di / Do
Ai = Qr / (Ui * F12 * dTlm)
Ao = Qr / (Uo * F12 * dTlm)

# 압력 강하(Pressure Drop) 계산
xt = Pt / Do
xl = xt
Ks = 1.1 * L / b
Nr = math.ceil(Ds / Pt / 2)
fp = (0.044 + 0.08 * xl / ((xt - 1)**(0.43 + 1.13 / xl))) * ((Do * Go) / muh)**(-0.15)
dPo = 2 * Ks * Nr * fp * (Go / 3600)**2 / (gc * rhoh * 144)

fD = 1 / (1.82 * math.log10(Nrei) - 1.64)**2
Gi = rhoc * ui
dPi = 0.6 * Np * fD * Gi**2 * L / (gc * rhoc * Di * 144)

# 4. 결과 출력 (MATLAB fprintf 형식 재현)
print(f"Number of tube passes = {Np}")
print(f"Baffle spacing = {b*12:.6f} in")
print(f"Shell-side overall heat transfer coefficient = {Uo:.6f} Btu/(ft^2*h*F)")
print(f"Tube-side overall heat transfer coefficient = {Ui:.6f} Btu/(ft^2*h*F)")
print(f"Tube inside heat transfer area = {Ai:.6f} ft^2")
print(f"Tube outside heat transfer area = {Ao:.6f} ft^2")
print(f"Shell-side pressure drop = {dPo:.6f} psi")
print(f"Tube-side pressure drop = {dPi:.6f} psi")