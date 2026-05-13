import numpy as np
from scipy.optimize import fsolve
from coolhxdat import coolhxdat  # 데이터 정의 함수 임포트
from coolhxf import coolhxf      # 잔차 계산 함수 임포트

# 1. 데이터 로드
# coolhxdat()는 열교환기 설계에 필요한 물성치와 사양을 딕셔너리로 반환합니다.
data = coolhxdat()

# 2. fsolve를 이용한 배플 간격(b) 결정
# MATLAB의 fsolve(@coolhxf, b0)와 동일한 동작을 수행합니다.
# b0는 배플 간격의 초기 추정값입니다.
b0 = 0.5 
# coolhxf 함수는 (b, data)를 인자로 받으므로 args를 통해 data를 전달합니다.
b_solution = fsolve(coolhxf, b0, args=(data,))[0]

# 3. 최종 설계 결과 계산 (ex813.m의 후반부 로직 반영)
b = b_solution
Ds, cl, Pt, mc = data['Ds'], data['cl'], data['Pt'], data['mc']
De, muc, Cpc, kc, muwc = data['De'], data['muc'], data['Cpc'], data['kc'], data['muwc']
Di, rhoh, ui, muh, kh, muwh = data['Di'], data['rhoh'], data['ui'], data['muh'], data['kh'], data['muwh']
Do, dx, Dlm, kw, Rfi, Rfo = data['Do'], data['dx'], data['Dlm'], data['kw'], data['Rfi'], data['Rfo']
Q, F12, dTlm, At, Nt, gc, rhoc = data['Q'], data['F12'], data['dTlm'], data['At'], data['Nt'], data['gc'], data['rhoc']

# 쉘 측(Shell-side) 열전달 계수 계산
Acf = Ds * cl * b / Pt
Go = mc / Acf
Nreo = De * Go / (muc * 3600 / 1488)
Npro = Cpc * (muc * 3600 / 1488) / kc
Nuo = 0.36 * Nreo**0.55 * Npro**(1/3) * (muc / muwc)**0.14
ho = Nuo * kc / De

# 튜브 측(Tube-side) 열전달 계수 계산
Nrei = Di * rhoh * ui / (muh / 1488)
Npri = Cph = data['Cph'] * (muh * 3600 / 1488) / kh
Nui = 0.027 * Nrei**0.8 * Npri**(1/3) * (muh / muwh)**0.14
hi = Nui * kh / Di

# 총괄 열전달 계수 및 면적 계산
Ui = 1 / ((Di/Do)/ho + (Di*dx)/(Dlm*kw) + 1/hi + Rfi + (Di/Do)*Rfo)
Uo = Ui * Di / Do
Ai = Q / (Ui * F12 * dTlm)
Ao = Q / (Uo * F12 * dTlm)
Np = np.ceil(Ai / (At * Nt))

# 압력 강하(Pressure Drop) 계산
xt = Pt / Do
xl = xt
Ks = 1.1 * data['L'] / b
Nr = np.ceil(Ds / Pt / 2)
fp = (0.044 + 0.08 * xl / ((xt - 1)**(0.43 + 1.13 / xl))) * ((Do * Go) / muc)**(-0.15)
dPo = 2 * Ks * Nr * fp * (Go / 3600)**2 / (gc * rhoc * 144)

fD = 1 / (1.82 * np.log10(Nrei) - 1.64)**2
Gi = rhoh * ui
dPi = 0.6 * Np * fD * Gi**2 * data['L'] / (gc * rhoh * Di * 144)

# 4. 결과 출력 (MATLAB fprintf 재현)
print(f"Baffle spacing = {b*12:.6f} in")
print(f"Shell-side overall heat transfer coefficient = {Uo:.6f} Btu/(ft^2*h*F)")
print(f"Tube-side overall heat transfer coefficient = {Ui:.6f} Btu/(ft^2*h*F)")
print(f"Tube inside heat transfer area = {Ai:.6f} ft^2")
print(f"Tube outside heat transfer area = {Ao:.6f} ft^2")
print(f"Shell-side pressure drop = {dPo:.6f} psi")
print(f"Tube-side pressure drop = {dPi:.6f} psi")