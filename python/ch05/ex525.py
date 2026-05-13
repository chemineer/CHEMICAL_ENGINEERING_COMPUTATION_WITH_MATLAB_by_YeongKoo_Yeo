import numpy as np
from twophreg import twophreg

# 하위 유동 양식 함수들 정의
def strat(x, Wa):
    rind = 1
    Yg = (15400 * x / (Wa**0.8))**2
    return Yg, rind

def wave(x, Wa):
    rind = 2
    Yg = 0
    return Yg, rind

def plug(x, Wa):
    rind = 3
    Yg = (27.315 * (x**0.855) / (Wa**0.17))**2
    return Yg, rind

def slug(x, Wa):
    rind = 4
    Yg = (1190 * (x**0.815) / np.sqrt(Wa))**2
    return Yg, rind

def bubb(x, Wa):
    rind = 5
    Yg = (14.2 * (x**0.75) / (Wa**0.1))**2
    return Yg, rind

def annul(x, d):
    rind = 6
    dx = d
    if d > 12:
        dx = 10
    Yg = ((4.8 - 0.3125 * dx) * x**(0.343 - 0.021 * dx))**2
    return Yg, rind

def dispr(x):
    rind = 7
    a0, a1, a2, a3 = 1.4659, 0.49138, 0.04887, -0.000349
    log_x = np.log(x)
    Yg = (np.exp(a0 + a1 * log_x + a2 * (log_x**2) + a3 * (log_x**3)))**2
    return Yg, rind

def twophmod(rhol, rhog, mul, sigl, Wl, Wg, delPl, delPg, d):
    D = d / 12
    A = np.pi * D**2 / 4
    # Baker parameter Bx calculation
    Bx = 531 * (Wl / Wg) * (np.sqrt(rhol * rhog) / (rhol**(2/3))) * ((mul**(1/3)) / sigl)
    By = 2.16 * (Wg / A) / np.sqrt(rhol * rhog)
    
    # twophreg.py에서 C1~C6 가져오기
    C = twophreg(Bx)
    
    x = np.sqrt(delPl / delPg)
    Wa = Wl / A
    
    # 유동 양식 분류 로직 (MATLAB 원본 구조 유지)
    if By <= C[0]: # MATLAB C(1)
        if By <= C[1]: # MATLAB C(2)
            Yg, rind = strat(x, Wa)
        else:
            Yg, rind = wave(x, Wa)
    else:
        if By < C[4]: # MATLAB C(5)
            if By < C[5]: # MATLAB C(6)
                Yg, rind = plug(x, Wa)
            else:
                if By < C[3]: # MATLAB C(4)
                    Yg, rind = slug(x, Wa)
                else:
                    if By <= C[2]: # MATLAB C(3)
                        Yg, rind = annul(x, d)
                    else:
                        Yg, rind = dispr(x)
        else:
            if Bx > 150:
                Yg, rind = bubb(x, Wa)
            else:
                if By <= C[2]: # MATLAB C(3)
                    Yg, rind = annul(x, d)
                else:
                    Yg, rind = dispr(x)
    return Yg, rind

# --- 메인 실행 부분 ---
# 데이터 설정
rhol, rhog, mul, mug, sigl = 66.7, 2.98, 1.0, 0.02, 70.0
L, d, Wl, Wg = 26400.0, 6.065, 77956.0, 12434.0
D = d / 12
eD = 0.00015 / D

# 마찰 계수(friction factor) 계산
Nrel = 6.31 * Wl / (d * mul)
Nreg = 6.31 * Wg / (d * mug)

if Nrel <= 2100:
    fL = 64 / Nrel
else:
    Av = eD / 3.7 + (6.7 / Nrel)**0.9
    fL = 4.0 / (-4 * np.log10(eD / 3.7 - 5.02 * np.log10(Av) / Nrel))**2

if Nreg <= 2100:
    fG = 64 / Nreg
else:
    Av = eD / 3.7 + (6.7 / Nreg)**0.9
    fG = 4.0 / (-4 * np.log10(eD / 3.7 - 5.02 * np.log10(Av) / Nreg))**2

fr = ['stratified', 'wave', 'plug', 'slug', 'bubble', 'annular', 'dispersed']

delPl = 3.66e-4 * fL * Wl**2 / (d**5 * rhol)
delPg = 3.66e-4 * fG * Wg**2 / (d**5 * rhog)

# 유동 양식 및 모듈러스 계산
Yg, rind = twophmod(rhol, rhog, mul, sigl, Wl, Wg, delPl, delPg, d)

if rind == 2: # wave flow인 경우 별도 계산
    fH = np.exp(0.2111 * np.log(Wl * mul / (Wg * mug)) - 3.993)
    delPt = 3.66e-4 * fH * Wg**2 * L / (d**5 * rhog * 100)
else:
    delPt = Yg * delPg * L / 100

# 결과 출력
print(f"Flow regime: {fr[rind-1]}")
print(f"Two-phase flow modulus (Yg): {Yg:.6f}")
print(f"Total pressure drop(psi): {delPt:.6f}")