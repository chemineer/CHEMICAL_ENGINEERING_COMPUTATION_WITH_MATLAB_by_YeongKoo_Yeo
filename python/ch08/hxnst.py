import math
import numpy as np
from hxndat import HeatExchangerData
from hxcp import hxcp
from hxvis import hxvis
from hxrho import hxrho
from hxthc import hxthc
from htcshell import htcshell
from htctube import htctube
from dptube import dptube
from dpshell import dpshell

def hxnst():
    # 데이터 로드
    d = HeatExchangerData()
    
    # 편의를 위해 변수 추출
    ptype = d.ptype
    Do, Di, Xkw = d.Do, d.Di, d.Xkw
    L, Ls, Nt = d.L, d.Ls, d.Nt
    Ti1, Ti2, Wi = d.Ti1, d.Ti2, d.Wi
    Ts1, Ts2, Ws = d.Ts1, d.Ts2, d.Ws
    cpreft, Trt, cprefs, Trs = d.cpreft, d.Trt, d.cprefs, d.Trs
    mureft, rhoreft, xkreft = d.mureft, d.rhoreft, d.xkreft
    murefs, rhorefs, xkrefs = d.murefs, d.rhorefs, d.xkrefs
    fsT, fsS, Npass = d.fsT, d.fsS, d.Npass
    Rds, Rdt, rf = d.Rds, d.Rdt, d.rf
    Ds, Lbc, Lbin, Lbout, Lc, Dotl, Dsb, Pt, Nss, Layout = (
        d.Ds, d.Lbc, d.Lbin, d.Lbout, d.Lc, d.Dotl, d.Dsb, d.Pt, d.Nss, d.Layout
    )

    # Nvar 결정
    Nvar = 0
    if ptype > 0:
        if ptype <= 3:
            Nvar = 1  # Ti1 계산
        elif ptype <= 6:
            Nvar = 2  # Ti2 계산
        else:
            Nvar = 3  # Wi 계산

    Rw = 1e-3 * (Do - Di) / (2 * Xkw)
    A = math.pi * 1e-3 * Do * (L - 2e-3 * Ls) * Nt

    if ptype == 0:
        # hxcp, hxvis 등 외부 함수는 환경에 맞게 정의되어 있다고 가정합니다.
        Cpim = hxcp(cpreft, Trt, (Ti1 + Ti2) / 2)
        Cpsm = hxcp(cprefs, Trs, (Ts1 + Ts2) / 2)
        Qi = Wi * Cpim * (Ti2 - Ti1)
        Qs = Ws * Cpsm * (Ts2 - Ts1)
        Qr = abs(Qi / Qs)
        if abs(1 - Qr) > 0.1:
            print('Tube-side and shell-side heat duty differ by more than 10%.')
            print('Tube-side heat duty is used as total duty.')
        Q = Qi

    varC = 10
    iter_count = 0
    while varC > 1e-3:
        if ptype > 0:
            Cpim = hxcp(cpreft, Trt, (Ti1 + Ti2) / 2)
            Cpsm = hxcp(cprefs, Trs, (Ts1 + Ts2) / 2)
            Qs = Ws * Cpsm * (Ts2 - Ts1)
            Qi = -Qs
            
            if Nvar == 1:
                Ti1new = Ti2 - Qi / (Wi * Cpim)
                crT = 10
                while crT >= 1e-3:
                    Ti1 = Ti1new
                    Cpim = hxcp(cpreft, Trt, (Ti1 + Ti2) / 2)
                    Ti1new = Ti2 - Qi / (Wi * Cpim)
                    crT = abs((Ti1new - Ti1) / Ti1new)
            elif Nvar == 2:
                Ti2new = Ti1 + Qi / (Wi * Cpim)
                crT = 10
                while crT >= 1e-3:
                    Ti2 = Ti2new
                    Cpim = hxcp(cpreft, Trt, (Ti1 + Ti2) / 2)
                    Ti2new = Ti1 + Qi / (Wi * Cpim)
                    crT = abs((Ti2new - Ti2) / Ti2new)
            elif Nvar == 3:
                Wi = abs(Qs / (Cpim * (Ti2 - Ti1)))

        Tib = (Ti1 + Ti2) / 2
        Tsb = (Ts1 + Ts2) / 2
        mus = hxvis(murefs, Trs, Tsb, fsS)
        rhos = hxrho(rhorefs, Trs, Tsb, fsS)
        Cps = hxcp(cprefs, Trs, Tsb)
        xks = hxthc(xkrefs, Trs, Tsb)
        
        Tw = (Tsb + Tib) / 2
        Twnew = Tw
        crT_inner = 10
        while crT_inner >= 1e-3:
            if fsS == 1:
                Phis = (mus / hxvis(murefs, Trs, Tw, 1))**0.14
            elif fsS == 2:
                Phis = (Tsb / Tw)**0.25
            
            Hs = htcshell(Do, Ds, L - 2e-3 * Ls, Lbc, Lbin, Lbout, Lc, Dotl, Dsb, Pt, Nt, Nss, Ws, mus, Cps, xks, Phis, Layout)
            
            mut = hxvis(mureft, Trt, Tib, fsT)
            rhot = hxrho(rhoreft, Trt, Tib, fsT)
            Cpt = hxcp(cpreft, Trt, Tib)
            xkt = hxthc(xkreft, Trt, Tib)
            
            Ui_vel = 4e6 * Wi * Npass / (math.pi * rhot * Di**2 * Nt)
            Rei = 1e-3 * Di * Ui_vel * rhot / mut
            Pri = Cpt * mut / xkt
            
            if fsT == 1:
                Phit = (mut / hxvis(mureft, Trt, Tw, 1))**0.14
            elif fsT == 2:
                Phit = (Tw / Tib)**0.25
                
            Ht = htctube(Rei, Pri, Di, L, xkt, Phit)
            Tw = Tib + Hs / (Hs + Ht) * (Tsb - Tib)
            crT_inner = abs((Tw - Twnew) / Tw)
            Twnew = Tw

        # 열전달 방정식 계산
        U = 1 / (Do / (Di * Ht) + 1 / Hs + Rw + Rds + Rdt)
        Dt1 = Ts1 - Ti2
        Dt2 = Ts2 - Ti1
        if Dt1 <= 0 or Dt2 <= 0:
            break
            
        # LMTD 계산 안전 장치 추가  
        if abs(Dt1 - Dt2) < 1e-5:
            Delt = Dt1
        else:
            Delt = (Dt1 - Dt2) / math.log(Dt1 / Dt2)
        Ft = 1
        if Npass > 1:
            R_val = (Ts1 - Ts2) / (Ti2 - Ti1)
            P_val = (Ti2 - Ti1) / (Ts1 - Ti1)
            tm = math.sqrt(R_val**2 + 1)
            # Ft 공식 (분모 0 체크 생략)
            num = (1 - P_val) / (1 - R_val * P_val)
            den1 = (2 - P_val * (R_val + 1 - tm)) / (2 - P_val * (R_val + 1 + tm))
            Ft = tm * math.log(num) / ((R_val - 1) * math.log(den1))
            
        Deltm = Delt * Ft
        if ptype == 0:
            Areq = Qi / (U * Deltm)
            Da = (A - Areq) / A * 100
            break
            
        Q = U * A * Deltm
        Sgn = -1 if Qs < 0 else 1
        
        if ptype in [1, 4, 7]:
            Ts1new = Ts2 - Sgn * Q / (Ws * Cps)
            varC = abs((Ts1 - Ts1new) / Ts1new)
            Ts1 = Ts1new
        elif ptype in [2, 5, 8]:
            Ts2new = Ts1 + Sgn * Q / (Ws * Cps)
            varC = abs((Ts2 - Ts2new) / Ts2new)
            Ts2 = Ts2new
        elif ptype in [3, 6]:
            Wsnew = abs(Q / ((Ts2 - Ts1) * Cps))
            varC = abs((Ws - Wsnew) / Wsnew)
            Ws = Wsnew
        
        iter_count += 1

    # 압력 강하 및 결과 출력
    DPs = dpshell(Do, Ds, L - 2e-3 * Ls, Lbc, Lbin, Lbout, Lc, Dotl, Dsb, Pt, Nt, Nss, Ws, rhos, mus, Phis, Layout)
    DPt = dptube(Di, L, rf, Ui_vel, rhot, mut, Phit, Npass)

    print(f"Overall heat transfer coefficient: U = {U:g}(W/m^2/K)")
    print(f"Heat transfer coefficient: tube-side = {Ht:g}(W/m^2/K), shell-side = {Hs:g}(W/m^2/K)")
    print(f"Heat duty: Q = {Q:g}(W)")
    print(f"Pressure drop: tube-side = {DPt:g}(Pa), shell-side = {DPs:g}(Pa)")
    print(f"Tube-side: Ti1 = {Ti1:g}(K), Ti2 = {Ti2:g}(K), flow rate = {Wi:g}(kg/sec)")
    print(f"Shell-side: Ts1 = {Ts1:g}(K), Ts2 = {Ts2:g}(K), flow rate = {Ws:g}(kg/sec)")
# 참고: htcshell, htctube, dpshell, dptube 함수는 별도로 정의되어 있어야 합니다.

# 테스트 실행 
if __name__ == "__main__":
    hxnst()