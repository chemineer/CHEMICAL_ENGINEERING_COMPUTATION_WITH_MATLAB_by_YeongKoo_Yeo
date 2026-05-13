import numpy as np
from scipy.optimize import fsolve

def frfactor(eD, Nre):
    """
    MATLAB frfactor.m의 상관식들을 파이썬으로 구현.
    수렴 안정성을 위해 Haaland 식의 결과를 초기값으로 활용합니다.
    """
    # eD가 0일 때 로그 오류 방지를 위한 안전 장치
    eD_safe = max(eD, 1e-15)
    
    # 4. Haaland eqn. (먼저 계산하여 다른 식의 초기값 f0로 활용)
    # fHA = 1./(-3.6*log10(6.9/Nre + (eD/3.7)^(10/9)))^2
    fHA = 1.0 / (-3.6 * np.log10(6.9/Nre + (eD_safe/3.7)**(10/9)))**2
    f0 = fHA # 수렴성이 검증된 값을 초기값으로 설정

    # 1. Shacham eqn.
    # fSC = 1./(log10(eD/3.7 - (5.02./Nre).*log10(eD/3.7+14.5/Nre))).^2 /16;
    term_sh1 = eD_safe / 3.7
    term_sh2 = (5.02 / Nre) * np.log10(term_sh1 + 14.5 / Nre)
    # 로그 내부가 양수가 되도록 처리
    fSC = 1.0 / (np.log10(np.abs(term_sh1 - term_sh2)))**2 / 16.0

    # 2. Colebrook eqn. (fsolve 사용)
    # funCB = @(f) (1/sqrt(f) + 1.7372*log(eD/3.7 + 1.255/Nre/sqrt(f)));
    def funCB(f):
        if f <= 0: return 1e6 # f가 음수일 때 solver가 튕기지 않도록 방어
        return (1.0/np.sqrt(f)) + 1.7372 * np.log(eD_safe/3.7 + 1.255/(Nre * np.sqrt(f)))
    fCB = fsolve(funCB, f0)[0]

    # 3. Colebrook-White eqn.
    # funCBW = @(f) (1/sqrt(f) + 4*log10(eD + 4.67/Nre/sqrt(f)) - 2.28);
    def funCBW(f):
        if f <= 0: return 1e6
        return (1.0/np.sqrt(f)) + 4 * np.log10(eD_safe + 4.67/(Nre * np.sqrt(f))) - 2.28
    fCBW = fsolve(funCBW, f0)[0]

    # 5. Chen eqn.
    Av = eD_safe/3.7 + (6.7/Nre)**0.9
    fCH = 1.0 / (-4 * np.log10(eD_safe/3.7 - 5.02 * np.log10(Av) / Nre))**2

    # 6. Nikuradse eqn.
    # funNK = @(f) (1./sqrt(f) - 4*log10(Nre*sqrt(f)) + 0.4);
    def funNK(f):
        if f <= 0: return 1e6
        return (1.0/np.sqrt(f)) - 4 * np.log10(Nre * np.sqrt(f)) + 0.4
    fNK = fsolve(funNK, f0)[0]

    # 7. Blasius eqn.
    fBS = 0.0791 * Nre**(-0.25)

    # MATLAB 형식 출력
    print(f"--- Results for eD={eD}, Nre={Nre:.1e} ---")
    print(f'Shacham eqn.: f = {fSC:.6g}')
    print(f'Colebrook eqn.: f = {fCB:.6g}')
    print(f'Colebrook-White eqn.: f = {fCBW:.6g}')
    print(f'Haaland eqn.: f = {fHA:.6g}')
    print(f'Chen eqn.: f = {fCH:.6g}')
    print(f'Nikuradse eqn.: f = {fNK:.6g}')
    print(f'Blasius eqn.: f = {fBS:.6g}\n')