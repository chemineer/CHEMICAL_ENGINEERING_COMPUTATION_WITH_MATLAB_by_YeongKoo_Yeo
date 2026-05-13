import numpy as np
from scipy.optimize import fsolve
from phiHeu import phiHeu  # phiHeu.py 모듈이 있다고 가정

def distBPeu(opdat, mxdat):
    """
    BP 방법을 이용한 다성분 증류 계산 (영국 단위계)
    """
    # 데이터 설정
    eos = opdat['eos']
    N = opdat['N']
    nc = opdat['nc']
    F = opdat['F']
    Tf = opdat['Tf']
    Pf = opdat['Pf']
    P = opdat['P']
    V = opdat['V'].copy()
    L = opdat['L'].copy()
    criv = opdat['criv']
    U = opdat['U']
    W = opdat['W']
    Q = opdat['Q'].copy()
    T = opdat['T0'].copy()
    nf = opdat['nf']
    z = opdat['z']
    
    Pc = mxdat['Pc']
    Ant = mxdat['Ant']
    
    hF = np.zeros(N)
    hV = np.zeros(N)
    hL = np.zeros(N)
    x = np.zeros((nc, N))
    y = np.zeros((nc, N))
    K = np.zeros((nc, N))
    
    # 초기화: V(2) 설정
    # 파이썬 인덱스 주의 (MATLAB V(2) -> Python V[1])
    V[1] = L[0] + V[0] + U[0] + W[0] - F[0]

    # 피드 엔탈피 계산
    fstate_upper = opdat['fstate'].upper()
    for j_val in nf:
        # nf가 MATLAB 인덱스라면 1을 빼서 사용 (이미 0-index로 제공되었다고 가정)
        Z, H, phi = phiHeu(z[:, j_val], Pf[j_val], Tf[j_val], fstate_upper, eos, opdat, mxdat)
        hF[j_val] = H

    Told = T.copy()
    criT = 10.0
    iter_count = 1

    # 초기 K(i,j) 계산 (Raoult의 법칙 기반)
    for j in range(N):
        pv0 = Pc * np.exp(Ant[:, 0] - Ant[:, 1] / (T[j] + Ant[:, 2]))
        K[:, j] = pv0 / P[j]

    # 반복 계산 시작
    while criT > criv:
        # 1. 조성 x(i,j) 계산 (Thomas Algorithm/Tridiagonal Matrix Solver)
        for i in range(nc):
            Pj = np.zeros(N)
            Qj = np.zeros(N)
            Rj = np.zeros(N)
            Sj = np.zeros(N)
            
            # 중간 단 방정식 계수 설정
            for j in range(1, N-1):
                Pj[j] = V[j] - V[0] + np.sum(F[:j] - U[:j] - W[:j])
                Qj[j] = V[0] - (V[j] + W[j]) * K[i, j] - U[j] - V[j+1] - np.sum(F[:j+1] - U[:j+1] - W[:j+1])
                Rj[j] = V[j+1] * K[i, j+1]
                Sj[j] = -F[j] * z[i, j]
            
            # 첫 번째 단과 마지막 단 계수
            Qj[0] = V[0] - (V[0] + W[0]) * K[i, 0] - U[0] - V[1] - (F[0] - U[0] - W[0])
            Rj[0] = V[1] * K[i, 1]
            Sj[0] = -F[0] * z[i, 0]
            
            Pj[N-1] = V[N-1] - V[0] + np.sum(F[:N-1] - U[:N-1] - W[:N-1])
            Qj[N-1] = V[0] - (V[N-1] + W[N-1]) * K[i, N-1] - U[N-1] - np.sum(F - U - W)
            Sj[N-1] = -F[N-1] * z[i, N-1]
            
            # TDMA 풀이
            r_diag = np.zeros(N)
            s_diag = np.zeros(N)
            r_diag[0] = Rj[0] / Qj[0]
            s_diag[0] = Sj[0] / Qj[0]
            
            for j in range(1, N-1):
                denom = Qj[j] - r_diag[j-1] * Pj[j]
                r_diag[j] = Rj[j] / denom
                s_diag[j] = (Sj[j] - s_diag[j-1] * Pj[j]) / denom
            
            s_diag[N-1] = (Sj[N-1] - s_diag[N-2] * Pj[N-1]) / (Qj[N-1] - r_diag[N-2] * Pj[N-1])
            x[i, N-1] = s_diag[N-1]
            
            for j in range(N-2, -1, -1):
                x[i, j] = s_diag[j] - r_diag[j] * x[i, j+1]

        # 정규화 (Normalization)
        for j in range(N):
            x[:, j] = x[:, j] / np.sum(x[:, j])

        # 2. BP 방법을 이용한 새로운 T(j) 및 y(i,j) 계산
        for j in range(N):
            y[:, j] = K[:, j] * x[:, j]
            y[:, j] = y[:, j] / np.sum(y[:, j])

        # 버블점 온도 T 계산 (fsolve 사용)
        for j in range(N):
            def bubble_point_obj(Tv):
                pv = Pc * np.exp(Ant[:, 0] - Ant[:, 1] / (Tv + Ant[:, 2]))
                return np.sum(pv * x[:, j] / P[j]) - 1.0
            
            T[j] = fsolve(bubble_point_obj, T[j])[0]

        # K 값 업데이트
        for j in range(N):
            _, hL_j, phiL = phiHeu(x[:, j], P[j], T[j], 'L', eos, opdat, mxdat)
            hL[j] = hL_j
            _, hV_j, phiV = phiHeu(y[:, j], P[j], T[j], 'V', eos, opdat, mxdat)
            hV[j] = hV_j
            K[:, j] = phiL / phiV

        # 3. 콘덴서 및 리보일러 열부하 계산
        Q[0] = (L[0] + V[0] + U[0] + W[0] - F[0]) * hV[1] - V[0] * hV[0] - (L[0] + U[0]) * hL[0] + F[0] * hF[0]
        Q[N-1] = np.sum(F * hF - U * hL - W * hV) - np.sum(Q[:N-1]) - V[0] * hV[0] - L[N-1] * hL[N-1]

        # 4. 유량 V(j) 및 L(j) 계산
        a_vec = np.zeros(N)
        b_vec = np.zeros(N)
        d_vec = np.zeros(N)
        
        for j in range(1, N):
            a_vec[j] = hL[j-1] - hV[j]
            b_vec[j-1] = hV[j] - hL[j-1]
            d_vec[j] = (hL[j] - hL[j-1]) * np.sum(F[:j] - U[:j] - W[:j]) + \
                       F[j] * (hL[j] - hF[j]) + W[j] * (hV[j] - hL[j]) + Q[j]
        
        a_vec[0] = -hV[0]
        b_vec[N-1] = -hL[N-1]
        d_vec[0] = F[0] * (hL[0] - hF[0]) + W[0] * (hV[0] - hL[0]) + Q[0]
        
        for j in range(1, N-1):
            V[j+1] = (d_vec[j] - a_vec[j] * V[j]) / b_vec[j]
            L[j] = V[j+1] - V[0] + np.sum(F[:j+1] - U[:j+1] - W[:j+1])

        # 수렴 확인
        criT = np.sum(np.abs(Told - T))
        Told = T.copy()
        iter_count += 1

    return x, y, T, L, V, iter_count