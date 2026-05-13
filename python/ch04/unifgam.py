import numpy as np

def unifgam(k, R, Q, nu, amn, Nc, x, T):
    """
    UNIFAC 방법을 사용하여 활동도 계수를 추정합니다.
    
    입력:
    k: 기능기(functional group)의 개수
    R, Q: 각 기능기의 부피 및 표면적 (벡터)
    nu: 각 성분에 포함된 기능기의 수 (행: 기능기 k, 열: 성분 i)
    amn: 그룹 상호작용 매개변수 행렬 (k x k)
    Nc: 성분 수
    x: 액체상 몰 분율 (벡터)
    T: 온도 (K)
    
    출력:
    gam: 각 성분에 대한 활동도 계수 (벡터)
    """
    R = np.array(R)
    Q = np.array(Q)
    nu = np.array(nu)
    amn = np.array(amn)
    x = np.array(x)
    
    # 1. 분자 부피 및 표면적 계산 (r_i, q_i)
    r = np.sum(R[:, np.newaxis] * nu, axis=0)
    q = np.sum(Q[:, np.newaxis] * nu, axis=0)
    
    # 2. 그룹 분율 계산
    # ek: 성분 i 내의 그룹 k의 표면적 분율
    ek = (nu * Q[:, np.newaxis]) / q
    
    # 3. 상호작용 매개변수 (tau)
    tau = np.exp(-amn / T)
    
    # 4. 혼합물 내 그룹 매개변수 (beta)
    # beta[i, j] = sum_k (ek[k, i] * tau[k, j])
    beta = ek.T @ tau
    
    # 5. 혼합물 내 평균 그룹 표면적 분율 (theta) 및 s
    theta = (x * q * ek) / np.sum(x * q)
    theta = np.sum(theta, axis=1) # 각 그룹별 평균
    
    s = theta @ tau
    
    # 6. 조합적(Combinatorial) 기여분 (gamc)
    J = r / np.sum(r * x)
    L = q / np.sum(q * x)
    gamc = 1 - J + np.log(J) - 5 * q * (1 - J / L + np.log(J / L))
    
    # 7. 잔류(Residual) 기여분 (gamr)
    gamr = np.zeros(Nc)
    for i in range(Nc):
        # sumb = sum_j (theta_j * beta_ij / s_j - ek_ji * ln(beta_ij / s_j))
        term1 = (theta * beta[i, :]) / s
        term2 = ek[:, i] * np.log(beta[i, :] / s)
        sumb = np.sum(term1 - term2)
        gamr[i] = q[i] * (1 - sumb)
        
    # 8. 최종 활동도 계수
    gam = np.exp(gamc + gamr)
    
    return gam