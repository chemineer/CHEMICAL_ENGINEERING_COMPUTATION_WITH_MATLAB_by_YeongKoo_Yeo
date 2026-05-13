import numpy as np

def maxstep(delfun, n, m, nc, x, x0, d, A, b, nca):
    """
    제약 조건 기반 최대 단계 크기 결정
    :param delfun: 목적 함수의 그래디언트 함수
    :param n: 변수의 개수
    :param m: 제약 조건의 총 개수와 관련된 인덱스 범위 상한
    :param nc: 비활성 제약 조건 시작 인덱스
    :param x: 현재 지점
    :param x0: 초기 지점
    :param d: 탐색 방향
    :param A: 제약 조건 행렬
    :param b: 제약 조건 상수 벡터
    :param nca: 제약 조건 인덱스 배열
    :return: alphak (최대 단계 크기)
    """
    nq = 0
    alphak = 0.0
    
    # 제약 조건을 만족하는 최대 이동 거리 계산
    for k in range(nc, m):
        cq = int(nca[k])
        c = -b[cq]
        aq = 0.0
        
        for j in range(n):
            c += A[cq, j] * x[j]
            aq += A[cq, j] * d[j]
            
        if aq != 0:
            am = -c / aq
            if am > 0:
                nq += 1
                if nq == 1:
                    alphak = am
                else:
                    if alphak > am:
                        alphak = am
                        
    # 제약 조건이 없는 경우 단계 크기 확장
    if nq == 0:
        alphak = 1.0
        while True:
            x_test = x0 + alphak * d
            # delfun(x)의 그래디언트와 방향 d의 내적 계산
            fp = np.dot(delfun(x_test), d)
            
            if fp > 0:
                break
            alphak *= 2
            
    return alphak