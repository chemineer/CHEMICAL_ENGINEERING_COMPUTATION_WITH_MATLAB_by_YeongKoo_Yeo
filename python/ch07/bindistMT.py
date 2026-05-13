import numpy as np
import matplotlib.pyplot as plt

def bindistMT(alpha, q, zf, xd, xb, R):
    """
    McCabe/Thiele 방법을 이용한 이성분 증류 계산
    alpha: 상대 휘발도
    q: 원료 공급 상태 파라미터
    zf, xd, xb: 원료, 탑상, 탑저 제품의 몰 분율
    R: 환류비
    """
    # 1. 초기화 및 평형 곡선 계산
    ye = np.linspace(0, 1, 100)
    xe = ye / (alpha + (1 - alpha) * ye)
    
    # 원료선(q-line)과 농축 조작선의 교점 (xq, yq)
    xq = ((R + 1) * zf + (q - 1) * xd) / (R + q)
    yq = (R * zf + q * xd) / (R + q)
    
    plt.figure(figsize=(8, 8))
    plt.plot(xe, ye, 'r', label='Equilibrium Curve') # 평형 곡선
    plt.plot([0, 1], [0, 1], 'k') # y=x 선
    
    # 조작선 그리기 (q-line, rectifying, stripping)
    plt.plot([xd, xq], [xd, yq], 'm') # 농축부 조작선 일부
    plt.plot([zf, xq], [zf, yq], 'm') # 원료선 (q-line)
    plt.plot([xb, xq], [xb, yq], 'm') # 회수부 조작선 일부
    
    # 계단 그리기 변수 초기화
    x_steps = [xd]
    y_steps = [xd]
    curr_y = xd
    curr_x = xd
    i = 0
    
    # 2. 농축부 (Rectifying section)
    while curr_x > xq:
        # 수평 이동: 평형 곡선 도달
        next_x = curr_y / (alpha + (1 - alpha) * curr_y)
        plt.hlines(y=curr_y, xmin=next_x, xmax=curr_x, colors='b')
        
        curr_x = next_x
        # 수직 이동: 조작선 도달 (xq보다 클 때만 농축부 조작선 사용)
        if curr_x > xq:
            next_y = R * curr_x / (R + 1) + xd / (R + 1)
            plt.vlines(x=curr_x, ymin=next_y, ymax=curr_y, colors='b')
            curr_y = next_y
            i += 1
            feedn = i
        else:
            # 원료 단을 지남
            break

    # 3. 회수부 (Stripping section)
    # 회수부 조작선 기울기 및 절편 계산
    c1 = (yq - xb) / (xq - xb)
    c2 = (yq - xq) / (xq - xb) # MATLAB의 c2는 (yq-xq)/(xq-xb)로 계산됨
    
    # 원료 단에서의 수직 하강 (농축부에서 회수부 조작선으로 전환)
    next_y = c1 * curr_x - (c1 * xb - xb) # y - xb = c1(x - xb) 형태
    plt.vlines(x=curr_x, ymin=next_y, ymax=curr_y, colors='b')
    curr_y = next_y
    i += 1

    while curr_x > xb:
        # 수평 이동
        next_x = curr_y / (alpha + (1 - alpha) * curr_y)
        plt.hlines(y=curr_y, xmin=next_x, xmax=curr_x, colors='b')
        curr_x = next_x
        
        if curr_x > xb:
            # 수직 이동: 회수부 조작선
            next_y = c1 * curr_x - (c1 * xb - xb)
            plt.vlines(x=curr_x, ymin=next_y, ymax=curr_y, colors='b')
            curr_y = next_y
            i += 1
        else:
            break
            
    totaln = i
    
    plt.xlabel('x (Liquid Mole Fraction)')
    plt.ylabel('y (Vapor Mole Fraction)')
    plt.title('McCabe-Thiele Method')
    plt.grid(True)
    plt.axis([0, 1, 0, 1])
    
    print(f"Feed stage = {feedn}")
    print(f"Number of stages = {totaln}")
    
    plt.show()
    return feedn, totaln

# 사용 예시:
# feedn, totaln = bindistMT(alpha=2.5, q=1, zf=0.5, xd=0.9, xb=0.1, R=2)