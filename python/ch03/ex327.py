# -*- coding: utf-8 -*-
from gfree import gfree

# 온도 설정 (섭씨 60도를 절대온도 켈빈으로 변환)
T = 60 + 273.15

# 깁스 자유 에너지 계산 및 출력
# 매틀랩 코드에서 gfree(T, 'Benzene')을 호출하므로 동일하게 구성합니다.
gf=gfree(T, 'Benzene')
# 결과 출력
print(f"온도 T = {T}K 일 때, 벤젠의 자유에너지 (gf) = {gf}")