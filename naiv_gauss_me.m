function x = naiv_gauss_me(A,b)
% 벡터화 버전(부분 피벗팅 없음): 내부 i-루프 제거
n = length(b);
x = zeros(n,1);

for k = 1:n-1
    if A(k,k) == 0
        error('zero pivot at k=%d (피벗팅 없이 진행 불가)', k);
    end
    invAk = 1 / A(k,k);
    xmult = A(k+1:n,k) * invAk;                       % (n-k)×1
    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - xmult * A(k,k+1:n);  % 랭크-1 업데이트
    b(k+1:n)       = b(k+1:n)       - xmult * b(k);
    A(k+1:n,k) = 0;                                   % 선택: 0으로 정리
end

x(n) = b(n)/A(n,n);
for i = n-1:-1:1
    x(i) = (b(i) - A(i,i+1:n)*x(i+1:n)) / A(i,i);     % 후진대입
end
end
