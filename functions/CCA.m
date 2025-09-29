function [pos , eleNum] = CCA(k,M,N)

if N>M && gcd(M,N) == 1
    ktimesExtSubArray     = 0:M:(N-1)*M;
    CoprimeSubArray       = 0:N:(k*M-1)*N;
    CompSubArray          = (k*M-1)*N-1:-1:(k*M-1)*N-(M-1);
    pos                   = unique([ktimesExtSubArray CoprimeSubArray CompSubArray]);

else
    % errordlg('Choose N greater than M or M&N are not co-prime integers','Warning')
    % return
    pos = [];
    % eleNum = [];

end

eleNum = (k+1)*M+N-2;

end


