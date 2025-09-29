function indices = MISC(N)

P = 2 * floor(N/4) + 2;

InterSpace = [1 P-3 P*ones(1,N-P) 2*ones(1,(P-4)/2) 3 2*ones(1,(P-4)/2)];

indices = [0 cumsum(InterSpace)];

end