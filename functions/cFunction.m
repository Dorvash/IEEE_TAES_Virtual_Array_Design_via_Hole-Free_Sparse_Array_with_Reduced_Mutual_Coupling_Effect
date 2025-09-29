function ci = cFunction(dist,B)

ci = 1 * (dist == 0) + (0.5 ./ dist).^(1/1.6) .* (dist <= B & dist > 0) + 0 * (dist > B);
ci(dist == 0) = 1;

end