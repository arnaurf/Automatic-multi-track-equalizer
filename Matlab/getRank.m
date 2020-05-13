function ranking = getRank(x, nF)
    ranking = zeros(nF, 1);
    for i = 1:nF
        [value, index] = max(x);
        x(index) = -Inf;
        ranking(index) = i;
    end

end
