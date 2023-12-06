function circulant_matrix = CreateUpperToeplitz(column)
    n = length(column);
    circulant_matrix = zeros(n);

    for i = 1:n
        for j = 1:n
            circulant_matrix(i, j) = column(mod((i - j), n) + 1);
        end
    end


    for i = 1:n
        for j = 1:n
            if i>j
            circulant_matrix(i, j) = inf;
            end
        end
    end
end
