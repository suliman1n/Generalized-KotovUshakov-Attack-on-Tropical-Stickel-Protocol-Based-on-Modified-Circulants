function circulant_matrix = CreateModifiedSCirculant(column, s)
    n = length(column);
    circulant_matrix = zeros(n);

    for i = 1:n
        for j = 1:n
            circulant_matrix(i, j) = column(mod((i - j), n) + 1);
        end
    end

    for i = 1:n
        for j = (i + 1):n
            circulant_matrix(i, j) = circulant_matrix(i, j) + s;
        end
    end
end
