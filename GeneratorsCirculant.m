function generator_matrix = GeneratorsCirculant(A, alpha, s)
    n = length(A);
    generator_matrix = zeros(n);

    for i = 1:n
        for j = 1:n
            if i >= j && mod((i - j), n) == alpha
                generator_matrix(i, j) = 0;
            elseif i < j && mod((i - j), n) == alpha
                generator_matrix(i, j) = s;
            else
                generator_matrix(i, j) = Inf;
            end
        end
    end
end