function generator_matrix = GeneratorsToeplitz(A, alpha, isUpper)
    n = length(A);
    generator_matrix = zeros(n);

    if isUpper
        for i = 1:n
            for j = 1:n
                if i <= j && mod((i - j), n) == mod((n - alpha), n) && alpha<n
                    generator_matrix(i, j) = 0;
                else
                    generator_matrix(i, j) = Inf;
                end
            end
        end
    else  
        for i = 1:n
            for j = 1:n
                if i >= j && mod((i - j), n) == alpha
                    generator_matrix(i, j) = 0;
                else
                    generator_matrix(i, j) = Inf;
                end
            end
        end
    end
end
