function min_indices = MinIndices(A)

    min_value = min(A(:));

    
    [row_indices, col_indices] = find(A == min_value);

    min_indices = cell(length(row_indices), 1);

    for i = 1:length(row_indices)
        min_indices{i} = [row_indices(i), col_indices(i)];
    end
end
