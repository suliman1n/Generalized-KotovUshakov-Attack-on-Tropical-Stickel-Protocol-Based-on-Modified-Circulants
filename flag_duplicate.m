function f = flag_duplicate(S)
    f = zeros(size(S));
    for idx = 1:numel(S)
        is_unique = true;
        for other_idx = 1:numel(S)
            if other_idx ~= idx && isequal(S{idx}, S{other_idx})
                is_unique = false;
                break;
            end
        end
        if ~is_unique
            f(idx) = 1;
        end
    end
end
