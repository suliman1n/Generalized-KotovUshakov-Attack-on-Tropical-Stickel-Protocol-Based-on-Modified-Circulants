function [S, c] = Find_c_S(A, B, U, d)
    
    S = cell(d + 1, d + 1);
    c = zeros(d + 1, d + 1);

    for alpha = 0:d
        for beta = 0:d
            T = minplusMulti(minplusMatPower(A,alpha), minplusMatPower(B,beta)) - U;

            c_ab = min(T(:));
            S_ab = MinIndices(T);

            S{alpha+1,beta+1} = S_ab;
            c(alpha+1,beta+1) = c_ab;
        end
    end
end