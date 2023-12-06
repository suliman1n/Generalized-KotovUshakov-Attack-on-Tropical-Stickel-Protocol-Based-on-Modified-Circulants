function [S, c] = Find_c_S_Circulant(A, B, U,s,t)
    
    n = size(A, 1);
    d=n-1;
 
    S = cell(d + 1, d + 1);
    c = zeros(d + 1, d + 1);

    for alpha = 0:d
        for beta = 0:d
            T = minplusMulti(GeneratorsCirculant(A,alpha,s), GeneratorsCirculant(B,beta,t)) - U;

            c_ab = min(T(:));
            S_ab = MinIndices(T);

            S{alpha+1,beta+1} = S_ab;
            c(alpha+1,beta+1) = c_ab;
        end
    end
end

