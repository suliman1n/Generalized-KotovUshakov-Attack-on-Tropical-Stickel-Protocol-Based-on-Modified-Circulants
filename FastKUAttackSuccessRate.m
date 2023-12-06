function successRate = FastKUAttackSuccessRate(n, d, mm, mM, pm, pM, trials)
    
    successfulTrials = 0;

    for trial = 1:trials
        [key,U,V,A,B] = GenerateKeyStickels(n, mm, mM, d, pm, pM);

        [S, c] = Find_c_S(A, B, U, d);
        
        cover = findMinimumCover(S,n,d);
 
        [x, y] = SolveSystemGivenCover(cover, c);

        
        if ~isempty(x)
            K_attack = minplusMulti(ApplyPolynomialMinPlus([0:d; x'], A), minplusMulti(V, ApplyPolynomialMinPlus([0:d; y'], B)));

            
            if isequal(K_attack, key)
                successfulTrials = successfulTrials + 1;
            end
        end
    end

    
    successRate = successfulTrials / trials;
end