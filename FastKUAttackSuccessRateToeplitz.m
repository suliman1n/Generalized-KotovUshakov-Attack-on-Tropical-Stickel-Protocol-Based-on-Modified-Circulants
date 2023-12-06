function successRate = FastKUAttackSuccessRateToeplitz(n,mm, mM,trials)
    
    d=n-1;
    successfulTrials = 0;

    for trial = 1:trials
        [key,U,V] = GenerateKeyStickelsToeplitz(n, mm, mM);

        [S, c] = Find_c_S_Toeplitz(minplusTropId(n), minplusTropId(n), U);

        
        cover = findMinimumCoverCirculant(S);
         
        [x, y] = SolveSystemGivenCover(cover, c);

        
        if ~isempty(x)
      
            K_attack = minplusMulti(ApplyGeneratorsToeplitzMinPlus([0:d; x'], minplusTropId(n),1),minplusMulti(V,ApplyGeneratorsToeplitzMinPlus([0:d; y'], minplusTropId(n),0)));


            
            if isequal(K_attack, key)
                successfulTrials = successfulTrials + 1;
            end
        end
    end

    
    successRate = successfulTrials / trials;
end
