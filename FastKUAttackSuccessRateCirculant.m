function successRate = FastKUAttackSuccessRateCirculant(n,mm, mM,trials)
    
    d=n-1;
    successfulTrials = 0;

    for trial = 1:trials
        s=randi([1 1000]);
        t=randi([1 1000]);
        [key,U,V] = GenerateKeyStickelsCirculant(n, mm, mM, s, t);

        [S, c] = Find_c_S_Circulant(minplusTropId(n), minplusTropId(n), U,s,t);

        
        cover = findMinimumCoverCirculant(S);

        
        [x, y] = SolveSystemGivenCover(cover, c);

        
        if ~isempty(x)
      
            K_attack = minplusMulti(ApplyGeneratorsMinPlus([0:d; x'], minplusTropId(n),s),minplusMulti(V,ApplyGeneratorsMinPlus([0:d; y'], minplusTropId(n),t)));


            
            if isequal(K_attack, key)
                successfulTrials = successfulTrials + 1;
            end
        end
    end

    
    successRate = successfulTrials / trials;
end
