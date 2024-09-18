function successRate = FastKUAttackSuccessRateCirculantn(n,mm, mM,trials)
    
    d=n-1;
    successfulTrials = 0;

    for trial = 1:trials
        s=randi([0 mM]);
        t=randi([0 mM]);
        [key,U,V,W] = GenerateKeyStickelsCirculant(n, mm, mM, s, t);

        [S, c] = Find_c_S_Circulant1(minplusTropId(n), minplusTropId(n),W, U,s,t);

        
        f=flag_duplicate(S);
        non_unique_indices = find(f);
        
        cover = findMinimumCover1(S);

        cover = setdiff(cover, non_unique_indices);
        
        

        [dPlusOne, ~] = size(S);
        linearIndex = cover;
        [alpha, beta] = ind2sub([dPlusOne, dPlusOne], linearIndex);
        cover_alpha_beta = [alpha', beta'];

        
        [x, y] = SolveSystemGivenCover(cover_alpha_beta, c);

        
        if ~isempty(x)
      
            K_attack = minplusMulti(ApplyGeneratorsMinPlus([0:d; x'], minplusTropId(n),s),minplusMulti(V,ApplyGeneratorsMinPlus([0:d; y'], minplusTropId(n),t)));


            
            if isequal(K_attack, key)
                successfulTrials = successfulTrials + 1;
            end
        end
    end

    
    successRate = successfulTrials / trials;
end
