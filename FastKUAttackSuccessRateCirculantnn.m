function successRate = FastKUAttackSuccessRateCirculantnn(n,mm, mM,trials)
    
    d=n-1;
    successfulTrials = 0;

    for trial = 1:trials
        s=randi([mm mM]);
        t=randi([mm mM]);
        [key,U,V,W] = GenerateKeyStickelsCirculant(n, mm, mM, s, t);

        [S, c] = Find_c_S_Circulant1(minplusTropId(n), minplusTropId(n),W, U,s,t);

        cover = findMinimumCover1(S);
        
        globalVariables=[];
        for k=find(sum(cat(n,S{:}),n) == 1)'
            for index = 1:numel(S)
                if S{index}(k)==1
                    if ~ismember(index, globalVariables)
                        globalVariables=[globalVariables,index];
                    end
                end  
            end
        end
        
        

        [dPlusOne, ~] = size(S);
        linearIndex = globalVariables;
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
