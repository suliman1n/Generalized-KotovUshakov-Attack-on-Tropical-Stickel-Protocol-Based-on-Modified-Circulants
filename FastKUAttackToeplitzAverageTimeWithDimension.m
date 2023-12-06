function FastKUAttackToeplitzAverageTimeWithDimension(max_n, mm, mM, trials)
    averageTimes = zeros(1, max_n);

    for n = 1:max_n
        totalElapsedTime = 0;
        d=n-1;

        for trial = 1:trials
            [key,U,V] = GenerateKeyStickelsToeplitz(n, mm, mM);
            tic;

            [S, c] = Find_c_S_Toeplitz(minplusTropId(n), minplusTropId(n), U);
            
            
            cover = findMinimumCoverCirculant(S);
  
            [x, y] = SolveSystemGivenCover(cover, c);

            
            if ~isempty(x)
                K_attack = minplusMulti(ApplyGeneratorsToeplitzMinPlus([0:d; x'], minplusTropId(n),1),minplusMulti(V,ApplyGeneratorsToeplitzMinPlus([0:d; y'], minplusTropId(n),0)));

                
                if isequal(K_attack, key)
                    totalElapsedTime = totalElapsedTime + toc;
                end
            end
        end

        
        averageTimes(d+1) = totalElapsedTime / trials;
    end

    
    figure;
    plot(1:max_n, averageTimes, '-o');
    title('Average Time vs Dimension');
    xlabel('Dimension (n)');
    ylabel('Average Time (seconds)');
    grid on;
end
