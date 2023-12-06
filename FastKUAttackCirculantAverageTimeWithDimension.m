function FastKUAttackCirculantAverageTimeWithDimension(max_n, mm, mM, trials)
    averageTimes = zeros(1, max_n);

    for n = 1:max_n
        totalElapsedTime = 0;
        d=n-1;

        for trial = 1:trials
            s=randi([1 1000]);
            t=randi([1 1000]);
            [key,U,V] = GenerateKeyStickelsCirculant(n, mm, mM, s, t);
            tic;

            [S, c] = Find_c_S_Circulant(minplusTropId(n), minplusTropId(n), U,s,t);
            
            
            cover = findMinimumCoverCirculant(S);

        
            [x, y] = SolveSystemGivenCover(cover, c);

            
            if ~isempty(x)
                K_attack = minplusMulti(ApplyGeneratorsMinPlus([0:d; x'], minplusTropId(n),s),minplusMulti(V,ApplyGeneratorsMinPlus([0:d; y'], minplusTropId(n),t)));

                
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
