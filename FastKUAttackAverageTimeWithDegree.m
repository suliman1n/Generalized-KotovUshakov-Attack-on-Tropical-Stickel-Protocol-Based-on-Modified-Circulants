function FastKUAttackAverageTimeWithDegree(n, mm, mM, pm, pM, maxDegree, trials)
    averageTimes = zeros(1, maxDegree);

    for d = 1:maxDegree
        totalElapsedTime = 0;

        for trial = 1:trials
            [key, U, V, A, B] = GenerateKeyStickels(n, mm, mM, d, pm, pM);
            tic;
            [S, c] = Find_c_S(A, B, U, d);
       
            cover = findMinimumCover(S,n,d);
            
            [x, y] = SolveSystemGivenCover(cover, c);

            
            if ~isempty(x)
                K_attack = minplusMulti(ApplyPolynomialMinPlus([0:d; x'], A), minplusMulti(V, ApplyPolynomialMinPlus([0:d; y'], B)));

                
                if isequal(K_attack, key)
                    totalElapsedTime = totalElapsedTime + toc;
                end
            end
        end

        
        averageTimes(d) = totalElapsedTime / trials;
    end

    
    figure;
    plot(1:maxDegree, averageTimes, '-o');
    title('Average Time vs Degree');
    xlabel('Degree (d)');
    ylabel('Average Time (seconds)');
    grid on;
end
