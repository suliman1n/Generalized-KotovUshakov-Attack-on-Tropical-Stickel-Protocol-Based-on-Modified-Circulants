function FastKUAttackSuccessRateWithDegree(n, mm, mM, pm, pM, maxDegree, trials)
    successRates = zeros(1, maxDegree);

    for d = 1:maxDegree
        successRates(d) = FastKUAttackSuccessRate(n, d, mm, mM, pm, pM, trials);
    end

    
    figure;
    plot(1:maxDegree, successRates, '-o');
    title('Success Rate vs Degree');
    xlabel('Degree (d)');
    ylabel('Success Rate');
    grid on;
end