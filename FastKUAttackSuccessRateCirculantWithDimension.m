function FastKUAttackSuccessRateCirculantWithDimension(max_n, mm, mM, trials)
    successRates = zeros(1, max_n);

    for n = 1:max_n
        successRates(n) = FastKUAttackSuccessRateCirculant(n,mm, mM,trials);
    end

    
    figure;
    plot(1:max_n, successRates, '-o');
    title('Success Rate vs Dimension');
    xlabel('Dimension (n)');
    ylabel('Success Rate');
    grid on;
end