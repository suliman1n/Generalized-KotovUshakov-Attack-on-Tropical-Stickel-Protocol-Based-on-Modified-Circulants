function [key,U,V,A,B] = GenerateKeyStickels(n, mm, mM, D, pm, pM)
    
    A = randi([mm, mM], n);
    B = randi([mm, mM], n);

    
    p1 = GenerateRandomPolynomial(D, pm, pM);
    p2 = GenerateRandomPolynomial(D, pm, pM);
    q1 = GenerateRandomPolynomial(D, pm, pM);
    q2 = GenerateRandomPolynomial(D, pm, pM);

    
    U = minplusMulti(ApplyPolynomialMinPlus(p1, A), ApplyPolynomialMinPlus(p2, B));
    V = minplusMulti(ApplyPolynomialMinPlus(q1, A), ApplyPolynomialMinPlus(q2, B));

    
    KA = minplusMulti(ApplyPolynomialMinPlus(p1, A), minplusMulti(V, ApplyPolynomialMinPlus(p2, B)));
    KB = minplusMulti(ApplyPolynomialMinPlus(q1, A), minplusMulti(U, ApplyPolynomialMinPlus(q2, B)));

    
    if isequal(KA, KB)
        key = KA;
    else
        key = [];
    end
end
