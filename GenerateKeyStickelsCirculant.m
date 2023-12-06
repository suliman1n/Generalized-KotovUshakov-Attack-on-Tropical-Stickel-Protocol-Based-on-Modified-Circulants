function [key,U,V] = GenerateKeyStickelsCirculant(n, mm, mM, s, t)
    A1 = CreateModifiedSCirculant(RandomColumnVector(n, mm, mM), s);
    A2 = CreateModifiedSCirculant(RandomColumnVector(n, mm, mM), t);

    B1 = CreateModifiedSCirculant(RandomColumnVector(n, mm, mM), s);
    B2 = CreateModifiedSCirculant(RandomColumnVector(n, mm, mM), t);

    U = minplusMulti(A1, A2);
    V = minplusMulti(B1, B2);

    KA = minplusMulti(A1, minplusMulti(V, A2));
    KB = minplusMulti(B1, minplusMulti(U, B2));

   if isequal(KA, KB)
        key = KA;
    else
        key = [];
    end
end
