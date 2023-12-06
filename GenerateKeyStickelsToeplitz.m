function [key,U,V] = GenerateKeyStickelsToeplitz(n, mm, mM)
    A1=CreateUpperToeplitz(randi([mm,mM],n,1));
    A2=CreateLowerToeplitz(randi([mm,mM],n,1));

    B1=CreateUpperToeplitz(randi([mm,mM],n,1));
    B2=CreateLowerToeplitz(randi([mm,mM],n,1));

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
