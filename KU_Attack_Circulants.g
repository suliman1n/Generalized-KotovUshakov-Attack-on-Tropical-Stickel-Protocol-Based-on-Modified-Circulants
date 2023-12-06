
# Generates a random matrix A (interval matrices)
GenerateRandomMatrix := function(n, mm)
  return List([1..n], i -> List([1..n], function(j) 
if i=j then
  return 0;
else
  return  Random(mm, 2*mm );
  fi;
  end));
end;

GenerateRandDiagonal := function(n, mm,mM)
  return List([1..n], i -> List([1..n], function(j) 
if i=j then
  return Random(mm, mM );
else
  return infinity  ;
  fi;
  end));
end;


# Generates a random polynomial p in ZZ[x] of degree d, d in [1, D], p_i in [pm, pM].
GenerateRandomPolynomial := function(D, pm, pM)
  return List([1..D], i -> [i, Random(pm, pM)]);
end;


PowerOfMatrixMinPlus:=function(A,alpha,s)
return PowersMinPlus(A,alpha,s);
end;


PowersMinPlus := function(A, alpha, s)
  local n;
  n := Length(A);
  return List([1..n], i -> List([1..n], function(j) 
      if i >= j and (i-j) mod n = alpha then
            return 0;
        elif i < j and (i-j) mod n = alpha then
            return s;
        else
            return infinity;
        fi;
    end));
end;

# Returns a \otimes b.
ProductOfTwoScalarMinPlus := function(a, b)
  if IsInfinity(a) or IsInfinity(b) then
    return infinity;
  else 
    return a + b;
  fi;
end;


# Retruns A \otimes b.
ProductOfTwoMatricesMinPlus := function(A, B)
  local i, j, n, result; 
  n := Length(A);
  result := [];
  for i in [1..n] do
    Add(result, []);
    for j in [1..n] do
      Add(result[i], Minimum(List([1..n], k -> ProductOfTwoScalarMinPlus(A[i][k], B[k][j]))));
    od;
  od;
  return result;
end;


# Returns zero matrix of size nxn over min-plus algebra.
ZeroMatrixMinPlus := function(n)
  return List([1..n], i -> List([1..n], j -> infinity));
end;


# Returns ident matrix of size nxn over min-plus algebra.
IdentityMatrixMinPlus := function(n)
  return List([1..n], i -> List([1..n], function(j) 
      if i = j then 
        return 0; 
      else 
        return infinity; 
      fi;
    end));
end;




# Returns A + B.
SumOfTwoMatricesMaxPlus := function(A, B)
  return ListN(A, B, function(a, b) return ListN(a, b, function(x, y) return Minimum(x, y); end); end);
end;


# Returns s \otimes A.
ProductOfScalarAndMatrixMinPlus := function(s, A)
  return List(A, a -> List(a, x -> ProductOfTwoScalarMinPlus(s, x)));
end;


# Retruns p(A).
ApplyPolynomialMinPlus := function(p, A,s)
  if Length(p) = 0 then
    return ZeroMatrixMinPlus(Length(A));
  fi;
  return Iterated(List(p, c -> ProductOfScalarAndMatrixMinPlus(c[2], PowerOfMatrixMinPlus(A, c[1],s))), SumOfTwoMatricesMaxPlus);
end;


# Returns A - B. Some elements of the matrices can be infinity.
MinusMatrixFromMatrix := function(A, B)
  return ListN(A, B, function(a, b) return ListN(a, b, function (x, y)
      if y = infinity then
        return fail;
      elif x = infinity then
        return infinity;
      else  
        return x - y;
      fi;
  end); end);
end;



# Returns minimum of a matrix A and the set of corresponding indexes.
GetMinimumOfMatrix := function(A)
  local i, j, best;
  best := rec(val := A[1][1], inds := []);
  for i in [1..Length(A)] do
    for j in [1..Length(A[i])] do
      if A[i][j] < best.val then
        best := rec(val := A[i][j], inds := [[i, j]]);
      elif A[i][j] = best.val then
        AddSet(best.inds, [i, j]);
      fi;
    od;
  od;
  return best;
end;


# Retruns the simpex-table constructed by a set of covers and a fixed cover.
MakeSimplexMatrix := function(F, D, h)
  local A, i, j, S,k;
  k:=D;
  A := NullMat((k + 1)^2 + 1, 2 * (k + 1) + (k + 1)^2 + 1);
  for i in [0..k] do
    for j in [0..k] do
      A[i * (k + 1) + j + 1][i + 1] := 1;
      A[i * (k + 1) + j + 1][(k + 1) + j + 1] := 1;
    od;
  od;
  for S in F do
    A[S.ijminval[1].i * (k + 1) + S.ijminval[1].j + 1][2 * (k + 1) + (k + 1)^2 + 1] := -S.ijminval[1].val;
  od;
  for i in [1..(k + 1)^2] do
    A[i][2 * (k + 1) + i] := -1;
  od;
  for S in h do
    A[S[1] * (k + 1) + S[2] + 1][2 * (k + 1) + S[1] * (k + 1) + S[2] + 1] := 0;
  od;
  for i in [1..Length(A[1])] do  
    S := 0;
    for j in [1..Length(A) - 1] do
      S := S + A[j][i];
    od;
    A[(k + 1)^2 + 1][i] := -S;
  od;
  return A;
end;


# Retruns index of pivot column.
FindPivotColumn := function(M)
  local i, best, best_i, t;
  best := 0;
  best_i := fail;
  for i in [1..Length(M[1]) - 1] do
    t := M[Length(M)][i];
    if t < best then
      best := t;
      best_i := i;
    fi;
  od;
  return best_i;
end;


# Retruns index of pivot row.
FindPivotRow := function(M, l)
  local i, a, b, best_i, best;
  best_i := fail;
  best := infinity;
  for i in [1..Length(M) - 1] do  
    a := M[i][l];
    b := M[i][Length(M[i])];
    if a > 0 and b / a < best then
      best_i := i;
      best := b / a;
    fi;
  od;
  return best_i;
end;


# Performs recalculation of a simplex-table. Changes M, Bs and Ns.
Recalc := function(M, k, l, Bs, Ns)
  local t, i, j;
  t := Ns[l];
  Ns[l] := Bs[k];
  Bs[k] := t;

  for i in [1..Length(M)] do
    for j in [1..Length(M[1])] do 
      if i <> k and j <> l then
        M[i][j] := M[i][j] - M[i][l] * M[k][j] / M[k][l];
      fi;
    od;
  od;

  for i in [1..Length(M)] do
    if i <> k then
      M[i][l] := -M[i][l] / M[k][l];
    fi;
  od;

  for i in [1..Length(M[1])] do
    if i <> l then
      M[k][i] := M[k][i] / M[k][l];
    fi;
  od;

  M[k][l] := 1/M[k][l];
end;


# Returns the solution constructed by a simplex-table.
WriteDownSolution := function(M, Bs)
  local i, result;
  result := List([1..Length(Bs)], i -> 0);
  for i in [1..Length(Bs)] do
    result[Bs[i]] := M[i][Length(M[1])];
  od;
  return result;
end; 


# Applies the simplex-method to a matrix M and returns the corresponding solution.
# If there is no solution, then returns fail.
ApplySimplex := function(M)
  local k, l, Bs, Ns;
  Bs := [Length(M[1]) - 1 + 1..Length(M[1]) - 1 + Length(M) - 1];
  Ns := [1..Length(M[1]) - 1];
  while true do
    l := FindPivotColumn(M);
    if l = fail then
      break;
    fi;
    k := FindPivotRow(M, l);
    if k = fail then
      return fail;
    fi;
    Recalc(M, k, l, Bs, Ns);
  od;
  if M[Length(M)][Length(M[1])] <> 0 then
    return fail;
  fi;
  return WriteDownSolution(M, Bs);
end;


# Returns the number of i-th indexes. 
NumberOfIndexes := function(S, i)
  return Length(Set(List(S, c -> c[i])));
end;


# Compresses a set of pair [a1, L1], [a2, L2], ..., [an, Ln] to the set of pair, where all first components are unique.
Rar := function(G)
  local Find, H, i, S;

  Find := function(H, S)
    for i in [1..Length(H)] do
      if S.inds = H[i].inds then
        return i;
      fi;
    od;
    return fail;
  end;

  H := [];
  for S in G do
    i := Find(H, S);
    if i = fail then
      Add(H, S);
    else 
      Append(H[i].ijminval, S.ijminval);
    fi;
  od;
  return H;
end;


# Returns the compressed set of covers. This set contains all the minimal covers.
GetCompressedCovers := function(F)
  local M, N, P, Z;
  if Length(F) = 0 then
    return [[]];
  fi;
  Z := Rar(F);
  M := Filtered(Z, S -> Length(Difference(S.inds, Union(List(Filtered(Z, T -> S <> T), T -> T.inds)))) <> 0);
  N := Union(List(M, S -> S.inds));
  P := Rar(List(Filtered(List(Z, S -> rec(ijminval := S.ijminval, inds := Difference(S.inds, N))), 
      S -> Length(S.inds) <> 0)));
  if Length(P) > 0 then
    Sort(P, function(S, T) return Length(S.inds) > Length(T.inds); end);
    return List(Concatenation(List(GetCompressedCovers(
        Filtered(List(P, S -> rec(ijminval := S.ijminval, inds := Difference(S.inds, P[1].inds))), S -> Length(S.inds) <> 0)),
        S -> Concatenation([P[1]], S)),
        GetCompressedCovers(P{[2..Length(P)]})), S -> Concatenation(M, S));
  fi; 
  return [M];
end;


# Applies our attack.
ApplyAttack := function(A, B, u, D, pm,s,t)
  local Repack, F, G, H, S, T,k;
  k:=D;
  Repack := function(ij, min)
    return rec(ijminval := [rec(i := ij[1], j := ij[2], val := min.val)], inds := min.inds);
  end;

  F := Filtered(List(
      Cartesian([0..k], [0..k]), 
      ij -> Repack(ij, GetMinimumOfMatrix(MinusMatrixFromMatrix(
              ProductOfTwoMatricesMinPlus(PowerOfMatrixMinPlus(A, ij[1],s), PowerOfMatrixMinPlus(B, ij[2],t)), 
              ProductOfScalarAndMatrixMinPlus(-2 * pm, u))))), S -> S.ijminval[1].val <= 0);
  G := GetCompressedCovers(F);
  H := Union(List(G, S -> Cartesian(List(S, T -> List(T.ijminval, c -> [c.i, c.j])))));
  Sort(H, function(a, b) 
      return Length(a) < Length(b) or (Length(a) = Length(b) and 
          NumberOfIndexes(a, 1) * NumberOfIndexes(a, 2) < NumberOfIndexes(b, 1) * NumberOfIndexes(b, 2)); end);

  for S in H do
    T := ApplySimplex(MakeSimplexMatrix(F, D, S));
    if T <> fail then
      return [List([0..k], i -> [i, T[i + 1] + pm]), List([0..k], i -> [i, T[k + 1 + i + 1] + pm])];
    fi;
  od;
  return fail;
end;


RandomColumnVector := function(length, minValue, maxValue)
    local vector, i;
    vector := [];
    for i in [1..length] do
        Add(vector, Random([minValue..maxValue]));
    od;
    return vector;
end;




CreateModifiedCirculant := function(column, s)
    local n, circulant_matrix, i, j;
    n := Length(column);
    circulant_matrix := List([1..n], i -> List([1..n], j -> 0));

    for i in [1..n] do
        for j in [1..n] do
            circulant_matrix[i][j] := column[(((i - j) mod n + n) mod n) + 1];
        od;
    od;

    for i in [1..n] do
        for j in [i+1..n] do
            circulant_matrix[i][j] := circulant_matrix[i][j] + s;
        od;
    od;

    return circulant_matrix;
end;





TestAttack := function(n, D, pm, pM, ApplyAttack,s,t)
  local A1, A2, B1, B2 , u, v, KA, KB, KC, attack_result;
  A1 := CreateModifiedCirculant(RandomColumnVector(n,pm,pM),s);
  A2 := CreateModifiedCirculant(RandomColumnVector(n,pm,pM),t);

  B1 := CreateModifiedCirculant(RandomColumnVector(n,pm,pM),s);
  B2 := CreateModifiedCirculant(RandomColumnVector(n,pm,pM),t);

  u := ProductOfTwoMatricesMinPlus(A1, A2);
  v := ProductOfTwoMatricesMinPlus(B1, B2);

  KA := ProductOfTwoMatricesMinPlus(A1, ProductOfTwoMatricesMinPlus(v,A2));
  KB := ProductOfTwoMatricesMinPlus(B1, ProductOfTwoMatricesMinPlus(u,B2));
  if KA <> KB then
    return false;
  fi;

  attack_result := ApplyAttack(IdentityMatrixMinPlus(n),IdentityMatrixMinPlus(n), u, n-1, pm,s,t);
  if attack_result = fail then
    return false;
  fi;
  KC := ProductOfTwoMatricesMinPlus(ApplyPolynomialMinPlus(attack_result[1], A1,s), 
      ProductOfTwoMatricesMinPlus(v, ApplyPolynomialMinPlus(attack_result[2], B1,t)));
  if KA <> KC then
    return false;
  fi;
  return true;
end;




# Runs a set of tests.
TestSuite := function(n,D, pm, pM, ApplyAttack, numberOfTests,s,t)
  local st, et, i, ok, fl, maxCovers;
  st := Runtime();
  ok := 0;
  fl := 0;
  for i in [1..numberOfTests] do
    if TestAttack(n,D, pm, pM, ApplyAttack,s,t) then
      Print("OK\n");
      ok := ok + 1;
    else
      Print("FAIL\n");
      fl := fl + 1;
    fi;
  od;
  et := Runtime();
  Print(et - st, "\n");
  Print("OK: ", ok, "\n");
  Print("FAIL: ", fl, "\n");
end;