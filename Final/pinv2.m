function [H] = pinv2(A,tolerance)
%PINV2 Calculates the pseudoinverse of a Matrix
%   With tolerance levels for the singular value
[U, S, V] = svd(A);
Sp = diag(S);
Sp_n = Sp/max(Sp);

n = max(find(Sp_n > tolerance));

U1 = U(:,1:n);
S1 = diag(Sp(1:n));
V1 = V(:,1:n);
H = V1*S1^-1*U1';


end

