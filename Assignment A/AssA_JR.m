% General Notes
% This script prints out answeres to all tasks in our Assignemnt. After
% each task it stops by pause command and wait for any key to be pressed to
% move to the next task. Also after every tasks it usess commands clear,
% clc hold off and close all.

% TASK 1
format rational; % changing format to rational for easier comparison of matrices with matrix shown in assignment
disp('===== TASK 1 =====');
fprintf('\nTesting Function that generates matricies as given in assignment, \nfor different values of N and x:')
fprintf('\nFor N=1, x=1 generating matrix A:\n')
A = Generate(1, 1)
fprintf('\nFor N=2, x=0 generating matrix B:\n')
B = Generate(2, 0)
fprintf('\nFor N=3, x=-2 generating matrix C:\n')
C = Generate(3, -2)
format short;
fprintf('\nFor N=4, x=sqrt(2) generating matrix D: (shown in short format since sqrt(2) is irrational)\n')
D = Generate(4, sqrt(2))
format rational;
fprintf('\nFor N=5, x=1, generating matrix E:\n')
E = Generate(5, 1)

pause;
clear;
clc;
close all;
hold off;
format short;
% =========================================================================

%TASK 2
disp('=== TASK 2 ===');
disp('');
% Precisions to be changed by user:
precision = 0.00001;
aprox = 0.00001;

fprintf('Searching for smalles alfa with condition that det < %d; and searching with steps = %d;\n', aprox, precision)
fprintf('For N=3:')
alfa3=SmallestAlfa(3, precision, aprox)

fprintf('For N=10:')
alfa10=SmallestAlfa(3, precision, aprox)

fprintf('For N=20:')
alfa20=SmallestAlfa(3, precision, aprox)

% Graphs:
low = 1 - 0.01;
hi = 1 + 0.01;
i=low;
k=1;
while (i<hi)
    X(1, k)=i;
    YD3(1, k) = det(Generate(3, log(i)));
    YC3(1, k) = cond(Generate(3, log(i)));
    
    YD10(1, k) = det(Generate(10, log(i)));
    YC10(1, k) = cond(Generate(10, log(i)));
    
    YD20(1, k) = det(Generate(20, log(i)));
    YC20(1, k) = cond(Generate(20, log(i)));
    
    i = i + 0.0001;
    k=k+1;
end

figure(1)
semilogy(X, YD3);
hold on
semilogy(X, YD10);
hold on
semilogy(X, YD20);
title('Determinant on alfa for N={3, 10, 20} and x=ln(alfa)')
xlabel('alfa')
ylabel('determinant')
legend ('N=3', 'N=10', 'N=20')
hold on

figure(2)
semilogy(X, YC3);
hold on
semilogy(X, YC10);
hold on
semilogy(X, YC20);
title('Cond on alfa for N={3, 10, 20} and x=ln(alfa')
xlabel('alfa')
ylabel('cond')
legend ('N=3', 'N=10', 'N=20')
hold on

pause;
hold off;
close all;
clear;
clc;
% =========================================================================

% TASK 3
A = [2 1; 1 3]
B = [2 -1 0; -1 2 -1; 0 -1 2]
C = [9 1 2 3; 1 8 0 4; 2 0 7 5; 3 4 5  6]
disp(' ');
disp('Inverting matricies using LU factorisation:');
Ainv = InverseOfMatrixLU(A)
Binv = InverseOfMatrixLU(B)
Cinv = InverseOfMatrixLU(C)
disp(' ');
disp('Checking corectness of inverse of matrix using LU factorisation (by compering multiplication of original')
disp('matrix and its inverse with identity matirx):');
disp('A 2x2 matrix: ');
round(A*Ainv, 5)==eye(2)
disp('B 3x3 matrix: ');
round(B*Binv, 5)==eye(3)
disp('C 4x4 matrix: ');
round(C*Cinv, 5)==eye(4)
disp(' ');

disp('Inverting matricies using LLT factorisation:');
A2inv = InverseOfMatrixLLT(A)
B2inv = InverseOfMatrixLLT(B)
C2inv = InverseOfMatrixLLT(C)
disp(' ');
disp('Checking corectness of inverse of matrix using LLT factorisation (by compering multiplication of original')
disp('matrix and its inverse with identity matirx):');
disp('A 2x2 matrix: ');
round(A*A2inv, 5)==eye(2)
disp('B 3x3 matrix: ');
round(B*B2inv, 5)==eye(3)
disp('C 4x4 matrix: ');
round(C*C2inv, 5)==eye(4)

pause;
clear;
clc;
% =========================================================================
% Task 4 and 5

% Legend:
% N - size of matrix, used to generate it
% AN - matrix generated for size N (A3, A10, A20)
% AN(row, column, k that was used to generate matrix)
% 
% AinvN - inverse of matrix AN
% AinvN(row, column, k that was used to generate original matrix, factorization that was used: LU=1, LLT=2)
% 
% RMSN - Root mean squared error for matricies of size N
% RMSN(k that was used to genearate original matrix AN, factorization that was used: LU=1, LLT=2)
% 
% ME3 = Maximum error of matrivies of size N - arguments same as in RMS

% Calculations:
for k=1:22
    x(1, k) = 2^(k-1)/300;
    A3(:, :, k) = Generate(3, x(1, k));
    Ainv3(:,  :, k, 1) = InverseOfMatrixLU(A3(:, :, k));
    Ainv3(:, :, k, 2) = InverseOfMatrixLLT(A3(:, :, k));
    
    RMS3(k, 1) = norm2((A3(:, :, k)*Ainv3(:, :, k, 1)-eye(3)));
    RMS3(k, 2) = norm2((A3(:, :, k)*Ainv3(:, :, k, 2)-eye(3)));
    ME3(k, 1) = normInf((A3(:, :, k)*Ainv3(:, :, k, 1)-eye(3)));
    ME3(k, 2) = normInf((A3(:, :, k)*Ainv3(:, :, k, 2)-eye(3)));

%N = 10:
    A10(:, :, k) = Generate(10, x(1, k));
    Ainv10(:,  :, k, 1) = InverseOfMatrixLU(A10(:, :, k));
    Ainv10(:, :, k, 2) = InverseOfMatrixLLT(A10(:, :, k));
    
    RMS10(k, 1) = norm2((A10(:, :, k)*Ainv10(:, :, k, 1)-eye(10)));
    RMS10(k, 2) = norm2((A10(:, :, k)*Ainv10(:, :, k, 2)-eye(10)));
    ME10(k, 1) = normInf((A10(:, :, k)*Ainv10(:, :, k, 1)-eye(10)));
    ME10(k, 2) = normInf((A10(:, :, k)*Ainv10(:, :, k, 2)-eye(10)));

%N = 20:

    A20(:, :, k) = Generate(20, x(1, k));
    Ainv20(:,  :, k, 1) = InverseOfMatrixLU(A20(:, :, k));
    Ainv20(:, :, k, 2) = InverseOfMatrixLLT(A20(:, :, k));
    
    RMS20(k, 1) = norm2((A20(:, :, k)*Ainv20(:, :, k, 1)-eye(20)));
    RMS20(k, 2) = norm2((A20(:, :, k)*Ainv20(:, :, k, 2)-eye(20)));
    ME20(k, 1) = normInf((A20(:, :, k)*Ainv20(:, :, k, 1)-eye(20)));
    ME20(k, 2) = normInf((A20(:, :, k)*Ainv20(:, :, k, 2)-eye(20)));
    
    %for norm comparison
    NormDiff2(k)=abs(norm2(A20(:, :, k))-norm(A20(:, :, k), 2));
    NormDiffInf(k)=abs(normInf(A20(:, :, k))-norm(A20(:, :, k), 'inf'));
end

figure(1)
semilogy(x(1, :), RMS3(:, 1));
hold on
semilogy(x(1, :), RMS10(:, 1));
hold on
semilogy(x(1, :), RMS20(:, 1));
hold on
semilogy(x(1, :), RMS3(:, 2));
hold on
semilogy(x(1, :), RMS10(:, 2));
hold on
semilogy(x(1, :), RMS20(:, 2));
hold on
title('Root Mean Square Error for inverses made with LU and LLT; x=2^k/300 (k={0, 1, ..., 21})')
xlabel('x')
ylabel('error value')
legend('LU, N=3', 'LU, N=10', 'LU, N=20', 'LLT, N=3', 'LLT, N=10', 'LLT, N=20')


figure(2)
semilogy(x(1, :), ME3(:, 1));
hold on
semilogy(x(1, :), ME10(:, 1));
hold on
semilogy(x(1, :), ME20(:, 1));
hold on
semilogy(x(1, :), ME3(:, 2));
hold on
semilogy(x(1, :), ME10(:, 2));
hold on
semilogy(x(1, :), ME20(:, 2));
hold on
title('Maximum Error for inverses made with LU and LLT; x=2^k/300 (k={0, 1, ..., 21})')
xlabel('x')
ylabel('error value')
legend('LU, N=3', 'LU, N=10', 'LU, N=20', 'LLT, N=3', 'LLT, N=10', 'LLT, N=20')

% Figures with LU and LLT factorisations drwed independantly on different
% graphs:
% figure(3)
% semilogy(x(1, :), RMS3(:, 2));
% hold on
% semilogy(x(1, :), RMS10(:, 2));
% hold on
% semilogy(x(1, :), RMS20(:, 2));
% hold on
% title('Root Mean Square Error for inverses made with LLT; x=2^k/300 (k={0, 1, ..., 21})')
% xlabel('x')
% ylabel('error value')
% legend('N=3', 'N=10', 'N=20')

% figure(4)
% semilogy(x(1, :), ME3(:, 2));
% hold on
% semilogy(x(1, :), ME10(:, 2));
% hold on
% semilogy(x(1, :), ME20(:, 2));
% hold on
% title('Maximum Error for inverses made with LLT; x=2^k/300 (k={0, 1, ..., 21})')
% xlabel('x')
% ylabel('error value')
% legend('N=3', 'N=10', 'N=20')

%Norms comparison
disp('Comparing bult in norm2 functions with these implemented by me (on all matrices with N=20 from Task 4):')
disp('Displaying biggest difference between them:')
max(NormDiff2(:))
disp('The same as previous but using norm inf:')
max(NormDiffInf(:))
%Matrix comparison
disp('Comparing matricies from task4 inverted by LU or LLT functions by my functions and build in inv function')
disp('by printing max(Ainv - inv(A)) which is maximal difference of single value in these matricies :')
disp(' ')
% Compering LU
disp('Compering inverses calculated using LU method: ')
disp('N=3:')
disp('k=0:')
max(max(Ainv3(:, :, 1, 1)-inv(A3(:, :, 1))))
disp('k=10:')
max(max(Ainv3(:, :, 11, 1)-inv(A3(:, :, 11))))
disp('k=20:')
max(max(Ainv3(:, :, 21, 1)-inv(A3(:, :, 21))))
disp(' ')
disp('N=10:')
disp('k=0:')
max(max(Ainv10(:, :, 1, 1)-inv(A10(:, :, 1))))
disp('k=10:')
max(max(Ainv10(:, :, 11, 1)-inv(A10(:, :, 11))))
disp('k=20:')
max(max(Ainv10(:, :, 21, 1)-inv(A10(:, :, 21))))
disp(' ')
disp('N=20:')
disp('k=0:')
max(max(Ainv20(:, :, 1, 1)-inv(A20(:, :, 1))))
disp('k=10:')
max(max(Ainv20(:, :, 11, 1)-inv(A20(:, :, 11))))
disp('k=20:')
max(max(Ainv20(:, :, 21, 1)-inv(A20(:, :, 21))))
disp(' ')

% Compering LLT
disp('Compering inverses calculated using LLT method: ')
disp('N=3:')
disp('k=0:')
max(max(Ainv3(:, :, 1, 2)-inv(A3(:, :, 1))))
disp('k=10:')
max(max(Ainv3(:, :, 11, 2)-inv(A3(:, :, 11))))
disp('k=20:')
max(max(Ainv3(:, :, 21, 2)-inv(A3(:, :, 21))))
disp(' ')
disp('N=10:')
disp('k=0:')
max(max(Ainv10(:, :, 1, 2)-inv(A10(:, :, 1))))
disp('k=10:')
max(max(Ainv10(:, :, 11, 2)-inv(A10(:, :, 11))))
disp('k=20:')
max(max(Ainv10(:, :, 21, 2)-inv(A10(:, :, 21))))
disp(' ')
disp('N=20:')
disp('k=0:')
max(max(Ainv20(:, :, 1, 2)-inv(A20(:, :, 1))))
disp('k=10:')
max(max(Ainv20(:, :, 11, 2)-inv(A20(:, :, 11))))
disp('k=20:')
max(max(Ainv20(:, :, 21, 2)-inv(A20(:, :, 21))))


pause;
% =========================================================================
% ------------------------------------------------------------------------
% Functions used:

% Function generating the matrix
function [A] = Generate(N,x) 
A = zeros(N);
A(1,1)=x^2;
for i=2:N
    A(1,i)=2*x/(3*((-1)^i));
    A(i,1)=2*x/(3*((-1)^i));
end
for i=2:N
    for j=2:N
        A(i,j)=j*4/(9*((-1)^(j-i)));
        A(j, i)=A(i, j);
    end
end
end

% Function finding smallest alpha with given precisions
function [alfa] = SmallestAlfa(N, precision, aprox) 
i=precision;
while(i<2)
        A=Generate(N, log(i));
        if (det(A)<aprox) break;
        end
        i=i+precision;
end
alfa=i;
end

% Function inverting matrices using LU factorisation 
function [B] = InverseOfMatrixLU(A)
N = size(A, 1);
Z = eye(N);
% [L, U] = FactorizationLU(A);  <- With this line active, function uses LU
% factorisation function written by me
[L, U, P] = lu(A);

for m=1:N
    
    %forward substitution - L*Y=Z  => Y = Z/L
    Y(1, m) = Z(1,m)/L(1, 1);
for i=2:N
    Y(i, m) = (Z(i, m)-(L(i, (1:(i-1)))* Y((1:(i-1)), m)))/L(i, i);
end

    %backward substitution - U*X=Y  => X = Y/U
    X(N, m) = (Y(N, m))/U(N, N);
for i=2:N
    X(N-i+1, m) = (Y(N-i+1, m)-(U(N-i+1, ((N-i+2):N))*X(((N-i+2):N), m)))/U(N-i+1, N-i+1);
end

end
%X - inverse of A with pivoted columns (becouse of built in LU function
% uses pivoting)
B=X*P;

end     

% Function inverting matrices using Cholesky-Banachiewicz factorisation 
function [B] = InverseOfMatrixLLT(A)

N = size(A, 1);
Z = eye(N);
LT = chol(A);
L = LT.';
for m=1:N
    
    %forward substitution - L*Y=Z  => Y = Z/L
    Y(1, m) = Z(1,m)/L(1, 1);
for i=2:N
    Y(i, m) = (Z(i, m)-(L(i, (1:(i-1)))* Y((1:(i-1)), m)))/L(i, i);
end

    %backward substitution - LT*X=Y  => X = Y/LT
    X(N, m) = (Y(N, m))/LT(N, N);
for i=2:N
    X(N-i+1, m) = (Y(N-i+1, m)-(LT(N-i+1, ((N-i+2):N))*X(((N-i+2):N), m)))/LT(N-i+1, N-i+1);
end

end
B=X;

end

% Norm functions
function [B] = norm2(A)
X = A*A.';
B = max(sqrt(eig(X)));
end
function [B] = normInf(A)
B = max((sum(abs(A), 2)));
end

%   Functions made before i knew that we can use built in function, finaly not used
%   in the project, but I decided to attach them anyway:
function [L, LT] = FactorizationLLT(A)
N=size(A,1);
L = zeros(N);
LT = zeros(N);
for n=1:N
   x = 0;
   for i=1:(n-1)
       x = x + L(n, i)*L(n, i);
   end
   if((A(n, n) - x)<0) 
       disp('Given Matrix is not symetric definite positive matrix, returning matrix filled with 0.')
       L = zeros(N);
       LT = zeros(N);
       return;
   end
   L(n, n) = sqrt(A(n, n) - x);
   for v=(n+1):N
   x = 0;
   for i=1:(n-1)
       x = x + L(v, i)*L(n, i);
   end

   L(v, n) = (A(v, n) - x)/L(n, n);
   end
end
LT = L.';
end
function [L, U] = FactorizationLU(A)
N=size(A,1);
B = ones(N,1);
L = eye(N);
U = zeros(N);
for step=1:(N-1)
    A(:, :, step+1)=A(:, :, step);
   for n=(2+step-1):N
       L(n, step) = A(n, step, step)/A(step, step, step);
       for v=(1+step-1):N
           A(n, v, step+1) = A(n, v, step) - L(n, step)*A(step, v, step);
           B(n, 1, step+1) = B(n, 1, step) - L(n, step)*B(step, 1, step);
       end
   end
end
U=A(:, :, N);

end
