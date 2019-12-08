% JAN RADZIMINSKI (293052)
% ==== ENUME - ASSIGNMENT B PROJECT ==== 
hold off;
close all;
clear;
clc;

% ========== TASK 1 =============
% Original function
Xarg = linspace(-1, 1, 1000); 
% this vector will be later used for making graphs in all cases
Y = fx(Xarg);

figure(1)
plot(Xarg, Y);
hold on;
title('Graph of original function');
legend('f(x)');

pause;
hold off;
close all;

%Generating points that will be used in aprox.
%N = 10
Xn1=linspace(-1, 1, 10);
Yn1=fx(Xn1);

%N = 20
Xn2=linspace(-1, 1, 20);
Yn2=fx(Xn2);

%N = 30
Xn3=linspace(-1, 1, 30);
Yn3=fx(Xn3);

% ========== TASK 2 =============

% Ploting aproximated funcitons

%N = 10
figure(1);
plot(Xarg , Y)
hold on;
plot(Xarg , AproxFx(10, 4, Xarg))
hold on;
plot(Xn1, Yn1, '*')
hold on;
title('Aproximations for N=10, K=4');
legend('Original f(x)', 'Aprox. for K=4', 'Points used for aproximation');

figure(2);
plot(Xarg , Y)
hold on;
plot(Xarg , AproxFx(10, 8, Xarg))
hold on;
plot(Xn1, Yn1, '*')
hold on;
title('Aproximations for N=10, K=8');
legend('Original f(x)', 'Aprox. for K=8', 'Points used for aproximation');

%N = 20
figure(4);
plot(Xarg , Y)
hold on;
plot(Xarg , AproxFx(20, 4, Xarg))
hold on;
plot(Xn2, Yn2, '*')
hold on;
title('Aproximations for N=20, K=4');
legend('Original f(x)', 'Aprox. for K=4', 'Points used for aproximation');

figure(5);
plot(Xarg , Y)
hold on;
plot(Xarg , AproxFx(20, 8, Xarg))
hold on;
plot(Xn2, Yn2, '*')
hold on;
title('Aproximations for N=20, K=8');
legend('Original f(x)', 'Aprox. for K=8', 'Points used for aproximation');

%N = 30

figure(7);
plot(Xarg , Y)
hold on;
plot(Xarg , AproxFx(30, 4, Xarg))
hold on;
plot(Xn3, Yn3, '*')
hold on;
title('Aproximations for N=30, K=4');
legend('Original f(x)', 'Aprox. for K=4', 'Points used for aproximation');

figure(8);
plot(Xarg , Y)
hold on;
plot(Xarg , AproxFx(30, 8, Xarg))
hold on;
plot(Xn3, Yn3, '*')
hold on;
title('Aproximations for N=30, K=8');
legend('Original f(x)', 'Aprox. for K=8', 'Points used for aproximation');

pause;

% ========== TASK 3 =============
hold off;
close all;
clc;

% Calculating RMS and Max error matrices
N = linspace(5, 50, 10);
for i=5:50
    for j=3:(i-1)
        S2tab(i, j) = S2(i, j, Xarg);
        SItab(i, j) = Sinf(i, j, Xarg);
    end
end

% Ploting 3D Graphs
figure(1)
surf(log10(S2tab));
title('3D graph of RMS errors for N={5,...,50}')
ylabel('N')
xlabel('K')
zlabel('RMS error value')

figure(2)
surf(log10(SItab));
title('3D graph of Max errors for N={5,...,50}')
ylabel('N')
xlabel('K')
zlabel('Max error value')

pause;
% ========== TASK 4 =============
hold off;
close all;
clc;

sig = linspace(-5, -1, 20);
k=1;
for p=sig
    RMSmin(k)=1000;
    Mmin(k)=1000;
for i=5:50
    [Xc1, Yc1] = Corrupted(i, p);
    for j=3:(i-1)
        S2ctab(i, j) = norm(AproxFx2(Yc1, j, Xarg) - fx(Xarg))/norm(fx(Xarg));
        SIctab(i, j) = norm(AproxFx2(Yc1, j, Xarg) - fx(Xarg), 'inf')/norm(fx(Xarg), 'inf');
        
%         finding smallest value:
        if(S2ctab(i, j)>0)
            if(S2ctab(i, j)<RMSmin(k)) RMSmin(k)=S2ctab(i, j);
            end
        end
        
        if(SIctab(i, j)>0)
            if(SIctab(i, j)<Mmin(k)) Mmin(k)=SIctab(i, j);
            end
        end
    end
    
end
k=k+1;
end

% Finding best fit to the points
sig=10.^sig;
space=logspace(-5, -1);
p = polyfit(sig, RMSmin, 3);
pm = polyval(p, space);

figure(1)
loglog(sig, RMSmin, '*');
hold on;
loglog(space, pm);
hold on;
xlabel('Standard deviation used (sigma)');
ylabel('Value of smallest RMS error found');

p = polyfit(sig, Mmin, 3);
pm = polyval(p, space);

figure(2)
loglog(sig, Mmin, '*');
hold on;
loglog(space, pm);
hold on;
xlabel('Standard deviation used (sigma)');
ylabel('Value of smallest Max error found');


% ================================
% Used functions:

function [Y] = AproxFx(N, K, Xarg)

Xn = linspace(-1, 1, N);
Yn(:) = fx(Xn(:));

Fi(:, :) = Gk(Xn, K);
A = ((Fi*Fi.')\(Fi*Yn.')).';

G(:, :) = Gk(Xarg, K);
Y=A*G;
end

function [Y] = AproxFx2(Yn, K, Xarg)

Xn = linspace(-1, 1, size(Yn,2));
Fi(:, :) = Gk(Xn, K);
A = ((Fi*Fi.')\(Fi*Yn.')).';
G(:, :) = Gk(Xarg, K);
Y=A*G;
end

function [Y] = fx(X)
Y(:)=-cos(pi*X(:)).*exp(-X(:)-1/3);
end

function [Y] = Gk(N, K)
xk=linspace(-1, 1, K);
Y(:, :) = 5/(2*pi).*exp(-25*(N(1, :)-xk(:)).^2);
end

function [Y] = S2(N, K, Xarg)
Y=norm(AproxFx(N, K, Xarg) - fx(Xarg))/norm(fx(Xarg));
end

function [Y] = Sinf(N, K, Xarg)
Y=norm(AproxFx(N, K, Xarg) - fx(Xarg), 'inf')/norm(fx(Xarg), 'inf');
end

function [X, Y] = Corrupted(N, sig_pow)
X=linspace(-1, 1, N);
Y = fx(X);
R=randn(1, N)*(10^sig_pow);
Y=Y+R;
end

% Made by Jan Radziminski
