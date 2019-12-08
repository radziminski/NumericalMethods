% JAN RADZIMINSKI (293052)
% ==== ENUME - ASSIGNMENT C PROJECT ==== 
hold off;
close all;
clear;
clc;

% ==== Task 1 ====
% Initial Data
y0 = [0 2];

% Note: Since I wanted to make my program more general, both ode113 and my
% Lobatto and Euler functions take A as argument, so they should work
% also for similar equations with different coefficients - one should just
% change A matrix
A = [0, 1; -10/9, -6/9];

h=0.01;
tspan = 0:h:10;

% Matlab ode113 Soluttion
options = odeset('AbsTol', 1*10^(-16), 'RelTol', 2.22045*10^-14);
[tm,ym] = ode113(@(tm,ym) odefcn(ym, A), tspan, y0, options);
figure(1)
plot(tm,ym(:,1),'-') %t,y(:,2),'.-')
hold on

% My Solution with Lobatto IID method
[X, Y] = Lobatto(tspan, y0, A, h);
figure(1)
plot(X, Y(1, :))
hold on
title('Solution to given ODE')
xlabel('T')
ylabel('Y')
legend('y(t) obtained with ode113', 'y(t) obtained with LobattoIIID')

pause;

% ==== Task 2 & 3 ====
clc;
hold off;
close all;

% Not: With first line of hi activated program works quite fast when with second
% line one should wait quite a while till calculations are finished
hi = logspace(-1, -3, 30);
% hi = logspace(-1, -5, 500);
k=1;

for i = hi
    tspan = 0:i:10;
    RMS_L(k) = S2(tspan, y0, A, i, options, 'L');
    MAX_L(k) = Sinf(tspan, y0, A, i, options, 'L');
    RMS_E(k) = S2(tspan, y0, A, i, options, 'E');
    MAX_E(k) = Sinf(tspan, y0, A, i, options, 'E');
    k=k+1;
end

% Plotting RMS Errors
figure(1)
loglog(hi, RMS_L, '-');
hold on;
loglog(hi, RMS_E, '-');
hold on;
title('Root Mean Squared Error')
xlabel('Integration step h')
ylabel('Error value')
legend('RMS for Lobatto IIID method', 'RMS for explicit Euler method')

% Plotting Maximum Errors
figure(2)
loglog(hi, MAX_L, '-');
hold on;
loglog(hi, MAX_E, '-');
hold on;
title('Maximum Error')
xlabel('Integration step h')
ylabel('Error value')
legend('MAX for Lobatto IIID method', 'MAX for explicit Euler method');




% ==== FUNCTIONS USED: ====
% Lobatto IIID implicit Runge Kutta method:
function [t, y] = Lobatto(tspan, y0, A, h)
y(1, 1) =  y0(1);
y(2, 1) = y0(2);
I = eye(size(A));
M = [I-((1/2)*h*A), -(1/2)*h*A; (1/2)*h*A, I-((1/2)*h*A)];
t=tspan;
for i=2:size(t, 2)
    R=[A*y(:, i-1); A*y(:, i-1)];
    f = (M\R);
    y(1, i) = y(1, i-1) + h/2*( f(1) + f(3) );
    y(2, i) = y(2, i-1) + h/2*( f(2) + f(4) ); 
end
end

% Euler explicit Runge Kutta method:
function [t, y] = Euler(tspan, y0, A, h)
y(1, 1) =  y0(1);
y(2, 1) = y0(2);
t=tspan;

for i=2:size(t, 2)
    y(:, i) = y(:, i-1) + h*A*y(:, i-1);
end

end

% RMS Error (opt are ode113 options and r is which error (for Lobatto 'L' for Euler sth else):
function [Y] = S2(t, y0, A, h, opt, r)

if r == 'L'
[t, y1] = Lobatto(t, y0, A, h);
else
[t, y1] = Euler(t, y0, A, h);
end

[~,ym] = ode113(@(tm,ym) odefcn(ym, A), t, y0, opt);
ym=ym.';
Y=norm((y1(1, :) - ym(1, :)))/norm(ym(1, :));
end

% Max Error (opt are ode113 options and r is which error (for Lobatto 'L' for Euler sth else):
function [Y] = Sinf(t, y0, A, h, opt, r)

if r == 'L'
[t, y1] = Lobatto(t, y0, A, h);
else
[t, y1] = Euler(t, y0, A, h);
end

[~,ym] = ode113(@(tm,ym) odefcn(ym, A), t, y0, opt);
ym=ym.';
Y=norm(y1(1, :) - ym(1, :), 'inf')/norm(ym(1, :), 'inf');
end

% Function used in ode113 method 
function dydt = odefcn(y, A)
dydt = zeros(2,1);
dydt(1) = A(1, 1)*y(1) + A(1, 2)*y(2);
dydt(2) = A(2, 1)*y(1)+ A(2, 2)*y(2);
end