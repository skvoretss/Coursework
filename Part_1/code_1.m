%% бифуркационная диаграмма
clc;

f = @(u,r) r*u.^(3/2).*(1-u);
n_r = 400;
n_dots = 100;
res = zeros(n_r * n_dots, 2);
r = 3.5;
for j = 1:n_r
    r = r + 0.005;
    u0 = 0.1;
    for i = 1:500
        u0 = f(u0, r);
    end
    for i = 1:n_dots
        u0 = f(u0, r);
        res(j * i, 1) = r;
        res(j * i, 2) = u0;
    end
end
hold on
plot(res(:,1), res(:,2), 'g.')
xlabel('r')
ylabel('u')
axis([3.5 5.5 0 1.5])
res = zeros(n_r * n_dots, 2);
r = 2.5;
for j = 1:n_r
    r = r + 0.005;
    u0 = 3;
    for i = 1:500
        u0 = f(u0, r);
    end
    for i = 1:n_dots
        u0 = f(u0, r);
        res(j * i, 1) = r;
        res(j * i, 2) = u0;
    end
end
hold on
plot(res(:,1), res(:,2), 'g.')
%% Показатель Ляпунова
r = 3.6;
res = zeros(300, 1);
for i = 1:300
    r = r + 0.005;
    f = @(u) r*u.^(3/2).*(1-u);
    f_u = @(u) r/2*sqrt(u).*(3-5*u);
    u_t = 0.1;
    n = 1000;
    sum = log(abs(f_u(u_t)))/n;
    for j = 1:n
        u_t = f(u_t);
        sum = sum + log(abs(f_u(u_t)))/n;
    end
    res(i) = sum;
end
plot(linspace(0.005, 2, 300), res, 'g', [0, 2], [0 0], 'k')
xlabel('r')
ylabel('h')
%% устойчивость (1,1) при r < 1/48
f = @(u,v,r) sqrt(u).*exp(r*(1-v.^3));
n = 50;
t = 1:n;
res = zeros(n, 1);
res(1) = 1.1;
res(2) = res(1);
r = 1/96;
for i = 3:n
    res(i) = f(res(i-1),res(i-2),r);
end
plot(t,res, 'r*', t,res, 'r')
xlabel('t')
ylabel('u_t')
el('u_t')

%% Неймарк - Сакер №1
clc;
b = 1;
r = 2.82;
f = @(u,v,r) r*u.^(3/2).*(1-b*v);
u0 = 0.54;
v0 = u0;
n = 5000;
res = zeros(n,2);
for i = 1:n
    t = f(u0, v0, r);
    v0 = u0;
    u0 = t;
    res(i,1) = u0;
    res(i,2) = v0;
end
plot(res(:,1), res(:,2), 'g.')
xlabel('u')
ylabel('v')

%% Неймарк - Сакер №2
clc;
b = 1;
r = 2.84;
f = @(u,v,r) r*u.^(3/2).*(1-b*v);
u0 = 0.54;
v0 = u0;
n = 5000;
res = zeros(n,2);
for i = 1:n
    t = f(u0, v0, r);
    v0 = u0;
    u0 = t;
    res(i,1) = u0;
    res(i,2) = v0;
end
plot(res(:,1), res(:,2), 'g.')
xlabel('u')
ylabel('v')