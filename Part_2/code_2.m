%% 1 фазовый 
clc;
figure;
f_x = @(x,y,z,a,b) a+b*x-x.*y./(1+x);
f_y = @(x,y,z,a,b) x.*y.*z./(1+x) - y;
f_z = @(x,y,z,a,b) 2 - z;
n = 5;
x = linspace(0,2,n);
y = linspace(0,3,n);
z = linspace(0,2,n);
a = 2;
b = 1;
[X,Y,Z] = meshgrid(x,y,z);
U = f_x(X,Y,Z,a,b);
V = f_y(X,Y,Z,a,b);
Q = f_z(X,Y,Z,a,b);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
newQ = V./(sqrt(U.^2 + V.^2));
quiver3(X,Y,Z,newU, newV, newQ)
hold on
plot3(a/b, 0, 2, 'r*')
plot3(1, 2*(a-b), 2, 'r*')
hold off
axis([0 2 0 3 0 2])
xlabel('x')
ylabel('y')
zlabel('z')
%% 2 фазовый
clc;
figure;
f_x = @(x,y,z,a,b) a+b*x-x.*y./(1+x);
f_y = @(x,y,z,a,b) x.*y.*z./(1+x) - y;
f_z = @(x,y,z,a,b) 2 - z;
n = 7;
x = linspace(-3,2,n);
y = linspace(-3,3,n);
z = linspace(-3,2,n);
a = 1;
b = 2;
[X,Y,Z] = meshgrid(x,y,z);
U = f_x(X,Y,Z,a,b);
V = f_y(X,Y,Z,a,b);
Q = f_z(X,Y,Z,a,b);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
newQ = V./(sqrt(U.^2 + V.^2));
quiver3(X,Y,Z,newU, newV, newQ)
hold on
plot3(a/b, 0, 2, 'r*')
plot3(1, 2*(a-b), 2, 'r*')
hold off
axis([-3 2 -3 3 -3 2])
xlabel('x')
ylabel('y')
zlabel('z')

%% 3 фазовый
clc;
figure;
f_x = @(x,y,z,a,b) a+b*x-x.*y./(1+x);
f_y = @(x,y,z,a,b) x.*y.*z./(1+x) - y;
f_z = @(x,y,z,a,b) 2 - z;
n = 5;
x = linspace(0,2,n);
y = linspace(0,3,n);
z = linspace(0,2,n);
a = 1;
b = 1;
[X,Y,Z] = meshgrid(x,y,z);
U = f_x(X,Y,Z,a,b);
V = f_y(X,Y,Z,a,b);
Q = f_z(X,Y,Z,a,b);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
newQ = V./(sqrt(U.^2 + V.^2));
quiver3(X,Y,Z,newU, newV, newQ)
hold on
plot3(a/b, 0, 2, 'r*')
plot3(1, 2*(a-b), 2, 'r*')
hold off
axis([0 2 0 3 0 2])
xlabel('x')
ylabel('y')
zlabel('z')

%% фазовый 4
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 30;
x = linspace(0,2,n);
y = linspace(0,5,n);
alp = 30;
gam = 0.7;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot([0, gam, 1], [0, alp*gam*(1-gam)/(1+gam), 0], 'r*')
axis([0 2 0 5])
xlabel('x')
ylabel('y')
%% фазовый 5
f_x = @(x,y,alp,gam) alp*x.^2.*(1-x)./(1+x) - x.*y;
f_y = @(x,y,alp,gam) -gam*y+x.*y;
n = 30;
x = linspace(0,2,n);
y = linspace(0,5,n);
alp = 20;
gam = 1.1;
[X,Y] = meshgrid(x,y);
U = f_x(X,Y,alp,gam);
V = f_y(X,Y,alp,gam);
newU = U./(sqrt(U.^2 + V.^2));
newV = V./(sqrt(U.^2 + V.^2));
quiver(X,Y,newU, newV, 0.5)
hold on
plot([0, 1], [0, 0], 'r*')
axis([0 2 0 5])
xlabel('x')
ylabel('y')

