clear all
close all
clc

% Drop-Wave
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2));

% Griewank
% f = @(x,y) ((x.^2/4000)+(y.^2/4000))-(cos(x).*cos(y/sqrt(2)))+1;

% Rastrigin
% f = @(x,y) 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y);

f = @(x,y) (x-2).^2 + (y-2).^2; %Esfera

xl = [-5 -5]';
xu = [5 5]';

G = 150;
N = 50;
D = 2;

c1 = 2.05; % Cognitivo
c2 = 2.05; % Social
phi = c1 + c2;
K = 2 / (abs(2-phi-sqrt(phi^2-4*phi)));

x = zeros(D,N);
v = zeros(D,N); % velocidades
xb = zeros(D,N); % xb mejores posiciones
fitness = zeros(1,N);
f_plot = zeros(1,G);

for i=1:N
    x(:,i) = xl + (xu-xl).*rand(D,1);
    v(:,i) = 0.5*randn(D,1);
    xb(:,i) = x(:,i);

    fitness(i) = f(x(1,i),x(2,i));
end

for g=1:G
    Plot_Contour(f,x,xl,xu);

    for i=1:N
        fx = f(x(1,i),x(2,i));

        if fx < fitness(i)
            xb(:,i) = x(:,i);
            fitness(i) = fx;
        end
    end

    [fx_best, I_best] = min(fitness);
    f_plot(g) = fx_best;

    for i=1:N
        v(:,i) = K*v(:,i) + rand()*c1*(xb(:,i)-x(:,i)) + rand()*c2*(xb(:,I_best)-x(:,i));
        x(:,i) = x(:,i) + v(:,i);
    end
end

%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)
title('Gráfica de convergencia','FontSize',15)

% figure
% Plot_Contour(f, xb, xl, xu)

figure
Plot_Surf(f, xb, xl, xu)
disp(['minimo global en: x=' num2str(xb(1,I_best)) ', y=' num2str(xb(2,I_best)) ', f(x,y)=' num2str(f(xb(1,I_best),xb(2,I_best)))])