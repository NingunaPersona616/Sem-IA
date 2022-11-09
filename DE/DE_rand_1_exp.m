clear all
close all
clc

% Drop-Wave
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2));

% Rastrigin
% f = @(x,y) 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y);

%Ackley
f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y))) + 20+exp(1); 

% f = @(x,y) (x-2).^2 + (y-2).^2;

xl = [-5 -5]';
xu = [5 5]';

D = 2;
G = 150;
N = 50;
F = 0.6;
CR = 0.9;

x = zeros(D, N);
fitness = zeros(1,N);
f_plot = zeros(1, G);

for i=1:N
    x(:, i) = xl+(xu-xl).*rand(D, 1);
    fitness(i) = f(x(1, i), x(2, i));
end  

for n=1:G
    for i=1:N   %Por cada generacion se actualizan las posibles soluciones si la evolucion es mejor que la solucion actual
        %Mutacion
        I = randperm(N);
        
        I(I == i) = [];

        r1 = I(1);
        r2 = I(2);
        r3 = I(3);
        r4 = I(4);
        r5 = I(5);
        
        %DE/rand/1/exp
        v = x(:, r1) + F*(x(:, r2) - x(:, r3)); %vector mutado
        
        %DE/rand/2/exp
%         v = x(:, r1) + F*(x(:, r2) - x(:, r3)) + F*(x(:, r4) - x(:, r5));
        
        %DE/best/1/exp
%         [~, best] = min(fitness);
%         v = x(:, best) + F*(x(:, r1) - x(:, r2));

        %Recombinacion Exponencial
        u = x(:, i);
        j = randi([1, D]);
        L=1; %Variable de control por si copia todo el arreglo
        
        while rand()<=CR && L<=D    %Mientras el ran no sea mayor a CR y no copie todo un arreglo
            u(j) = v(j);    %copia una pos del arreglo
            j = 1 + mod(j,D); %Se hace un mod a j para lograr un reinicio automatico de indice, por si sobrepasa D
            L = L+1;
        end
        
        %Seleccion
        fu = f(u(1), u(2));

        if fu < fitness(i)
            x(:, i) = u;
            fitness(i) = fu;
        end
    end
% %     Plot_Contour(f, x, xl, xu);
    [fx_best, I_best] = min(fitness);
    f_plot(n) = fx_best;
end    



figure
title('Gráfica en 2D','FontSize',15)
Plot_Contour(f, x(:, I_best), xl, xu)

figure
title('Gráfica en 3D','FontSize',15)
Plot_Surf(f, x(:, I_best), xl, xu)

%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)
title('Gráfica de convergencia','FontSize',15)

display(['Minimo global en x=' num2str(x(1,I_best)) ', y=' num2str(x(2, I_best)) ', f(x,y)=' num2str(fx_best)]);