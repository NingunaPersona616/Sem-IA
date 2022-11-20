clear all
close all
clc

% Drop-Wave
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2));

% Rastrigin
% f = @(x,y) 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y);

%McCormick
% f = @(x,y) sin(x+y)+(x-y).^2-1.5*x+2.5*y+1;
% xl = [-1.5 -3]';
% xu = [4 4]';

% f = @(x,y) (x-2).^2 + (y-2).^2;

%Ackley
% f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y))) + 20+exp(1); 

%Sistema 1
% A = [2 3; 5 1];
% b = [7; 11];

%Sistema 2
% A = [2 3 4; 1 2 3; 5 1 0; 3 4 1];
% b = [19; 13; 11; 13];

%Sistema 3
% A = [2 3; 5 4; 2 5; 4 1; 0.5 0.5];
% b = [-5; 5; -15; 15; 0];

%Sistema 4
A = [1 -2 1; 2 2 0; 5 -3 4];
b = [5; 7; 1];

%Funcion Objetivo para sistemas de ecuaciones lineales(enfoque cuadratico)
M = numel(b);   %use M poruqe m ya la usa el algoritmo
f = @(x) (1/(2*M))*sum((b-A*x).^2);

% limites inferior y superior, modificar cada que se necesiten mas o menos variables
xl = [-10 -10 -10]';
xu = [10 10 10]';
D = 3;
G = 150;
N = 250;

L = 70;
Pf = 190;
Po = N - Pf;

l = zeros(1, Pf);

x = zeros(D, Pf);
fitness = zeros(1,Pf);
aptitud = zeros(1,Pf);

f_plot = zeros(1, G);

for i=1:Pf
    x(:, i) = xl+(xu-xl).*rand(D, 1);
    fitness(i) = f(x(:,i));

    if fitness(i) >= 0
        aptitud(i) = 1/(1+fitness(i));
    else
        aptitud(i) = 1+abs(fitness(i));
    end
end  

for g=1:G
    %Abejas empleadas
    for i=1:Pf
        k=i;

        while k ==i
            k = randi([1,Pf]);
        end

        j = randi([1 D]);
        phi = 2*rand()-1;
        v = x(:,i);
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));

        fv = f(v);

        if fv < fitness(i)
             x(:, i) = v;
             fitness(i) = fv;
             l(i) = 0;
        else
            l(i) = l(i)+1;
        end
        
        %Actualiza la aptitud para no generar mas problemas
        if fitness(i) >= 0
            aptitud(i) = 1/(1+fitness(i));
        else
            aptitud(i) = 1+abs(fitness(i));
        end
    end
    
    %Abejas Observadoras
    for i=1:Po
        %Seleccion de abeja trabajadora  
        m = Ruleta(aptitud);

        k=m;

        while k == m
            k = randi([1,Pf]);
        end

        j = randi([1 D]);
        phi = 2*rand()-1;

        v = x(:,m);
        v(j) = x(j,m) + phi*(x(j,m)-x(j,k));

        fv = f(v);

        if fv < fitness(m)
             x(:, m) = v;
             fitness(m) = fv;
             l(m) = 0;
        else
            l(m) = l(m)+1;
        end
        
        %Actualiza la aptitud para no generar mas problemas
        if fitness(m) >= 0
            aptitud(m) = 1/(1+fitness(m));
        else
            aptitud(m) = 1+abs(fitness(m));
        end
    end
    
    %Abejas exploradoras
    for i=1:Pf
        if l(i) > L     %Si ya tuvo muchos intentos para mejorar y no mejoro
            x(:, i) = xl+(xu-xl).*rand(D, 1);   %se mueve a un punto al azar
            fitness(i) = f(x(:,i));
        
            if fitness(i) >= 0
                aptitud(i) = 1/(1+fitness(i));
            else
                aptitud(i) = 1+abs(fitness(i));
            end

            l(i)=0; %Se reinicia su contador de intentos
        end
    end


    [fx_best, I_best] = min(fitness);
    f_plot(g) = fx_best;
end    

%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)
title('Gráfica de convergencia','FontSize',15)

display(['Valor de las variables:'])
for j=1:D
    display(['X_' num2str(j) ' = ' num2str(x(j, I_best))] );
end

display([' ']);
display(['---------------------']);

display(['Error = ' num2str(fitness(I_best))])

display([' ']);
display(['---------------------']);
display(['Valor de b candidata:'])
A*x(:, I_best)

%Funcioneeees
function k = Ruleta(aptitud)
    p = aptitud/sum(aptitud);   %Se calculan las probabilidades de cada individuo
    r = rand();                 %Se calcula un porcentaje random para que lo busque la ruleta

    p_sum = 0;                  %Variable donde se guardara la acumulacion de porcentajes

    N = numel(aptitud);
    k = N;

    for i=1: N
        p_sum = p_sum + p(i);   %Se va aumentando con el porcentaje de cada individuo

        if(p_sum >= r)      %Si el el porcentaje del ultimo individuo hace que la sumatoria sea mayor que el % random
            k = i;          %Se selecciona ese individuo
            return
        end
    end
end