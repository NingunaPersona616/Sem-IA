clear all
close all
clc

% Drop-Wave
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2));

% Rastrigin
% f = @(x,y) 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y);

%McCormick
f = @(x,y) sin(x+y)+(x-y).^2-1.5*x+2.5*y+1;
xl = [-1.5 -3]';
xu = [4 4]';

% f = @(x,y) (x-2).^2 + (y-2).^2;

%Ackley
% f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y))) + 20+exp(1); 

% xl = [-5 -5]';
% xu = [5 5]';

D = 2;
G = 150;
N = 50;

L = 40;
Pf = 30;
Po = N - Pf;

l = zeros(1, Pf);

x = zeros(D, Pf);
fitness = zeros(1,Pf);
aptitud = zeros(1,Pf);

f_plot = zeros(1, G);

for i=1:Pf
    x(:, i) = xl+(xu-xl).*rand(D, 1);
    fitness(i) = f(x(1, i), x(2, i));

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

        fv = f(v(1), v(2));

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

        fv = f(v(1), v(2));

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
            fitness(i) = f(x(1, i), x(2, i));
        
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

    Plot_Contour(f, x, xl, xu)
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