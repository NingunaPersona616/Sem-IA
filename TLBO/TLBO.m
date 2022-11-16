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

f = @(x,y) (x-2).^2 + (y-2).^2;

%Ackley
% f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y))) + 20+exp(1); 

xl = [-5 -5]';
xu = [5 5]';

D = 2;
G = 150;
N = 50;

x = zeros(D, N);
fitness = zeros(1,N);

f_plot = zeros(1, G);

for i=1:N
    x(:, i) = xl+(xu-xl).*rand(D, 1);
    fitness(i) = f(x(1, i), x(2, i));
end  

for g=1:G
    
    for i=1:N
        %Fase de ensenanza
        [~, t] = min(fitness); %se selecciona un maestro
        Tf = randi([1,2]);  %Factor de ensenanza
        c = zeros(D, 1);    %Nueva solucion generada

        for j=1:D
            x_mean = mean(x(j, :)); %saca el promedio de todos los elementos de cada materia j, siendo j la materia
             r = rand();
             c(j) = x(j,i) + r*(x(j,t) - Tf*x_mean); %calcula la nueva solucion usando (la pos del maestro - el promedio) + la sol actual
        end

        fc = f(c(1), c(2));

        if fc < fitness(i)  %Si la nueva solucion candidata es mejor a la solucion actual
            x(:, i) = c;    %Cambia la sol actual por la sol candidata
            fitness(i) = fc;
        end


        %Fase de Aprendizaje
        k = i;

        while k==i
            k = randi([1, N]);
        end

        c = zeros(D, 1);

         if fitness(i) < fitness(k) %si el alumno puede aprender de otro alumno, acerca al alumno a la solucion objetivo
            for j=1:D
                r = rand();
                c(j) = x(j, i) + r*(x(j, i) - x(j, k)) ; %calcula la nueva sol con la sol actual y suma la diferencia con el otro alumno
            end
         else   %Aleja ala alumno de la solucion objetivo
             for j=1:D
                r = rand();
                c(j) = x(j, i) + r*(x(j, k) - x(j, i)) ; %calcula la nueva sol con la sol actual y resta la diferencia con el otro alumno
            end
         end
        
        fc = f(c(1), c(2));

        if fc < fitness(i)  %Si la nueva solucion candidata es mejor a la solucion actual
            x(:, i) = c;    %Cambia la sol actual por la sol candidata
            fitness(i) = fc;
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
