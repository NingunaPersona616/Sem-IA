clear all
close all
clc
%-----------------------------(mu+1)-ES----------------------------------

% f = @(x,y) x.*exp(-x.^2 -y.^2);
% xl = [-2 -2]';
% xu = [2 2]';

f = @(x,y) (x-2).^2 + (y-2).^2;
xl = [-5 -5]';
xu = [5 5]';

G = 1500;
D = 2;
mu = 30;

%mu + 1 porque se reserva la ultima celda para el hijo que se generara
x = zeros(D, mu+1);
sigma = zeros(D, mu+1);
fitness = zeros(1, mu+1);
f_plot = zeros(1, G);

for i=1:mu
    x(:,i) = xl + (xu -xl).*rand(D,1);
    sigma(:,i) = 0.2.*rand(D,1);   %Mientras mas grande sea el escalar mayor sera la exploracion, mientras sea menor la explotacion aumentara
end

for t=1:G
    %Seleccion
    r1 = randi([1 mu]);
    r2 = r1;

    while r2 == r1
        r2 = randi([1 mu]);
    end
    
%   --------Recombinacion---------

    %Sexual Discreta
%     for j=1:D   %Toma un valor del padre o la madre hasta completar el total de dimensiones
%         if randi([0 1])
%             x(j, mu+1) = x(j, r1);  %Mu+1 es la celda del hijo
%             sigma(j, mu+1) = sigma(j, r1);
%         else
%             x(j, mu+1) = x(j, r2);
%             sigma(j, mu+1) = sigma(j, r2);
%         end
%     end

%     Sexual Intermedia
    x(:, mu+1) = (x(:, r1)+ x(:, r2))/2;  %Mu+1 es la celda del hijo
    sigma(:, mu+1) = (sigma(:, r1)+sigma(:, r2))/2;

%     Sexual Global
%     for j=1:D   %Toma un valor de cualquier padre hasta completar las dimensiones del vector del hijo
%         r1 = randi([1 mu]);
%         x(j, mu+1) = x(j, r1);  %Mu+1 es la celda del hijo
%         sigma(j, mu+1) = sigma(j, r1);
%     end
    

%---------Mutacion---------
    r = normrnd(0, sigma(:, mu+1)); %Toma un vector aleatorio de una distribucion normal respecto a las varianzas del hijo
    x(:, mu+1) = x(:, mu+1) + r;    %se suma el vector al hijo
    
    %Eliminacion del peor elemento
    for i=1:mu+1
        fitness(i) = f(x(1,i), x(2,i)); %Se calculan las soluciones de todos los elementos
    end
    %Se sortean del mejor al peor y se recupera la lista de indices ya ordenados
    [~, I] = sort(fitness);

    x = x(:, I);   %Se sortean las listas utliziando la lista de indices anteriormente ordenada
    sigma = sigma(:, I);
    fitness = fitness(I);
%     Plot_Contour(f, x(:, 1:mu), xl, xu);
    f_plot(t) = fitness(1);
end

%El mejor de todos, es el primero en la lista ordenada de x
xb = x(:, 1);  

figure
title('Gráfica en 2D','FontSize',15)
Plot_Contour(f, xb, xl, xu)

figure
title('Gráfica en 3D','FontSize',15)
Plot_Surf(f, xb, xl, xu)

%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)
title('Gráfica de convergencia','FontSize',15)

display(['Minimo global en x=' num2str(xb(1)) ', y=' num2str(xb(2)) ', f(x,y)=' num2str(f(xb(1), xb(2)))]);