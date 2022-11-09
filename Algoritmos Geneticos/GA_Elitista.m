clear all
close all
clc

% f = @(x,y) x.*exp(-x.^2 -y.^2);
% xl = [-2 -2]';
% xu = [2 2]';

f = @(x,y) (x-2).^2 + (y-2).^2;
xl = [-5 -5]';
xu = [5 5]';

D = 2;  %Dimension del problema
G = 100; %Numero de generaciones
N = 100; %Tamaño de la poblacion
pm=0.08;  %Probabilidad de mutacion
E = 10;  %Cantidad de elitistas

%Un poco de inicializacion
x = zeros(D, N);    %Se llena una matriz de posibles soluciones(individuos) en este caso de 2 renglones por 10 columnas
                    %Primer renglon para los valore de x y el segundo para los valores de y

fitnes = zeros(1, N);   %Arreglo donde se guardan las evaluaciones f(xi) de los individuos

aptitud = zeros(1, N);  %Arreglo donde se guarda la aptitud de los individuos 

f_plot = zeros(1, D);

for i=1: N
    x(:, i) = xl+(xu-xl).*rand(D,1);

end

%Algoritmo Genetico Clasico
for g=1:G
    %Plot_Contour(f, x, xl, xu)
    for i=1: N  %Inicializa las aptitudes de los padres
        fitnes(i) = f(x(1, i), x(2, i));    %Evalua f(xi,yi) en la columna i

        if(fitnes(i)>=0)
            aptitud(i) = 1 / (1 + fitnes(i));
    
        else
            aptitud(i) = 1 + abs(fitnes(i));
        end 
    end
    
    f_plot(g) = min(fitnes);

    y = zeros(D, N-E);    %Inicializa en ceros los hijos
    

    for i=1:2:N-E %Llenar por parejas de hijos el arreglo, osea de 2 en 2
%Seleccion
        r1 = Ruleta(aptitud);   %se busca un primer padre
%         r1 = Ranking(aptitud);
%         r1 = Torneo(aptitud);
        r2 = r1;
        
        while(r2 == r1)     %se busca un segundo padre distitno al primero
            r2 = Ruleta(aptitud);
%              r2 = Ranking(aptitud);
%              r2 = Torneo(aptitud);
        end

        padre1 = x(:, r1);
        padre2 = x(:, r2);

         %cruza
        hijo1 = padre1;
        hijo2 = padre2;
        
        pc = randi([1 D]);   %Punto de cruce toma un valor desde el indice inicial hasta el maximo de poblacion
        
        hijo1(pc:D) = padre2(pc:D); %Desde el indice pc(punto de cruce) en adelante intercambia genes del padre2
        hijo2(pc:D) = padre1(pc:D);  %Desde el indice pc(punto de cruce) en adelante intercambia genes del padre1
        
        %EJEMPLO 
        % padre1 = [0 1 2 3];
        % padre2 = [9 8 7 6];
        % 
        % pc = 2;
        % 
        % hijo1 = [0 1 7 6];
        % hijo2 = [9 8 2 3];

        %Cruza plana y Aritmetica
%         r = rand();
%         hijo1 = r*padre1 + (1-r)*padre2;
%         hijo2 = (1-r)*padre1 + r*padre2;
    
        y(:, i) = hijo1;    %Se guarda el hijo1 
        y(:, i+1) = hijo2;  %Y se guarda el hijo2 consecutivamente
    end
    
    %Mutacion
    for i=1:N-E
        for j=1:D
            
            if rand() < pm
                %y(j, i) = xl(j)+(xu(j)-xl(j))*rand();
                %y(j, i) = y(j, i) + normrnd(0,0.3); %Distribucion Normal para la primer funcion
                y(j, i) = y(j, i) + normrnd(0,1);   %Distribucion normal 2da funcion
            end
        end
    end    
    
    %Seleccion de elitistas
    [~, I] = sort(aptitud,'descend');

    %x(:, I(1,E));  todos los renglones del array de X en las columnas
    %dadas por el array de elitistas
    
    x=[y x(:, I(1:E))];    %Los hijos y los elitistas sustituyen a la generacion anterior
end



%Calculamos de nuevo las aptitudes de los ultimos padres
for i=1: N  %Inicializa las aptitudes de los padres
    fitnes(i) = f(x(1, i), x(2, i));    %Evalua f(xi,yi) en la columna i

    if(fitnes(i)>=0)
        aptitud(i) = 1 / (1 + fitnes(i));

    else
        aptitud(i) = 1 + abs(fitnes(i));
    end 
end

[~, i_mejor] = max(aptitud);

figure
title('Gráfica en 2D','FontSize',15)
Plot_Contour(f, x, xl, xu)


figure
title('Gráfica en 3D','FontSize',15)
Plot_Surf(f, x, xl, xu)


%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
title('Gráfica de convergencia','FontSize',15)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)

display(['Min Global en: x=', num2str(x(1, i_mejor)), ',  y=', num2str(x(2, i_mejor)), ', f(x,y)=' num2str(f(x(1, i_mejor), x(2, i_mejor)))])

