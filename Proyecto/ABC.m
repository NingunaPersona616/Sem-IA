clear all
close all
clc

img = imread('Image_3.bmp');
temp = imread('Template.bmp');

img_g = rgb2gray(img);
temp_g = rgb2gray(temp);

[img_H,img_W] = size(img_g);
[temp_H,temp_W] = size(temp_g);

xl = [1 1]';
xu = [img_W-temp_W img_H-temp_H]';

D = 2;
G = 90;
N = 80;

L = 50;
Pf = 70;
Po = N - Pf;

l = zeros(1, Pf);

x = zeros(D, Pf);
fitness = zeros(1,Pf);
aptitud = zeros(1,Pf);

f_plot = zeros(1, G);
tic
for i=1:Pf
    x(:, i) = round(xl+(xu-xl).*rand(D, 1));
    fitness(i) = NCC(img_g, temp_g, x(1,i), x(2,i));
    
%     Modifique un poco el calculo de la aptitud debido a que daba valores muy pequenos para fitness muy buenos(altos)
    if fitness(i) >= 0
        aptitud(i) = 1+fitness(i);
    else
        aptitud(i) = abs(fitness(i));
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
        
        %Esquema de recalculo de penalizacion PD:Ocupa usar la funcion original
        for j=1:D
            if xl(j) < v(j) && v(j) < xu(j)
                %La solucion es la misma
            else
                v(j) = xl(j) + (xu(j) - xl(j))*rand();  %Recalculas el eje que se sale del espacio de busqueda
            end
        end

        v = round(v);
        fv = NCC(img_g, temp_g, v(1), v(2));

        if fv > fitness(i)
             x(:, i) = v;
             fitness(i) = fv;
             l(i) = 0;
        else
            l(i) = l(i)+1;
        end
        
        %Actualiza la aptitud para no generar mas problemas
        if fitness(i) >= 0
            aptitud(i) = 1+fitness(i);
        else
            aptitud(i) = abs(fitness(i));
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

        %Esquema de recalculo de penalizacion PD:Ocupa usar la funcion original
        for j=1:D
            if xl(j) < v(j) && v(j) < xu(j)
                %La solucion es la misma
            else
                v(j) = xl(j) + (xu(j) - xl(j))*rand();  %Recalculas el eje que se sale del espacio de busqueda
            end
        end

        v = round(v);
        fv = NCC(img_g, temp_g, v(1), v(2));

        if fv > fitness(m)
             x(:, m) = v;
             fitness(m) = fv;
             l(m) = 0;
        else
            l(m) = l(m)+1;
        end
        
        %Actualiza la aptitud para no generar mas problemas
        if fitness(i) >= 0
            aptitud(i) = 1+fitness(i);
        else
            aptitud(i) = abs(fitness(i));
        end
    end
    
    %Abejas exploradoras
    for i=1:Pf
        if l(i) > L     %Si ya tuvo muchos intentos para mejorar y no mejoro
            x(:, i) = round(xl+(xu-xl).*rand(D, 1));   %se mueve a un punto al azar
            fitness(i) = NCC(img_g, temp_g, x(1,i), x(2,i));

            %Recalculo de la aptitud
            if fitness(i) >= 0
                aptitud(i) = 1+fitness(i);
            else
                aptitud(i) = abs(fitness(i));
            end

            l(i)=0; %Se reinicia su contador de intentos
        end
    end


    [fx_best, I_best] = max(fitness);
    f_plot(g) = fx_best;
end    
toc

figure
fitness(I_best)
imshow(img)
line([x(1, I_best) x(1, I_best)+temp_W], [x(2, I_best) x(2, I_best)],'Color','g','LineWidth',3);
line([x(1, I_best) x(1, I_best)], [x(2, I_best) x(2, I_best)+temp_H],'Color','g','LineWidth',3);
line([x(1, I_best)+temp_W x(1, I_best)+temp_W], [x(2, I_best) x(2, I_best)+temp_H],'Color','g','LineWidth',3);
line([x(1, I_best) x(1, I_best)+temp_W], [x(2, I_best)+temp_H x(2, I_best)+temp_H],'Color','g','LineWidth',3);

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