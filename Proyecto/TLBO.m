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
G = 120;
N = 50;

x = zeros(D, N);
fitness = zeros(1,N);

f_plot = zeros(1, G);
tic
for i=1:N
    x(:, i) = round(xl+(xu-xl).*rand(D, 1));
    fitness(i) = NCC(img_g, temp_g, x(1,i), x(2,i)); %Aqui va el NCC
end  

for g=1:G
    
    for i=1:N
        %Fase de ensenanza
        [~, t] = max(fitness); %se selecciona un maestro
        Tf = randi([1,2]);  %Factor de ensenanza
        c = zeros(D, 1);    %Nueva solucion generada

        for j=1:D
            x_mean = mean(x(j, :)); %saca el promedio de todos los elementos de cada materia j, siendo j la materia
             r = rand();
             c(j) = x(j,i) + r*(x(j,t) - Tf*x_mean); %calcula la nueva solucion usando (la pos del maestro - el promedio) + la sol actual
        end
        
        %Esquema de recalculo de penalizacion PD:Ocupa usar la funcion original
        for j=1:D
            if xl(j) < c(j) && c(j) < xu(j)
                %La solucion es la misma
            else
                c(j) = xl(j) + (xu(j) - xl(j))*rand();  %Recalculas el eje que se sale del espacio de busqueda
            end
        end

        c = round(c);
        fc = NCC(img_g, temp_g, c(1), c(2));

        if fc > fitness(i)  %Si la nueva solucion candidata es mejor a la solucion actual
            x(:, i) = c;    %Cambia la sol actual por la sol candidata
            fitness(i) = fc;
        end


        %Fase de Aprendizaje
        k = i;

        while k==i
            k = randi([1, N]);
        end

        c = zeros(D, 1);

         if fitness(i) > fitness(k) %si el alumno puede aprender de otro alumno, acerca al alumno a la solucion objetivo
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

        %Esquema de recalculo de penalizacion PD:Ocupa usar la funcion original
        for j=1:D
            if xl(j) < c(j) && c(j) < xu(j)
                %La solucion es la misma
            else
                c(j) = xl(j) + (xu(j) - xl(j))*rand();  %Recalculas el eje que se sale del espacio de busqueda
            end
        end
        
        c = round(c);
        fc = NCC(img_g, temp_g, c(1), c(2));

        if fc > fitness(i)  %Si la nueva solucion candidata es mejor a la solucion actual
            x(:, i) = c;    %Cambia la sol actual por la sol candidata
            fitness(i) = fc;
        end

    end
    
    [fx_best, I_best] = max(fitness);
    f_plot(g) = fx_best;
end
toc

%Grafica de convergencia
% figure
% plot(f_plot, 'b-', 'LineWidth', 2)
% xlabel('iteracion','FontSize',15)
% ylabel('fx','FontSize',15)
% title('Gráfica de convergencia','FontSize',15)



figure
fitness(I_best)
imshow(img)
line([x(1, I_best) x(1, I_best)+temp_W], [x(2, I_best) x(2, I_best)],'Color','g','LineWidth',3);
line([x(1, I_best) x(1, I_best)], [x(2, I_best) x(2, I_best)+temp_H],'Color','g','LineWidth',3);
line([x(1, I_best)+temp_W x(1, I_best)+temp_W], [x(2, I_best) x(2, I_best)+temp_H],'Color','g','LineWidth',3);
line([x(1, I_best) x(1, I_best)+temp_W], [x(2, I_best)+temp_H x(2, I_best)+temp_H],'Color','g','LineWidth',3);