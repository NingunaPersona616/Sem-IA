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
G = 150;
N = 50;

x = zeros(D, N);
fitness = zeros(1,N);

f_plot = zeros(1, G);

p = 0.8;
lambda = 1.5;
sigma2 = (((gamma(1+lambda))/(lambda*gamma((1+lambda)/2))) * ((sin((pi*lambda)/2))/(2^((lambda-1)/2))))^(1/lambda); %Sigma pal vuelo de levy

tic
for i=1:N
    x(:, i) = round(xl+(xu-xl).*rand(D, 1));
    fitness(i) = NCC(img_g, temp_g, x(1,i), x(2,i)); %Aqui va el NCC
end  

for g=1:G

    [~, igb] = max(fitness);    %Saca la mejor sol global

    for i=1:N
        %Polinizacion global
        if rand() < p
            u = normrnd(0, sigma2, [D, 1]) ; %Genera un vector v (con D reng y 1 col) usando una dist normal
            v = normrnd(0, 1, [D,1]);   %genera un vector igual al enterior pero con desv estdr=1 
            L = u./(abs(v).^(1/lambda));    %Paso de polinizacion con el vuelo de Levy
            
            y = x(:, i ) + L.*(x(:, igb) - x(:,i)); %Calcula una sol candidata con un error del global y el actual multiplicadas por el paso de polinizacion

        else  %Polinizacion local      
            j=i;
            while j==i
                j = randi([1, N]);
            end
            
            k=j;
            while k==i && k==j
                k = randi([1, N]);
            end

            y = x(:, i) + rand()*(x(:, j) - x(:, k));
        end

        %Esquema de recalculo de penalizacion PD:Ocupa usar la funcion original
        for j=1:D
            if xl(j) < y(j) && y(j) < xu(j)
                %La solucion es la misma
            else
                y(j) = xl(j) + (xu(j) - xl(j))*rand();  %Recalculas el eje que se sale del espacio de busqueda
            end
        end
        
        %Redondeamos y para no generar errores con los indices
        y = round(y);

        fy = NCC(img_g, temp_g, y(1), y(2)); %Aqui va el NCC

        if fy > fitness(i)
            x(:, i) = y;
            fitness(i) = fy;
        end
    end

    f_plot(g) = max(fitness);
end
toc

[~, I_best] = max(fitness);


figure
hold on
fitness(I_best)
imshow(img)
line([x(1, I_best) x(1, I_best)+temp_W], [x(2, I_best) x(2, I_best)],'Color','g','LineWidth',3);
line([x(1, I_best) x(1, I_best)], [x(2, I_best) x(2, I_best)+temp_H],'Color','g','LineWidth',3);
line([x(1, I_best)+temp_W x(1, I_best)+temp_W], [x(2, I_best) x(2, I_best)+temp_H],'Color','g','LineWidth',3);
line([x(1, I_best) x(1, I_best)+temp_W], [x(2, I_best)+temp_H x(2, I_best)+temp_H],'Color','g','LineWidth',3);