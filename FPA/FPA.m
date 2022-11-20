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
m = numel(b);
f = @(x) (1/(2*m))*sum((b-A*x).^2);

% limites inferior y superior, modificar cada que se necesiten mas o menos variables
xl = [-10 -10 -10]';
xu = [10 10 10]';

D = 3;
G = 150;
N = 250;

x = zeros(D, N);
fitness = zeros(1,N);

f_plot = zeros(1, G);

for i=1:N
    x(:, i) = xl+(xu-xl).*rand(D, 1);
    fitness(i) = f(x(:,i));
end  

p = 0.8;
lambda = 1.5;
sigma2 = (((gamma(1+lambda))/(lambda*gamma((1+lambda)/2))) * ((sin((pi*lambda)/2))/(2^((lambda-1)/2))))^(1/lambda); %Sigma pal vuelo de levy


for g=1:G
    [~, igb] = min(fitness);    %Saca la mejor sol global

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

        fy = f(y);

        if fy < fitness(i)
            x(:, i) = y;
            fitness(i) = fy;
        end
    end

    f_plot(g) = min(fitness);
end

[~, I_best] = min(fitness);



%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)
title('Gráfica de convergencia','FontSize',15)

% display(['Minimo global en x=' num2str(x(1,I_best)) ', y=' num2str(x(2, I_best)) ', f(x,y)=' num2str(fitness(I_best))]);


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

