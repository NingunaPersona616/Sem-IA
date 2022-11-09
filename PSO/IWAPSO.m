%--------------------------------- IWAPSO --------------------------
clear all
close all
clc

% Drop-Wave
% f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2));

% Rastrigin
% f = @(x,y) 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y);

%Ackley
f = @(x,y) -20*exp(-0.2*sqrt(0.5*(x.^2 + y.^2))) - exp(0.5*(cos(2*pi*x)+cos(2*pi*y))) + 20+exp(1); 

xl = [-5 -5]';
xu = [5 5]';

D = 2;
G = 150;
N = 50;

x = zeros(D, N);
xb = zeros(D, N);
v = zeros(D, N);
fitness = zeros(1,N);
f_plot = zeros(1, N);

w_max = 0.8;
w_min = 0.1;
c1 = 2; %Ajusta el conocimiento cognitivo(el de la propia particula)
c2 = 2; %Ajusta el conocimiento social(el global del enjambre)

for i=1:N
    x(:, i) = xl+(xu-xl).*rand(D, 1);
    v(:, i) = 0.3*randn(D, 1);  %Numero randon de una distribucion normal multiplicado para que no tome valores grandes
    xb(:, i) =  x(:, i);
    fitness(i) = f(x(1, i), x(2, i));
end    


for g=1:G
    %Plot_Contour(f, x, xl, xu);
    
    for i=1:N
        fx = f(x(1, i), x(2, i));

        if(fx < fitness(i)) %si el fx actual de la particula mejora el mejor fx de la misma particula
            xb(:, i) = x(:, i); %fxi y xi se convierten ahora en los mejores
            fitness(i) = fx;
        end
    end

    %Elegir mejor particula del enjambre

    [fx_best, I_best] = min(fitness);
    f_plot(g) = fx_best;
    
    %Inertia Weight
    w = w_max - (g/G)*(w_max-w_min);

    for i=1:N
        v(:, i) = w*v(:, i) + rand()*c1*(xb(:, i) - x(:, i)) + rand()*c2*(xb(:, I_best) - x(:, i)); %Actualiza la velocidad de la particula
        x(:, i) =  x(:, i) + v(:, i);
    end
end  

figure
title('Gráfica en 2D','FontSize',15)
Plot_Contour(f, xb(:, I_best), xl, xu)

figure
title('Gráfica en 3D','FontSize',15)
Plot_Surf(f, xb(:, I_best), xl, xu)

%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)
title('Gráfica de convergencia','FontSize',15)

display(['Minimo global en x=' num2str(xb(1,I_best)) ', y=' num2str(xb(2, I_best)) ', f(x,y)=' num2str(fx_best)]);