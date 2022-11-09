clear all
close all
clc
%------------------------------(1+1)-ES-----------------------------

% f = @(x,y) x.*exp(-x.^2 -y.^2);
% xl = [-2 -2]';
% xu = [2 2]';

f = @(x,y) (x-2).^2 + (y-2).^2;
xl = [-5 -5]';
xu = [5 5]';

G = 300;
D = 2;
sigma = 0.5;

x = xl + (xu -xl).*rand(D,1);
f_plot = zeros(1, G);

for t=1:G
    r = normrnd(0, sigma, [D, 1]);  %Genera un arreglo aleatorio de D renglones y 1 columa (vector vertical), usando valores de uuna distribucion normal N(Media, Varianza)
    
    %Mutacion
    y = x + r;

    if f(y(1), y(2)) < f(x(1), x(2))    %si la solucion del hijo mejora a la del padre
        x = y;  %Entonces ahora el hijo se conviete ahora en el padre
    end
    f_plot(t) = f(x(1), x(2));
    %Plot_Contour(f, [x y], xl, xu);
end


figure
title('Gráfica en 2D','FontSize',15)
Plot_Contour(f, x, xl, xu)

figure
title('Gráfica en 3D','FontSize',15)
Plot_Surf(f, x, xl, xu)

%Grafica de convergencia
figure
plot(f_plot, 'b-', 'LineWidth', 2)
xlabel('iteracion','FontSize',15)
ylabel('fx','FontSize',15)
title('Gráfica de convergencia','FontSize',15)


display(['Minimo global en x=' num2str(x(1)) ', y=' num2str(x(2)) ', f(x,y)=' num2str(f(x(1), x(2)))]);