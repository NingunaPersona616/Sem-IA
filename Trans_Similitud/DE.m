clear all
close all
clc

fp = @(x, xl, xu) f(x(1), x(2)) + 1000*Penalty(x, xl, xu);

%%%
img_ref = imread("ref_3.png");
[~, ~, P] = readBarcode(img_ref, "QR-CODE");
[Nu, Mu, ~] = size(img_ref);    %N upper y M upper

img_des = imread("letsgous.jpg");
[n, m, ~] = size(img_des);
%%%

%Puntos de referencia del QR
X1 = P(1, :)';
X2 = P(2, :)';
X3 = P(3, :)';

%Puntos de referencia de la imagen a sobreponer
x1 = [1, n]';   %esquina izq inferior
x2 = [1, 1]';   %Esquina izq superior
x3 = [m, 1]';   %Esquina derecha superior

% Parametros para optimizar
dx = 1;
dy = 1;
theta = 0;
s = 0.5;

xl = [1; 1; -pi; 0];
xu = [Nu; Mu; pi; 10];

D = 4;
G = 150;
N = 50;
F = 0.6;
CR = 0.9;

x = zeros(D, N);
fitness = zeros(1,N);
f_plot = zeros(1, G);

for i=1:N
    x(:, i) = xl+(xu-xl).*rand(D, 1);
    q = x;

    xp1 = Transformacion_Similitud(q,x1);
    xp2 = Transformacion_Similitud(q,x2);
    xp3 = Transformacion_Similitud(q,x3);
    
    e1 = Distancia_Euclidiana(X1,xp1);
    e2 = Distancia_Euclidiana(X2,xp2);
    e3 = Distancia_Euclidiana(X3,xp3);
    
    f = (1/6)*(e1^2+e2^2+e3^2);
    fitness(i) = f;
end  

for n=1:G
    for i=1:N   %Por cada generacion se actualizan las posibles soluciones si la evolucion es mejor que la solucion actual
        %Mutacion
        r1 = i;
        while r1 == i
            r1 = randi([1,N]);
        end
    
        r2 = r1;
        while r2==i || r2==r1
            r2 = randi([1,N]);
        end
    
        r3 = r2;
        while r3==i || r3==r1 || r3==r2
            r3 = randi([1,N]);
        end

        v = x(:, r1) + F*(x(:, r2) - x(:, r3)); %vector mutado

        %Recombinacion
        u = zeros(D, 1);

        for j=1:D
            if rand() <= CR %Si tiene probabilidad recombinarse con el vector mutado
                u(j) = v(j);
            else
                u(j) = x(j, i);
            end
        end
        
        %Seleccion
        q = u;

        xp1 = Transformacion_Similitud(q,x1);
        xp2 = Transformacion_Similitud(q,x2);
        xp3 = Transformacion_Similitud(q,x3);
        
        e1 = Distancia_Euclidiana(X1,xp1);
        e2 = Distancia_Euclidiana(X2,xp2);
        e3 = Distancia_Euclidiana(X3,xp3);
        
        fu = (1/6)*(e1^2+e2^2+e3^2) + 1000*(Penalty(q, xl, xu)); %Funcion de penalizacion

        if fu < fitness(i)
            x(:, i) = u;
            fitness(i) = fu;
        end
    end
end    

[f_best, I_best] = min(fitness);
display(['Error:' num2str(f_best)]);
display(['Minimo global en']);
x(:,I_best)
Imprimir_Imagenes(x(:,I_best),img_des,img_ref)

%Metodo 1 de penalizacion
function z = Penalty(x, xl, xu)
    z = 0;
    D = numel(xl);

    for j=1:D
        if xl(j) < x(j) && x(j) < xu(j)
            z = z + 0;
        else
            z = z + 1;
        end
    end
end

%Metodo 2 de penalizacion
% function z = Penalty(x, xl, xu)
%     z = 0;
%     D = numel(xl);
% 
%     for j=1:D
%         if xl(j) < x(j)
%             z = z + 0;
%         else
%             z = z + (xl(j) - x(j))^2;
%         end
% 
%         if x(j) < xu(j)
%             z = z + 0;
%         else
%             z = z + (xu(j) - x(j))^2;
%         end
%     end
% end

%Funciones
function xp = Transformacion_Similitud (qi,xi)
    dx = qi(1);
    dy = qi(2);
    theta = qi(3);
    s = qi(4);
    
    xp = [s*cos(theta) -s*sin(theta); s*sin(theta) s*cos(theta)]*xi + [dx dy]';
end

function e = Distancia_Euclidiana (X,x)
    e = sqrt((X(1)-x(1))^2+(X(2)-x(2))^2);
end

function Imprimir_Imagenes (q,img_des,img_ref)
    dx = q(1);
    dy = q(2);
    theta = q(3);
    s = q(4);
    
    T = [s*cos(theta) -s*sin(theta) dx; s*sin(theta) s*cos(theta) dy; 0 0 1];
    Tp = projective2d(T');
    
    [N,M,~] = size(img_ref);
    [n,m,~] = size(img_des);

    panoramaView = imref2d([N M]);
    Iwarp = imwarp(img_des,Tp,'OutputView',panoramaView);
    Imask = imwarp(true(n,m),Tp,'OutputView',panoramaView);
    
    blender = vision.AlphaBlender('Operation','Binary mask','MaskSource','Input port');
    panorama = step(blender,img_ref,Iwarp,Imask);
    
    imshow(panorama)
end