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