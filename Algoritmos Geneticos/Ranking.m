function k = Ranking(aptitud)
    [~, I] = sort(aptitud, 'descend');  %Se ordenan los individuos de mejor a peor
    N = numel(aptitud);     %se obtiene la cantidad de individuos que estaran en el ranking

    rank = N:-1:1;  %Lista de pesos de cada individuo(relacionada a la probabilidad por su rank)
    rank_total = sum(rank); %Total (ayuda para sacar la probabilidad)

    p_sum = 0;
    r = rand();

    for i=1:N   %Se suman las probabilidades hasta que se llegue o se pase del porcentaje r
        p_sum = p_sum + rank(i) / rank_total;   %suma la probabilidad acumulada + la probabilidad por ranking del individuo actual

        if p_sum >= r   %Si se llega o se pasa del porcentaje  r se regresa el individuo actual
            k = I(i);   
            return
        end
    end
end