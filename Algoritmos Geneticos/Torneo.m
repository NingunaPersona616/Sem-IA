function k = Torneo (aptitud)
    N = numel(aptitud); %Total de individuos

    tao = round(N*0.3); %porcentaje de individuos seleccionados

    I = randi(N, [1 tao]);  %Lista de inidividuos seleccionados(se escoge la cantidad representada por el porcentaje)

    [~,i] = max(aptitud(I)); %De los individuos seleccionados se escoge el de mayor aptitud

    k = I(i);
end
