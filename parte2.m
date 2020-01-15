for j = 1:tot
    
    % Con eq.11
    freqs = psd{j}.get_mode_freq(); %#ok<*SUSENS>
    h = [];
    for i = 1:min(length(freqs), 3) % Recorre cada modo
        wi = freqs(i) * 2 * pi;
        h = [h, (wi * cable_largo / (i * pi))^2 * cable_densidad]; %#ok<*AGROW>
    end
    h = h .* 0.10197162; % kg -> kgf
    hprom = mean(h);
    
    h11 = h;
    
    % Con eq.12
    g = 9.80665;
    A = (cable_diametro / 2)^2 * pi; % m2
    E = 27.67507 / (0.7 * A) * g; % kgf/m
    m = cable_densidad;
    L = cable_largo;
    
    h4 = [];
    for i = 1:min(length(freqs), 3)
        if 2 * i - 1 > length(freqs)
            break
        end
        
        % Prueba varios valores de h, cercanos al anterior
        eps = 0.25;
        h = linspace(hprom*(1 - eps), hprom*(1 + eps), 1000);
        err = zeros(1, 1000);
        for k = 1:length(h)
            wi = freqs(2*i-1) * 2 * pi;
            d = m * g * cable_largo^2 / (8 * h(k) * g);
            lamb2 = ((8 * d / L)^2 * E * A) / (h(k) * g) * 10^3;
            win = wi * L / sqrt(h(k)*g/m); % Normalizado
            err(k) = tan(win/2) - (win / 2) + (4 / lamb2) * (win / 2)^3;
        end
        
        % Busca el maximo gradiente de cambio, ese punto corresponde al
        % valor de h
        grad = 0;
        indx = -1;
        for k = 2:length(h)
            gk = abs(err(k)-err(k-1));
            if gk > grad
                grad = gk;
                indx = k;
            end
        end
        
        h4 = [h4, (h(indx) + h(indx-1)) / 2];
        
    end
    
    h4
    h11
    
end