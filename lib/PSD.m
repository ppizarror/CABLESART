function [f_fft, f_psd, psd, fft, fftcomp, envFormaModal, tlocMean, tlocStd, locMean, ...
    locStd, locFreq_fft, maxlocs, pks, betaNodo, betaFreqNodo, fftmean, fftstd, ...
    psdmean, psdstd] = PSD(a, fs, gdl, varargin)
% PSD: Power Spectral Density. Calcula la FFT de un registro sismico
% analizado para varios nodos con distintos grados de libertad. La funcion
% permite calcular distintos tipos de periodos naturales, arrojando un
% vector. Adicionalmente se calcula una razon de amortiguamiento con cada
% modo.
%
% Input:
%   a               Registro de aceleracion
%   fs              Factor de muestreo del registro
%   gdl             Vector con los grados de libertad de los nodos en el analisis
%
% Parametros opcionales:
%   betaFFT         Realiza el calculo del amortiguamiento con FFT o PSD
%   betaFFTMax      El amortiguamiento se calcula con el maximo FFT de todos los nodos
%   peakFFT         Realiza el calculo de peaks con FFT o PSD
%   peakMinDistance Distancia minima entre peaks requerida
%   peakThreshold   Limite del peak
%   peakUseMean     Usa el promedio
%   tmin            Tiempo minimo de analisis
%   tmax            Tiempo maximo de analisis
%   tukeywinr       Factor de la ventana de tukey
%   zerofill        Indica relleno de ceros para FFT
%   usewindows      Usa wventanas en vez de [tmin, tmax]
%   usepwelch       Usa pwelch para la identificacion del PSD
%   ventana         Tamaño en segundos de cada ventana
%   factor          Factor del tamaño en PSD con pwelch
%   noverlap        Tamaño de translape entre ventanas pwelch
%
% Output:
%   f               Vector de frecuencias
%   psd             Cell PSD de cada registro de los nodos
%   fft             Cell FFT de cada registro de los nodos
%   fftcomp         Registro completo (Re+Im)
%   envFormaModal   Cell con formas modales para cada periodo
%   tlocMean        Periodos medios de cada modo obtenidos de los peaks
%   tlocStd         Desviacion estandar de cada modo obtenido de los peaks
%   locMean         Frecuencia media de cada modo obtenido en los peaks
%   locStd          Desviacion estandar de la frecuencia de cada modo
%   locFreq         Posicion en el vector de frecuencias de cada modo
%   maxlocs         Numero de modos encontrados en los peaks
%   pks             Vector de peaks
%   betaNodo        Vector de amortiguamientos por cada modo
%   betaFreqNodo    Valor del FFT objetivo por cada amortiguamiento modal
%   fftmean         Promedios FFT
%   fftstd          Desviacion estandar FFT
%   psdmean         Promedios PSD
%   psdstd          Desviacion estandar PSD

%% Parametros opcionales
p = inputParser;
p.KeepUnmatched = true;
addOptional(p, 'betaFFT', true);
addOptional(p, 'betaFFTMax', false); % Calcula el amortiguamiento con el maximo
addOptional(p, 'fmin', 0.5);
addOptional(p, 'peakFFT', true);
addOptional(p, 'peakMinDistance', 0.5); % Requerido para el calculo
addOptional(p, 'peakThreshold', 0);
addOptional(p, 'peakFreqThreshold', 0);
addOptional(p, 'peakUseMean', false);
addOptional(p, 'tmax', -1);
addOptional(p, 'tmin', 0);
addOptional(p, 'tukeywinr', 0.01);
addOptional(p, 'zerofill', 0);
addOptional(p, 'usewindows', false);
addOptional(p, 'windows', {});
addOptional(p, 'forcepeakmax', false);

% Para el PSD
addOptional(p, 'usepwelch', true);
addOptional(p, 'ventana', 10); % Segundos
addOptional(p, 'factor', 100);
addOptional(p, 'noverlap', 0.5); % Del tamaño de la ventana

parse(p, varargin{:});
r = p.Results;

if r.betaFFT
    fprintf('\tSe calculan amortiguamientos con FFT\n');
else
    fprintf('\tSe calculan amortiguamientos con PSD\n');
end
if r.peakFFT
    fprintf('\tSe calculan peaks periodos con FFT\n');
else
    fprintf('\tSe calculan peaks periodos con PSD\n');
end

% Se obtiene la ventana de tiempo
c1 = 1;
if r.tmin ~= 0
    if r.tmin < 0
        error('El tiempo inferior no puede ser cero');
    end
    c1 = fix(r.tmin*fs);
end
cend = false;
if r.tmax ~= -1
    if r.tmin >= r.tmax
        error('El tiempo inferior tmin no puede ser mayor a tmax');
    end
    c2 = fix(r.tmax*fs);
else
    c2 = -1;
    cend = true;
end

% Crea la ventana si es falso
if ~r.usewindows
    windows = {[c1, c2]};
else
    windows = r.windows;
end

% Numero de grados de libertad por el numero de ventanas
ng = length(gdl) * length(windows);

%% Calcula la FFT
fft = cell(1, ng);
psd = cell(1, ng);

% Calcula el tamaño maximo de las ventanas, FFT
nfft_fft = 0;
for win = 1:length(windows)
    a_win = a(windows{win}(1):windows{win}(2), :);
    for k = 1:length(gdl)
        fft_acc = a_win(:, gdl(k));
        fft_acc = [fft_acc; zeros(floor(r.zerofill*length(fft_acc)), 1)];
        nfft_fft = max(nfft_fft, length(fft_acc));
    end
end

% Calcula el tamaño maximo de las ventanas, PSD
window = r.ventana * fs; % Tamaño de ventana de analisis
nfft_psd = 2^nextpow2(floor(window*r.factor)); % Numero de puntos de la FFT
noverlap = floor(r.noverlap*window);

kk = 1; % Posicion en la que se guarda el registro
for win = 1:length(windows)
    
    % Corta la ventana
    a_win = a(windows{win}(1):windows{win}(2), :);
    
    % Recorre cada grado de libertad (columna del vector)
    for k = 1:length(gdl)
        
        % Obtiene la aceleracion del grado de libertad analizado
        fft_acc = a_win(:, gdl(k));
        
        % Limita la ventana
        if cend
            c2 = length(fft_acc);
        end
        if c1 > length(fft_acc)
            error('El tiempo inferior excede el largo del vector de aceleracion');
        end
        if c2 > length(fft_acc)
            error('El tiempo superior excede el largo del vector de aceleracion');
        end
        
        % Rellena con ceros
        fft_acc = [fft_acc; zeros(floor(r.zerofill*length(fft_acc)), 1)];
        tuck = tukeywin(length(fft_acc), r.tukeywinr);
        acctuck = fft_acc .* tuck;
        [f, fftt, ~] = DFT(fs, acctuck, nfft_fft);
        
        % Busca hasta el valor minimo
        ftarget = r.fmin;
        for i = 1:length(f)
            if f(i) > ftarget
                ftarget = f(i);
                break
            end
        end
        tf = find(f == ftarget);
        f_fft = f(tf:end);
        fftt = fftt(tf:end);
        fftcomp{kk} = fftt; % Guarda el registro complejo
        fftt = abs(fftt); % O si no plot reclama
        fft{kk} = fftt; % Guarda el registro
        
        % Calcula el PSD con pwelch
        if r.usepwelch
            psd_acc = a_win(:, gdl(k));
            tuck = tukeywin(length(psd_acc), r.tukeywinr);
            acctuck = psd_acc .* tuck;
            
            [psdacc, f_psd] = pwelch(acctuck, window, noverlap, nfft_psd, fs);
            
            % Busca hasta el valor minimo
            ftarget_psd = r.fmin;
            for i = 1:length(f_psd)
                if f_psd(i) > ftarget_psd
                    ftarget_psd = f_psd(i);
                    break
                end
            end
            tf = find(f_psd == ftarget_psd);
            f_psd = f_psd(tf:end); %#ok<*NASGU>
            psdacc = psdacc(tf:end);
            psd{kk} = psdacc;
            
        else % Metodología sencilla de cálculo
            psd{kk} = fftt.^2 / 2;
            f_psd = f_fft;
        end
        
        % Actualiza la posicion
        kk = kk + 1;
        
    end % for k
    
end % for win

%% Calcula el promedio y la desviacion estandar de los fft
fftmean = zeros(1, length(f_fft));
fftstd = zeros(1, length(f_fft));
fftdata = zeros(1, ng);
fftmax = zeros(1, ng);
for i = 1:length(f_fft)
    for j = 1:ng % Recorre cada grado de libertad
        fftdata(j) = fft{j}(i);
    end % for j
    fftmean(i) = mean(fftdata);
    fftstd(i) = std(fftdata);
    fftmax(i) = max(fftdata);
end % for i

%% Calcula el promedio y la desviacion estandar de los psd
psdmean = zeros(1, length(f_psd));
psdstd = zeros(1, length(f_psd));
psddata = zeros(1, ng);
psdmax = zeros(1, ng);
for i = 1:length(f_psd)
    for j = 1:ng % Recorre cada grado de libertad
        psddata(j) = psd{j}(i);
    end % for j
    psdmean(i) = mean(psddata);
    psdstd(i) = std(psddata);
    psdmax(i) = max(psddata);
end % for i

%% Calcula los peaks
maxlocs = 0;
if ~r.peakUseMean
    locs = cell(1, ng);
    for i = 1:length(locs)
        if r.peakFFT
            [~, ploc] = findpeaks(fft{i}, f_fft, ...
                'MinPeakDistance', r.peakMinDistance, 'MinPeakHeight', r.peakThreshold);
        else
            [~, ploc] = findpeaks(psd{i}, f_psd, ...
                'MinPeakDistance', r.peakMinDistance, 'MinPeakHeight', r.peakThreshold);
        end
        locs{i} = ploc;
    end % for i
else
    locs = cell(1, 1);
    if r.peakFFT
        [~, ploc] = findpeaks(fftmean, f_fft, ...
            'MinPeakDistance', r.peakMinDistance, 'MinPeakHeight', r.peakThreshold);
    else
        [~, ploc] = findpeaks(psdmean, f_psd, ...
            'MinPeakDistance', r.peakMinDistance, 'MinPeakHeight', r.peakThreshold);
    end
    locs{1} = ploc; % En las frecuencias
end

% Independiente de la frecuencia
for i = 1:length(locs)
    jlim = 0;
    for j = 1:length(locs{i})
        if locs{i}(j) >= r.peakFreqThreshold
            jlim = j;
            break
        end
    end
    locs{i} = locs{i}(jlim:end); % Corta el vector
    maxlocs = max(length(locs{i}), maxlocs);
end

%% Calcula el promedio y la desviacion estandar de las frecuencias
locMean = zeros(1, maxlocs);
locStd = zeros(1, maxlocs);

% Calcula datos pero en periodos
tlocMean = zeros(1, maxlocs);
tlocStd = zeros(1, maxlocs);
for i = 1:maxlocs
    locData = []; % Datos para la posicion i
    tlocData = [];
    for k = 1:length(locs) % Recorre cada nodo de analisis
        if i <= length(locs{k})
            locData = [locData, locs{k}(i)];
            tlocData = [tlocData, 1 / locs{k}(i)];
        end
    end % for k
    
    locMean(i) = mean(locData);
    locStd(i) = std(locData);
    
    % Estadistica para los periodos
    tlocMean(i) = mean(tlocData);
    tlocStd(i) = std(tlocData);
    
end % for i

%% Busca las posiciones de la frecuencia para locMean
locFreq_fft = []; % Frecuencias (posicion)
locFreq_psd = [];

j = 1;
for i = 1:length(f_fft)
    if f_fft(i) >= locMean(j)
        locFreq_fft(j) = i; %#ok<*AGROW>
        j = j + 1; % Avanza
        if j > maxlocs
            break;
        end
    end
end % for i

j = 1;
for i = 1:length(f_psd)
    if f_psd(i) >= locMean(j)
        locFreq_psd(j) = i;
        j = j + 1; % Avanza
        if j > maxlocs
            break;
        end
    end
end % for i

maxlocs = min(length(locFreq_fft), length(locFreq_psd));

% Peaks periodos
if ~r.peakUseMean || r.forcepeakmax
    if r.peakFFT
        pks = fftmax(locFreq_fft);
    else
        pks = psdmax(locFreq_psd);
    end
else
    if r.peakFFT
        pks = fftmean(locFreq_fft);
    else
        pks = psdmean(locFreq_psd);
    end
end

%% Calcula los amortiguamientos por cada periodo de cada nodo registrado
if r.betaFFT
    pksBeta = fftmax(locFreq_fft);
else
    pksBeta = psdmax(locFreq_psd);
end

betaNodo = cell(1, ng);
betaFreqNodo = cell(1, ng);

% Recorre cada registro
for k = 1:ng
    
    if ~r.betaFFTMax % Si se usan todos los registros
        if r.betaFFT
            ftNodo = fft{k};
            fbeta = f_fft;
        else
            ftNodo = psd{k};
            fbeta = f_psd;
        end
        pksNodo = ftNodo(locFreq_fft);
        pksObj = pksNodo ./ sqrt(2);
    else % Si se usa solo el maximo
        if r.betaFFT
            ftNodo = fftmax;
            fbeta = f_fft;
        else
            ftNodo = psdmax;
            fbeta = f_psd;
        end
        pksObj = pksBeta ./ sqrt(2);
    end
    beta = zeros(1, maxlocs);
    betaFreq = cell(1, maxlocs);
    lastj = 1;
    
    % Recorre cada peak del nodo registrado
    for i = 1:maxlocs
        for j = lastj:length(fbeta) - 1 % Al comenzar desde el punto anterior asegura que no se repitan frecuencias
            if (ftNodo(j) - pksObj(i)) * (ftNodo(j+1) - pksObj(i)) < 0 % Cruzo el objetivo en i
                % Si el ultimo que cruzo fue superior a la frecuencia del peak
                % objetivo entonces este corresponde a la frecuencia derecha, y
                % el anterior a la izquierda
                if j > locFreq_fft(i)
                    izq = fbeta(lastj);
                    der = fbeta(j);
                    lastj = j;
                    if r.betaFFT
                        beta(i) = (der - izq) / (der + izq);
                    else
                        beta(i) = (der^2 - izq^2) / (der^2 + izq^2);
                    end
                    betaFreq{i} = [izq, der, fbeta(locFreq_fft(i)), pksObj(i)];
                    break;
                end
                lastj = j + 1; % Ultimo en atravesar
            end
        end % for j
    end % for i
    
    % Guarda el resultado
    betaNodo{k} = beta;
    betaFreqNodo{k} = betaFreq;
    
    % Termina la ejecucion (k==1)
    if r.betaFFTMax
        betaNodo = betaNodo{1};
        betaFreqNodo = betaFreqNodo{1};
        break;
    end
    
end % for k

%% Calcula la envolvente de los peaks por cada una de las formas modales, usa FFT
envFormaModal = cell(1, maxlocs);
for k = 1:maxlocs
    envModo = zeros(1, ng);
    for i = 1:ng % Recorre cada nodo
        % Obtiene la fft asociada al periodo k del registro i
        envModo(i) = fft{i}(locFreq_fft(k));
    end % for i
    envModo = envModo ./ max(envModo);
    envFormaModal{k} = envModo;
end % for k

end % PSD function