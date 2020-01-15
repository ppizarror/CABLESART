cable = 1; %#ok<IJCL>

ex = psd{cable}; %#ok<*SUSENS>

% Valores de corte
threshold = [0, 1e-6, 1e-5, 6e-6, 1e-5];
peakmin = [2, 2.5, 2, 2, 2];
fthreshold = [2, 4, 4, 5, 5]; % Frecuencias minimas
maxpeaks = [5, 5, 5, 5, 5]; % Maximos peaks [16, 6, 8, 5, 5];

ex.exec('zeroFill', 2, 'fftLim', 35, 'maxPeaks', maxpeaks(cable), ...
    'psdPlot', false, ... # Grafica el PSD
    'fftPlot', false, ... # Grafica la FFT
    'fftMean', false, ... # Grafica el promedio en vez de el resto de datos
    'peaksT', true, ... # Grafica los peaks
    'betaPlot', true, ... # Grafica los beta de cada modo
    'peakUseMean', true, ... # El cálculo de peaks usa el promedio
    'peakMinDistance', peakmin(cable), 'peakThreshold', threshold(cable), ...
    'peakFreqThreshold', fthreshold(cable), ... # peaks
    'pwelch', true, 'tiempoVentanas', 20, 'overlap', 0.5, 'factor', 100, ... # Ajuste de las ventanas
    'fase', false, 'faseNodos', [1, 2], 'faseTLim', [0, 25], ... # Análisis de la fase
    'plotlog', true);