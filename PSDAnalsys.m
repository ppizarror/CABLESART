classdef PSDAnalsys < handle
    %PSDANALSYS Crea un analisis PSD de un cable
    
    properties
        cable % Indica que cable se analiza
        cable_prev % Datos previo al filtro
        etiqueta % Nombre del analisis
        fs % Factor de muestro de la señal
        dt % Tiempo de muestreo
        tiempo % Vector de tiempo
        lengthreg % Largo de los vectores de la señal
        numchannels % Numero de canales
        channelnames % Nombres de los registros
        channelcolors % Colores de cada canal
        Tn % Periodos reales de la estructura
        
        % Variables de estado
        prev_plot % Grafico previo
        
        % Ventanas
        window_fix
        window_color
        window_time
        
        % Resultados de analisis
        r_freqs % Frecuencias de cada modo
        r_beta % Amortiguamiento de cada modo
        r_nlfit % Ultimo fit
        r_freqnlfit % Frecuencias promedio de cada modo del ajuste no lineal
        
    end
    
    methods(Access = public)
        function obj = PSDAnalsys(cable_reg, etiqueta, filter_hampel)
            %PSDANALSYS Constructor de clase
            % cable_reg     registro de los distintos ensayos
            % etiqueta      nombre del ensayo
            % filter_hampel Numero de datos para realizar el filtro
            
            obj.fs = 200;
            
            segdel = 0; % Segundos a borrar primeros
            obj.cable = cable_reg(max(fix(segdel*obj.fs), 1):end, :);
            [~, n] = size(obj.cable);
            obj.numchannels = n;
            
            obj.channelnames = cell(1, obj.numchannels);
            for i = 1:obj.numchannels
                obj.channelnames{i} = sprintf('Canal %d', i);
            end
            obj.channelcolors = {}; % Almacena los colores
            
            obj.etiqueta = etiqueta;
            obj.dt = 1 / obj.fs;
            obj.lengthreg = length(obj.cable(:, 1));
            obj.tiempo = linspace(0, obj.lengthreg*obj.dt, obj.lengthreg);
            obj.Tn = ones(1, 10000);
            
            % Aplico filtros
            obj.apply_filter(filter_hampel);
            
            obj.r_freqs = [];
            obj.r_beta = [];
            obj.r_nlfit = {};
            
            % Pide las ventanas
            obj.close_prev_plot();
            obj.request_windows();
        end
        
        function f = get_mode_freq(obj)
            % Retorna las frecuencias modales del analisis del PSD
            f = obj.r_freqs;
        end
        
        function f = get_mode_freq_nlfit(obj)
            % Retorna las frecuencias modales del analisis del ajuste no
            % lineal
            f = obj.r_freqnlfit;
        end
        
        function b = get_mode_beta(obj)
            % Retorna los beta modales del analisis
            b = obj.r_beta;
        end
        
        function spectrogram_win(obj, window, channel, varargin)
            % Calcula el espectrograma de una ventana y un registro
            
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'ventana', 1);
            addOptional(p, 'factor', 2);
            addOptional(p, 'noverlap', 0.5);
            parse(p, varargin{:});
            r = p.Results;
            
            % Obtiene la ventana de datos
            wf = obj.window_fix{window};
            reg = obj.cable;
            % t = obj.tiempo(wf(1):wf(2));
            x = reg(:, channel);
            x = x(wf(1):wf(2));
            
            % Genera el espectrograma
            n_window = r.ventana * obj.fs; % Tamaño de ventana de analisis
            nfft_psd = 2^nextpow2(floor(n_window*r.factor)); % Numero de puntos de la FFT
            noverlap = floor(r.noverlap*n_window);
            
            ctitle = sprintf('Espectrograma ventana %d canal %d', window, channel);
            fig_title = sprintf('%s - %s', ...
                ctitle, obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            spectrogram(x, nfft_psd, noverlap, nfft_psd, obj.fs, 'yaxis');
            xlabel('Tiempo (min)');
            ylabel('Frecuencia (Hz)');
            title(fig_title);
            
        end % spectrogram_win function
        
        function spectrogram(obj, channel, varargin)
            % Calcula el espectrograma de una ventana y un registro
            
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'ventana', 1);
            addOptional(p, 'factor', 2);
            addOptional(p, 'noverlap', 0.9);
            parse(p, varargin{:});
            r = p.Results;
            
            % Obtiene la ventana de datos
            reg = obj.cable;
            % t = obj.tiempo(wf(1):wf(2));
            x = reg(:, channel);
            
            % Genera el espectrograma
            n_window = r.ventana * obj.fs; % Tamaño de ventana de analisis
            nfft_psd = 2^nextpow2(floor(n_window*r.factor)); % Numero de puntos de la FFT
            noverlap = floor(r.noverlap*n_window);
            
            ctitle = 'Espectrograma';
            fig_title = sprintf('%s - %s', ...
                ctitle, obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            spectrogram(x, nfft_psd, noverlap, nfft_psd, obj.fs, 'yaxis');
            xlabel('Tiempo (min)');
            ylabel('Frecuencia (Hz)');
            title(fig_title);
            
        end % spectrogram function
        
        function plot_diff(obj)
            % Grafica la diferencia entre el vector original y el
            % filtrado
            
            ctitle = 'Diferencia despues del filtrado';
            fig_title = sprintf('%s - %s', ...
                ctitle, obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            hold on;
            
            title(fig_title);
            for i = 1:obj.numchannels
                plot(obj.tiempo, obj.cable(:, i)-obj.cable_prev(:, i), 'color', obj.channelcolors{i});
            end
            legend(obj.channelnames, 'location', 'southeast');
            
            xlim([0, max(obj.tiempo)]);
            ylabel('Diferencia (-)');
            xlabel('Tiempo (s)');
            grid on;
            grid minor;
            obj.prev_plot = plt;
            
        end % plot_diff function
        
        function w = get_number_windows(obj)
            % Retorna el numero de ventanas del analisis
            w = length(obj.window_fix);
        end % get_number_windows function
        
        function out = get_reg_window(obj, window, original)
            % Retorna los registros de una ventana, el primer elemento es
            % el tiempo, el resto los canales
            
            if ~exist('original', 'var')
                original = true;
            end
            
            lw = length(obj.window_fix);
            if window > lw
                error('El numero de ventana %d no existe, max %d', window, lw);
            end
            wf = obj.window_fix{window};
            
            reg = obj.cable;
            if ~original
                reg = obj.cable_prev;
            end
            
            out = cell(1, obj.numchannels);
            t = obj.tiempo(wf(1):wf(2));
            
            out{1} = t;
            
            for i = 2:obj.numchannels
                data = reg(:, i);
                out{i} = data(wf(1):wf(2));
            end
            
        end % get_reg_window function
        
        function plot_reg_window(obj, window, original)
            % Grafica todos los registros de una ventana
            
            if ~exist('original', 'var')
                original = true;
            end
            
            lw = length(obj.window_fix);
            if window > lw
                error('El numero de ventana %d no existe, max %d', window, lw);
            end
            wf = obj.window_fix{window};
            
            ctitle = sprintf('Registros ventana %d', window);
            fig_title = sprintf('%s - %s', ...
                ctitle, obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            
            title(fig_title);
            hold on;
            
            % Grafica
            reg = obj.cable;
            if ~original
                reg = obj.cable_prev;
            end
            
            t = obj.tiempo(wf(1):wf(2));
            for i = 2:obj.numchannels
                data = reg(:, i);
                plot(t, data(wf(1):wf(2)), 'color', obj.channelcolors{i});
            end
            
            legend(obj.channelnames{2:obj.numchannels}, 'location', 'southeast');
            xlim([min(t), max(t)]);
            grid on;
            grid minor;
            ylabel('Registro (cm)');
            xlabel('Tiempo (s)');
            
            % Ajusta el eje
            ymax = max(abs(get(gca, 'ylim')));
            ylim([-ymax, ymax]);
            obj.prev_plot = plt;
            drawnow;
            
        end % plot_reg_window function
        
        function plot_reg_all(obj, original, fullscreen)
            % Grafica todos los registros
            
            if ~exist('original', 'var')
                original = true;
            end
            
            if ~exist('fullscreen', 'var')
                fullscreen = false;
            end
            
            ctitle = 'Registros';
            fig_title = sprintf('%s - %s', ...
                ctitle, obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            
            if ~fullscreen
                movegui(plt, 'center');
            else
                plt.WindowState = 'maximized';
            end
            
            title(fig_title);
            hold on;
            
            savecolors = isempty(obj.channelcolors);
            if savecolors
                obj.channelcolors{1} = [0, 0, 0];
            end
            
            % Grafica
            reg = obj.cable;
            if ~original
                reg = obj.cable_prev;
            end
            for i = 2:obj.numchannels
                if savecolors
                    pl = plot(obj.tiempo, reg(:, i));
                    obj.channelcolors{i} = get(pl, 'Color');
                else
                    plot(obj.tiempo, reg(:, i), 'color', obj.channelcolors{i});
                end
            end
            
            legend(obj.channelnames{2:obj.numchannels}, 'location', 'southeast');
            xlim([0, max(obj.tiempo)]);
            grid on;
            grid minor;
            ylabel('Registro (cm)');
            xlabel('Tiempo (s)');
            
            % Ajusta el eje
            ymax = max(abs(get(gca, 'ylim')));
            ylim([-ymax, ymax]);
            obj.prev_plot = plt;
            drawnow;
            
        end % plot_reg_all function
        
        function plot_reg_all_windows(obj, original, fullscreen)
            % Mismo que plot_reg_all pero grafica las ventanas
            
            if ~exist('original', 'var') % original o prev
                original = true;
            end
            
            if ~exist('fullscreen', 'var')
                fullscreen = false;
            end
            
            obj.plot_reg_all(original, fullscreen);
            
            % Dibuja las lineas de las ventanas
            for i = 1:length(obj.window_color)
                for j = 1:length(obj.window_time{i})
                    pl = drawVxLine(obj.window_time{i}(j), '--', 1, obj.window_color{i});
                    set(get(get(pl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                end
            end
            
        end % plot_reg_all_windows function
        
        function plot_reg(obj, reg)
            % Grafica el registro
            
            if reg == 1
                ctitle = 'Registro Aceleración';
                ytitle = 'cm/s^2';
            else
                ctitle = 'Registro Desplazamiento';
                ytitle = 'cm';
            end
            fig_title = sprintf('%s - %s', ...
                ctitle, obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            hold on;
            plot(obj.tiempo, obj.cable(:, reg), 'color', obj.channelcolors{reg});
            grid on;
            grid minor;
            title(fig_title);
            ylabel(ytitle);
            xlabel('Tiempo (s)');
            obj.prev_plot = plt;
            
        end % plot_reg function
        
        function calc_nl(obj, varargin)
            % calc_nl: calcula identificacion no lineal
            %
            % Parametros opcionales:
            %   betalim                 Limite inferior/superior amortiguamiento (2x1)
            %   betaRayleigh            Los amortiguamientos los calcula con Rayleigh
            %   functionTolerance       Tolerancia maxima (lsqnonlin)
            %   maxFunctionEvaluations  Numero maximo de evaluaciones (lsqnonlin)
            %   nmodos                  Numero de modos de analisis
            %   rholim                  Limite inferior/superior amplitud modo (2x1)
            %   thetalim                Limite inferior/superior fase (2x1)
            %   unidadL                 Unidad longitud
            %   wlim                    Limite inferior/superior frecuencia (2x1)
            
            % Inicia proceso
            tinicial = clock;
            fprintf('Identificacion no lineal:\n');
            
            % Recorre parametros opcionales
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'betalim', [0, Inf]);
            addOptional(p, 'betaRayleigh', true);
            addOptional(p, 'functionTolerance', 1e-12);
            addOptional(p, 'maxFunctionEvaluations', 5000);
            addOptional(p, 'modos', []);
            addOptional(p, 'load', false);
            addOptional(p, 'rholim', [-Inf, Inf]);
            addOptional(p, 'thetalim', [-Inf, Inf]);
            addOptional(p, 'unidadL', 'cm');
            addOptional(p, 'wlim', [0, Inf]);
            parse(p, varargin{:});
            r = p.Results;
            
            % Si está vacio retorna
            if isempty(r.modos)
                r.modos = 1:length(obj.r_freqs);
            end
            if min(r.modos) < 1
                error('No pueden haber modos negativos');
            end
            
            % Verifica variables
            r.modos = floor(r.modos);
            if max(r.modos) > length(obj.r_freqs)
                error('Numero modos excede el analisis');
            end
            if length(r.rholim) ~= 2 || length(r.thetalim) ~= 2 || ...
                    length(r.betalim) ~= 2 || length(r.wlim) ~= 2
                error('Parametros opcionales deben ser vectores de dos componentes');
            end
            
            % Ordena los limites
            r.rholim = sort(r.rholim);
            r.thetalim = sort(r.thetalim);
            r.betalim = sort(r.betalim);
            r.wlim = sort(r.wlim);
            
            xf_total = {};
            dfit_total = {};
            J_total = {};
            despl_total = {};
            t_total = {};
            
            % Obtengo el amortiguamiento y las frecuencias de los modos requeridos
            beta = obj.r_beta(r.modos);
            omega = 2 * pi * obj.r_freqs(r.modos);
            
            t_win = {};
            nventana = {};
            
            if r.load
                if isempty(obj.r_nlfit)
                    fprintf('\tError al cargar resultados guardados, se reiniciara el proceso\n');
                    r.load = false;
                else % Ya existe la entrada
                    fprintf('\tCargando resultados guardados\n');
                end
            end
            
            % Obtiene el registro
            k = 0;
            ncanal = 0;
            for window = 1:obj.get_number_windows()
                data = obj.get_reg_window(window);
                
                t = data{window}';
                t = linspace(0, max(t)-min(t), length(t))';
                t_win{window} = t;
                
                % Llamamos a la funcion
                if ~r.load
                    fprintf('\tOptimizando funcion, ventana %d\n', window);
                end
                for i = 2:length(data)
                    k = k + 1;
                    despl = data{i};
                    if ~r.load
                        [~, xf, dfit, J] = NLFIT(despl, t, omega, beta, length(r.modos), r.rholim, r.thetalim, ...
                            r.betalim, r.wlim, 'maxFunctionEvaluations', r.maxFunctionEvaluations, ...
                            'functionTolerance', r.functionTolerance);
                    else
                        xf = obj.r_nlfit{k, 1}; %#ok<*USENS>
                        dfit = obj.r_nlfit{k, 2};
                        J = obj.r_nlfit{k, 3};
                    end
                    
                    % Añade una columna a xf
                    xf = [xf, zeros(length(r.modos), 1)];
                    
                    % Calcula qué modo pertenece
                    for h = 1:length(r.modos)
                        modow = 0; % Modo ganador
                        dist = Inf;
                        for tmodo = 1:length(r.modos) % Recorre cada modo
                            dd = abs(omega(tmodo)-xf(h, 1)); % Distancia
                            if dd < dist % Actualiza
                                dist = dd;
                                modow = tmodo;
                            end
                        end
                        xf(h, 5) = modow;
                    end
                    
                    xf_total{k} = xf; %#ok<*AGROW>
                    dfit_total{k} = dfit;
                    J_total{k} = J;
                    nventana{k} = window;
                    despl_total{k} = despl;
                    t_total{k} = t;
                end
                
                if ncanal == 0
                    ncanal = k;
                end
            end
            
            % Guarda los resultados
            if ~r.load
                fprintf('\tGuardando resultados\n');
                obj.r_nlfit = cell(k, 3);
                for i = 1:k
                    obj.r_nlfit{i, 1} = xf_total{i};
                    obj.r_nlfit{i, 2} = dfit_total{i};
                    obj.r_nlfit{i, 3} = J_total{i};
                end
            end
            
            % Tabula los valores iniciales y finales de las iteraciones
            % fprintf('\tValores iniciales de la optimizacion:\n');
            % obj.tabularAnalisisIdentificacionNL(xo);
            
            % fprintf('\tValores finales de la optimizacion:\n');
            % obj.tabularAnalisisIdentificacionNL(xf);
            
            % =============================================================
            % Realiza las interpolaciones
            % =============================================================
            
            % Obtiene el t mayor
            indx = 0;
            lmax = 0;
            for window = 1:obj.get_number_windows()
                if length(t_win{window}) > lmax
                    indx = window;
                    lmax = length(t_win{window});
                end
            end
            t = t_win{indx};
            
            % Interpola el resto de datos de J, despl y dfit para
            J_interp = {};
            % despl_interp = {};
            % dfit_interp = {};
            
            for i = 1:length(J_total)
                J_interp{i} = J_total{i};
                J_interp{i}(end +1:length(t)) = 0;
                % despl_interp{i} = interp1(t_total{i}, despl_total{i}, t, 'linear', 'extrap');
                % dfit_interp{i} = interp1(t_total{i}, dfit_total{i}, t, 'linear', 'extrap');
            end
            
            % Calcula los std y las desviaciones estandar
            J_mean = zeros(1, length(t));
            J_std = zeros(1, length(t));
            data = zeros(1, length(J_interp));
            
            for i = 1:length(t)
                for j = 1:length(J_interp)
                    data(j) = J_interp{j}(i);
                end
                J_mean(i) = mean(data);
                J_std(i) = std(data);
            end
            J_mean = medfilt1(J_mean, 500);
            J_std = medfilt1(J_std, 500);
            
            % =============================================================
            % Analisis frecuencias por modo
            % =============================================================
            fig_title = sprintf('Análisis de frecuencias - %s', obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            hold on;
            ylim([0, ncanal]);
            
            % Grafica las frecuencias del PSD
            legplot = cell(1, length(omega));
            for i = 1:length(omega)
                drawVxLine(omega(i)/(2 * pi), '-', 2.5);
                legplot{i} = sprintf('Modo %d PSD', i);
            end
            
            marcador = {'+', '*', 'x', 'o', '.'};
            cn = 1; % Numero de canal
            
            freqs = cell(length(r.modos), 1);
            for i = 1:length(r.modos)
                freqs{i} = [];
            end
            
            for j = 1:k
                xf = xf_total{j};
                [n, ~] = size(xf);
                for w = 1:n
                    m = xf(w, 5); % Modo del resultado
                    f = xf(w, 1) / (2 * pi);
                    freqs{m} = [freqs{m}, f];
                    pl = plot(f, cn, marcador{m}, 'markersize', 7);
                    set(get(get(pl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                end
                cn = cn + 1; % Aumenta el canal
                if cn > ncanal
                    cn = 1;
                end
            end
            
            legend(legplot, 'location', 'northeast', 'fontsize', 7);
            title({fig_title, ''});
            xlabel('Frecuencia (s)');
            ylabel('Canal');
            grid on;
            grid minor;
            yTickInteger();
            
            % =============================================================
            % Realiza el histograma de cada frecuencia
            % =============================================================
            fig_title = sprintf('Análisis de frecuencias - %s', obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            
            hold on;
            freqtab = zeros(length(r.modos), 4);
            obj.r_freqnlfit = zeros(length(r.modos), 1);
            for m = 1:length(r.modos) % Recorre cada modo
                subplot(length(r.modos), 1, m);
                if m == 1
                    title({fig_title, ''});
                end
                hold on;
                hist(freqs{m});
                grid on;
                grid minor;
                f_psd = omega(m) / (2 * pi);
                f_prom = mean(freqs{m});
                f_median = median(freqs{m});
                
                freqtab(m, 1) = f_psd;
                freqtab(m, 2) = f_prom;
                freqtab(m, 3) = std(freqs{m});
                freqtab(m, 4) = f_median;
                obj.r_freqnlfit(m) = f_prom;
                
                drawVxLine(f_psd, '-', 4);
                drawVxLine(f_prom, '--', 1.5);
                drawVxLine(f_median, '--', 1.5);
                if m == 1
                    legend({'Histograma', 'PSD', 'Promedio', 'Mediana'}, 'location', 'best');
                end
                ylabel(sprintf('Modo %d', m));
            end
            xlabel('Frecuencia (Hz)');
            
            fprintf('\tAnálisis de frecuencias:\n');
            obj.tabularNLFreq(freqtab);
            
            % =============================================================
            % Análisis de amortiguamientos
            % =============================================================
            fig_title = sprintf('Análisis de amortiguamientos - %s', obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            hold on;
            
            % Betas por cada modo
            betas = cell(length(r.modos), 1);
            for i = 1:length(r.modos)
                betas{i} = [];
            end
            
            for j = 1:k
                xf = xf_total{j};
                [n, ~] = size(xf);
                for w = 1:n
                    m = xf(w, 5); % Modo del resultado
                    b = xf(w, 2); % beta
                    betas{m} = [betas{m}, b];
                    pl = plot(m, b, marcador{m}, 'markersize', 10);
                    set(get(get(pl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                end
            end
            
            for i = 1:length(omega)
                % drawVxLine(i, '--', 1, [0, 0, 0]);
            end
            
            grid on;
            grid minor;
            xlim([0, length(r.modos) + 1]);
            xTickInteger();
            title(fig_title);
            xlabel('Modo');
            ylabel('Amortiguamiento (-)');
            
            % =============================================================
            % Realiza el histograma de cada amortiguamiento
            % =============================================================
            fig_title = sprintf('Análisis de amortiguamientos - %s', obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            
            hold on;
            betatab = zeros(length(r.modos), 4);
            for m = 1:length(r.modos) % Recorre cada modo
                subplot(length(r.modos), 1, m);
                if m == 1
                    title({fig_title, ''});
                end
                hold on;
                hist(betas{m});
                grid on;
                grid minor;
                b_psd = beta(m);
                b_prom = mean(betas{m});
                b_median = median(betas{m});
                
                betatab(m, 1) = b_psd;
                betatab(m, 2) = b_prom;
                betatab(m, 3) = std(betas{m});
                betatab(m, 4) = b_median;
                
                % drawVxLine(b_psd, '-', 4);
                drawVxLine(b_prom, '-', 3);
                drawVxLine(b_median, '-', 3);
                if m == 1
                    legend({'Histograma', 'Promedio', 'Mediana'}, 'location', 'northeast');
                end
                ylabel(sprintf('Modo %d', m));
            end
            xlabel('Amortiguamiento (-)');
            
            fprintf('\tAnálisis de amortiguamientos:\n');
            obj.tabularNLAmort(betatab);
            
            % AJUSTE
            fig_title = sprintf('Respuesta de Desplazamiento Ajustada - %s', obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            hold on;
            for i = 1:length(despl_total)
                color = [0, 0, 0.5 * (1 + rand())];
                if i == 1
                    color = [0, 0, 1];
                end
                if i ~= 1
                    pl = plot(t_total{i}, despl_total{i}, 'color', color);
                    set(get(get(pl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                else
                    plot(t_total{i}, dfit_total{i}, 'color', color, 'linewidth', 2);
                end
            end
            for i = 1:length(dfit_total)
                color = [0.5 * (1 + rand()), 0, 0];
                if i == 1
                    color = [1, 0, 0];
                end
                if i ~= 1
                    pl = plot(t_total{i}, dfit_total{i}, 'color', color);
                    set(get(get(pl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                else
                    plot(t_total{i}, dfit_total{i}, 'color', color, 'linewidth', 2);
                end
            end
            title({fig_title, ''});
            xlabel('Tiempo (s)');
            ylabel(sprintf('Desplazamiento (%s)', r.unidadL));
            grid on;
            grid minor;
            xlim([min(t), max(t)]);
            ymax = max(abs(get(gca, 'ylim')));
            ylim([-ymax, ymax]);
            legend({'Respuesta real', 'Respuesta ajustada'}, 'location', 'northeast');
            
            % J
            fig_title = sprintf('Historial de Error - %s', obj.etiqueta);
            plt = figure('Name', fig_title, 'NumberTitle', 'off');
            movegui(plt, 'center');
            hold on;
            for i = 1:length(J_interp)
                plot(t, J_interp{i});
            end
            plot(t, J_mean, 'k-', 'linewidth', 2);
            plot(t, J_mean+3*J_std, 'k--', 'linewidth', 2);
            % plot(t, J_mean-3*J_std, 'k--', 'linewidth', 2);
            title({fig_title, ''});
            xlabel('Tiempo (s)');
            ylabel(sprintf('Función de error J (%s)', r.unidadL));
            grid on;
            grid minor;
            % ylim([0, max(get(cga, 'YTick'))]);
            
            % Finaliza proceso
            fprintf('\tProceso finalizado en %.2f segundos\n', etime(clock, tinicial));
            
        end % calcularIdentificacionNL function
        
        function exec(obj, varargin)
            % EXEC Ejecuta el analisis
            %
            % Parametros opcionales:
            %   betaFFT         El amortiguamiento se calcula con FFT en vez de PSD
            %   betaLineWidth   Ancho de linea de los graficos de amortiguamiento
            %   betaPlot        Grafica el calculo del amortiguamiento de cada modo
            %   betaPlotComp    Grafico amortiguamiento modal comparado con el real
            %   betaRayleigh    Los amortiguamientos los calcula con Rayleigh
            %   closeAll        Cierra todos los graficos antes del analisis
            %   fase            Realiza analisis de fases
            %   faseNodos       Nodos en los que se realiza la fase
            %   faseTlim        Limites periodo grafico fase
            %   fftLim          Limite de frecuencia en grafico FFT
            %   fftMeanStd      Grafica el promedio y desviacion estandar para FFT
            %   fftPlot         Muestra el grafico de la FFT simple
            %   filtMod         Realiza analisis de filtros
            %   filtNodo        Nodos de analisis de filtros
            %   filtRange       Rango de cada peak considerado en el analisis del filtro
            %   filtTlim        Limite periodo grafico filtros
            %   formaModal      Vector con periodos a graficar de las formas modales
            %   formaModalComp  Comparacion formas modales con las teoricas
            %   formaModalDir   Vector direccion de analisis formas modales (x,y,z)
            %   formaModalError Grafico error forma modal con la teorica
            %   formaModalLeg   Muestra la leyenda de las formas modales
            %   formaModalLw    Ancho de linea formas modales
            %   formaModalMark  Muestra un marcador en cada nodo de las formas modales
            %   formaModalMz    Tamano marcador nodo en formas modales
            %   formaModalPlot  Grafico de las formas modales
            %   legend          Muestra la leyenda
            %   legendloc       Ubicacion de la leyenda
            %   maxPeaks        Numero de peaks maximos calculados
            %   peakMinDistance Distancia minima entre peaks requerida
            %   peaksFFT        El calculo de peaks de periodos es con FFT en vez de PSD
            %   peaksT          Grafica los peaks
            %   peaksTComp      Grafico comparacion peaks teorico y FFT
            %   peaksTError     Error peaks periodos con respecto al teorico
            %   psdMeanStd      Grafica el promedio y desviacion estandar para PSD
            %   psdPlot         Grafica el PSD por cada frecuencia
            %   tmax            Tiempo maximo de analisis
            %   tmin            Tiempo minimo de analisis
            %   tukeywinr       Factor de la ventana de tukey
            %   unidadL         Unidad longitud
            %   zeroFill        Indica relleno de ceros para FFT
            
            % Recorre parametros opcionales
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'betaFFT', true);
            addOptional(p, 'betaLineWidth', 1.75);
            addOptional(p, 'betaPlot', false);
            addOptional(p, 'betaPlotComp', false);
            addOptional(p, 'betaRayleigh', true);
            addOptional(p, 'closeAll', false);
            addOptional(p, 'fase', false);
            addOptional(p, 'faseNodos', []);
            addOptional(p, 'faseTlim', [0, 1]);
            addOptional(p, 'fftLim', 0);
            addOptional(p, 'fftMean', false);
            addOptional(p, 'fftStd', false);
            addOptional(p, 'fftPlot', false);
            addOptional(p, 'filtMod', []);
            addOptional(p, 'filtNodo', {});
            addOptional(p, 'filtRange', 0.2);
            addOptional(p, 'filtTlim', [0, 1]);
            addOptional(p, 'formaModal', []);
            addOptional(p, 'formaModalColorFactor', 2.5);
            addOptional(p, 'formaModalComp', false);
            addOptional(p, 'formaModalDir', [0, 0, 0]); % Puede ser [0, 1, 0] (y)
            addOptional(p, 'formaModalError', false);
            addOptional(p, 'formaModalErrorLegPos', 'best');
            addOptional(p, 'formaModalErrorStdLegPos', 'best');
            addOptional(p, 'formaModalLeg', true);
            addOptional(p, 'formaModalLegPos', 'best');
            addOptional(p, 'formaModalLw', 1.5);
            addOptional(p, 'formaModalMark', true);
            addOptional(p, 'formaModalMz', 5);
            addOptional(p, 'formaModalPlot', true);
            addOptional(p, 'legend', false);
            addOptional(p, 'fmin', 0.5);
            addOptional(p, 'legendloc', 'best');
            addOptional(p, 'maxPeaks', -1); % El maximo se ajusta al dato
            addOptional(p, 'peakThreshold', 0); % Limite
            addOptional(p, 'peakFreqThreshold', 0);
            addOptional(p, 'peakUseMean', false);
            addOptional(p, 'plotlog', false);
            addOptional(p, 'peakMinDistance', 0.02); % Requerido para el calculo
            addOptional(p, 'peaksFFT', false); % El calculo de peaks es con FFT o PSD
            addOptional(p, 'peaksT', false);
            addOptional(p, 'peaksTComp', false);
            addOptional(p, 'peaksTCompErrorBar', true); % Nucleo, no se busca que se use en produccion
            addOptional(p, 'peaksTError', false);
            addOptional(p, 'psdMeanStd', false);
            addOptional(p, 'psdPlot', false);
            addOptional(p, 'tukeywinr', 0.05);
            addOptional(p, 'unidadL', 'm');
            addOptional(p, 'zeroFill', 0);
            
            % Adicionales
            addOptional(p, 'pwelch', true);
            addOptional(p, 'tiempoVentanas', 10);
            addOptional(p, 'factor', 1);
            addOptional(p, 'overlap', 0.5);
            parse(p, varargin{:});
            r = p.Results;
            
            if r.plotlog
                plotlog = @(x) log(x);
            else
                plotlog = @(x) x;
            end
            
            [f, f_psd, psd, fft, fftcomp, ~, tlocMean, tlocStd, locMean, ...
                ~, ~, maxlocs, pks, beta, betaFreq, fftmean, fftstd, ...
                psdmean, psdstd] = PSD(obj.cable, 200, 2:obj.numchannels, ...
                'peakMinDistance', r.peakMinDistance, ...
                'tukeywinr', r.tukeywinr, ...
                'zeroFill', r.zeroFill, ...
                'betaFFTMax', true, ...
                'peakFFT', r.peaksFFT, ...
                'betaFFT', r.betaFFT, ...
                'peakUseMean', r.peakUseMean, ...
                'peakThreshold', r.peakThreshold, ...
                'fmin', r.fmin, ...
                'forcepeakmax', ~r.fftMean, ...
                'usewindow', true, ...
                'windows', obj.window_fix, ...
                'usepwelch', r.pwelch, ...
                'ventana', r.tiempoVentanas, ...
                'factor', r.factor, ...
                'noverlap', r.overlap, ...
                'peakFreqThreshold', r.peakFreqThreshold);
            ctitle = '';
            ng = length(fft);
            obj.r_beta = beta;
            obj.r_freqs = locMean;
            
            % Grafica la fft de cada nodo
            if r.fftPlot
                
                fig_title = sprintf('%s %s - Analisis FFT', ...
                    ctitle, obj.etiqueta);
                plt = figure('Name', fig_title, 'NumberTitle', 'off');
                movegui(plt, 'center');
                hold on;
                for i = 1:ng
                    plot(f, plotlog(fft{i}), '-');
                end % for i
                
                % Grafica el promedio y desviacion estandar
                if r.fftMean
                    plot(f, plotlog(fftmean), 'k-', 'lineWidth', 2);
                end
                if r.fftStd
                    plot(f, plotlog(fftmean+fftstd), 'k--', 'lineWidth', 1);
                    plot(f, plotlog(fftmean-fftstd), 'k--', 'lineWidth', 1);
                end
                
                ylabel('FFT');
                if r.plotlog
                    ylabel('log (FFT)');
                end
                xlabel('Frecuencia (Hz)');
                title({fig_title, ''});
                if r.fftLim == 0
                    xlim([min(f), max(f)]);
                else
                    xlim([min(f), r.fftLim]);
                end
                ylim([0, max(get(gca, 'ylim'))]);
                grid on;
                if r.legend
                    legend(legnodos, 'location', r.legendloc);
                end
                obj.prev_plot = plt;
                
            end % fftPlot
            
            % Grafica el PSD por cada nodo
            if r.psdPlot
                
                fig_title = sprintf('%s %s - Analisis PSD', ...
                    ctitle, obj.etiqueta);
                plt = figure('Name', fig_title, 'NumberTitle', 'off');
                movegui(plt, 'center');
                hold on;
                for i = 1:ng
                    plot(f_psd, plotlog(psd{i}), '-');
                end % for i
                
                % Grafica el promedio y desviacion estandar
                if r.psdMeanStd
                    plot(f_psd, plotlog(psdmean), 'k-', 'lineWidth', 3);
                    plot(f_psd, plotlog(psdmean+psdstd), 'k--', 'lineWidth', 1.5);
                    plot(f_psd, plotlog(psdmean-psdstd), 'k--', 'lineWidth', 1.5);
                end
                
                ylabel('PSD');
                if r.plotlog
                    ylabel('log (PSD)');
                end
                xlabel('Frecuencia (Hz)');
                title({fig_title, ''});
                if r.fftLim == 0
                    xlim([min(f_psd), max(f_psd)]);
                else
                    xlim([min(f_psd), r.fftLim]);
                end
                if ~r.plotlog
                    ylim([0, max(get(gca, 'ylim'))]);
                end
                grid on;
                if r.legend
                    legend(legnodos, 'location', r.legendloc);
                end
                obj.prev_plot = plt;
                
            end % psdPlot
            
            % Tabla de periodos
            maxlocsDisp = maxlocs; % Puntos a mostrar
            if r.maxPeaks > 0
                maxlocsDisp = min(maxlocs, r.maxPeaks);
            end
            
            % Imprime en consola la tabla de los peaks
            peakMethod = 'FFT';
            if ~r.peaksFFT
                peakMethod = 'PSD';
            end
            
            peakMethodTitle = peakMethod;
            if r.plotlog
                peakMethodTitle = sprintf('log (%s)', peakMethodTitle);
            end
            if r.peaksFFT
                fplot = f;
            else
                fplot = f_psd;
            end
            
            fprintf('\tAnalisis de peaks (%s), periodos formas modales:\n', peakMethod);
            fprintf('\t\t|\tN\t|\tT peak\t\t\t|\tF peak\t|\tValor\t\t|\n');
            fprintf('\t\t---------------------------------------------------------\n');
            for i = 1:maxlocsDisp
                fprintf('\t\t|\t%d\t|\t%.2f +- %.2f\t|\t%.2f\t|\t%.2e\t|\n', ...
                    i, tlocMean(i), tlocStd(i), locMean(i), pks(i));
            end % for i
            fprintf('\t\t---------------------------------------------------------\n');
            
            % Grafica los peaks
            if r.peaksT
                
                % Grafico de peaks
                fig_title = sprintf('%s %s - Analisis %s peaks', ...
                    ctitle, obj.etiqueta, peakMethod);
                plt = figure('Name', fig_title, 'NumberTitle', 'off');
                movegui(plt, 'center');
                hold on;
                
                if r.peakUseMean && r.fftMean
                    if r.peaksFFT
                        plot(fplot, plotlog(fftmean), '-');
                    else
                        plot(fplot, plotlog(psdmean), '-');
                    end
                else
                    for i = 1:ng % Recorre cada nodo
                        if r.peaksFFT
                            plot(fplot, plotlog(fft{i}), '-');
                        else
                            plot(fplot, plotlog(psd{i}), '-');
                        end
                    end % for i
                    % legend(obj.channelnames{2:obj.numchannels}, 'location', 'northeast');
                end
                
                % Limita los peaks
                locMeanL = locMean(1:maxlocsDisp);
                pksL = pks(1:maxlocsDisp);
                text(locMeanL+0.25, plotlog(pksL.*1.0), num2str((1:numel(pksL))'));
                
                % Dibuja triangulitos
                pl = plot(locMeanL, plotlog(pksL), 'r^', 'markerfacecolor', [1, 0, 0], 'markersize', 4);
                set(get(get(pl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                ylabel(peakMethodTitle);
                xlabel('Frecuencia (Hz)');
                title({fig_title, ''});
                if r.fftLim == 0
                    xlim([min(fplot), max(fplot)]);
                else
                    xlim([min(fplot), r.fftLim]);
                end
                if ~r.plotlog
                    ylim([0, max(get(gca, 'ylim'))]);
                end
                grid on;
                grid minor;
                obj.prev_plot = plt;
                
            end % peaksT
            
            % Amortiguamientos
            betaMethod = 'FFT';
            if ~r.betaFFT
                betaMethod = 'PSD';
            end
            fprintf('\tAmortiguamiento por periodos (%s):\n', betaMethod);
            fprintf('\t\t|\tN\t|\t%%Beta\t|\n');
            fprintf('\t\t---------------------\n');
            for i = 1:maxlocsDisp
                if isempty(betaFreq{i}) % Si no se encontro el modo retorna
                    continue;
                end
                fprintf('\t\t|\t%d\t|\t%.3f\t|\n', i, beta(i));
            end % for i
            fprintf('\t\t---------------------\n');
            
            % Grafica los limites de las frecuencias
            if r.betaPlot
                
                % No aplicar plotlog
                plotlogprev = r.plotlog;
                r.plotlog = false;
                plotlog = @(x) x;
                
                % Crea la figura del amortiguamiento en cada FFT
                fig_title = sprintf('%s %s - Calculo de amortiguamientos %s', ...
                    ctitle, obj.etiqueta, betaMethod);
                plt = figure('Name', fig_title, 'NumberTitle', 'off');
                movegui(plt, 'center');
                hold on;
                
                if r.betaFFT
                    fplot = f;
                else
                    fplot = f_psd;
                end
                
                for i = 1:ng % Recorre cada nodo
                    if r.betaFFT
                        plot(fplot, plotlog(fft{i}), '-');
                    else
                        plot(fplot, plotlog(psd{i}), '-');
                    end
                end % for i
                
                for i = 1:maxlocsDisp
                    if isempty(betaFreq{i})
                        continue;
                    end
                    gc = plot([betaFreq{i}(1), betaFreq{i}(1)], plotlog([0, betaFreq{i}(4)]), ...
                        '-', 'lineWidth', r.betaLineWidth);
                    c = get(gc, 'Color');
                    plot([betaFreq{i}(2), betaFreq{i}(2)], plotlog([0, betaFreq{i}(4)]), ...
                        '-', 'lineWidth', r.betaLineWidth, 'color', c);
                    plot([betaFreq{i}(1), betaFreq{i}(2)], plotlog([betaFreq{i}(4), betaFreq{i}(4)]), ...
                        '-', 'lineWidth', r.betaLineWidth, 'color', c);
                    plot(betaFreq{i}(3), plotlog(betaFreq{i}(4)), '^', 'markerSize', 5, ...
                        'markerfacecolor', c, 'color', c);
                    text(betaFreq{i}(3).*1.05, plotlog(betaFreq{i}(4).*1.0), num2str(i));
                end % for i
                
                if r.plotlog
                    betaMethodTitle = sprintf('log (%)', betaMethod);
                else
                    betaMethodTitle = betaMethod;
                end
                ylabel(betaMethodTitle);
                xlabel('Frecuencia (Hz)');
                title({fig_title, ''});
                if r.fftLim == 0
                    xlim([min(fplot), max(fplot)]);
                else
                    xlim([min(fplot), r.fftLim]);
                end
                if ~r.plotlog
                    ylim([0, max(get(gca, 'ylim'))]);
                end
                grid on;
                grid minor;
                
                r.plotlog = plotlogprev;
                
            end % betaPlot
            
            % Grafica la fase
            if r.fase
                
                % Extrae los fft de los nodos
                faseNodos = sort(r.faseNodos);
                for k = 1:length(faseNodos)
                    if r.faseNodos(k) > length(fftcomp)
                        error('faseNodos %d excede el numero de nodos analizados por la funcion', k);
                    end
                end
                if length(r.faseNodos) ~= 2
                    error('faseNodos debe contener dos elementos');
                end
                for i = 1:length(fftcomp) - 1
                    division(:, i) = fftcomp{i} ./ fftcomp{end};
                end % for i
                
                Fftcomp(:, 1) = fftcomp{faseNodos(1)};
                Fftcomp(:, 2) = fftcomp{faseNodos(2)};
                
                % Crea la figura
                fig_title = sprintf('%s %s - Fase de la Transformada', ...
                    ctitle, obj.etiqueta);
                plt = figure('Name', fig_title, 'NumberTitle', 'off');
                movegui(plt, 'center');
                hold on;
                subplot(311);
                plot(f, abs(division(:, faseNodos(1))), 'lineWidth', 1.5);
                title('Division FFT - Parte Real');
                xlabel('Frecuencia (Hz)');
                ylabel('FFT');
                grid on;
                grid minor;
                xlim(r.faseTlim);
                
                subplot(312);
                plot(f, angle(division(:, faseNodos(1)))./pi, 'lineWidth', 1.5);
                title('Fase');
                xlabel('Frecuencia (Hz)');
                ylabel('Angulo');
                grid on;
                grid minor;
                xlim(r.faseTlim);
                
                % Modifica los ticks del grafico del angulo
                yTickInteger();
                ytickp = get(gca, 'YTick');
                if length(ytickp) == 3
                    ytick = {'-\pi', '', '\pi'};
                    set(gca, 'YTick', ytickp, 'yticklabel', ytick);
                end
                
                subplot(313);
                plot(f, abs(Fftcomp'));
                xlabel('Frecuencia (Hz)');
                ylabel('FFT');
                grid on;
                grid minor;
                zoom on;
                title('FFT Nodos');
                legend({sprintf('Nodo %s', num2str(faseNodos(1))), ...
                    sprintf('Nodo %s', num2str(faseNodos(2)))});
                xlim(r.faseTlim);
                
            end % fase
            
        end % exec function
        
    end
    
    methods(Access = private)
        
        function apply_filter(obj, filter_hampel)
            % Aplica los filtros a las señales
            
            % Guarda una copia
            obj.cable_prev = obj.cable;
            
            % Aplica hampel
            for i = 1:obj.numchannels
                if filter_hampel > 1
                    obj.cable(:, i) = hampel(obj.cable(:, i), filter_hampel, 3);
                end
            end
            
            % Elimina el primer valor
            for i = 1:obj.numchannels
                obj.cable(:, i) = obj.cable(:, i) - obj.cable(1, i);
            end
            
            % Aplica detrend
            for i = 1:obj.numchannels
                obj.cable(:, i) = detrend(obj.cable(:, i), 1);
            end
            
        end % apply_filter function
        
        function close_prev_plot(obj)
            % Cierra la última ventana
            
            if obj.prev_plot ~= -1
                close(obj.prev_plot);
                obj.prev_plot = -1;
            end
            
        end % close_prev_plot
        
        function request_windows(obj)
            % Pide las ventanas al usuario
            
            % Grafica las ventanas
            obj.plot_reg_all();
            
            % Pregunta si eliminar algun registro
            regs = inputdlg({'Eliminar canal(es) (N°)'}, 'Filtrar');
            if ~isempty(regs)
                obj.close_prev_plot();
                
                regs = regs{1};
                regs_del = strsplit(regs, ',');
                regs_int = [];
                for i = 1:length(regs_del)
                    ch = floor(str2num(regs_del{i})); %#ok<*ST2NM>
                    if ch <= 1
                        error('Canal no puede ser inferior o igual a 1');
                    end
                    if ch > obj.numchannels
                        error('Numero de canal %d no puede ser superior a la muestra', ch);
                    end
                    regs_int = [regs_int, ch];
                end
                regs_int = unique(regs_int);
                
                % Crea el nuevo vector
                j = 1;
                newchannels = obj.numchannels - length(regs_int);
                cablenew = zeros(length(obj.cable), newchannels);
                cablenew_prev = zeros(length(obj.cable), newchannels);
                newchannelnames = cell(1, newchannels);
                newchannelcolors = cell(1, newchannels);
                
                for i = 1:obj.numchannels
                    if ~ismember(i, regs_int)
                        cablenew(:, j) = obj.cable(:, i); % Copia los datos
                        cablenew_prev(:, j) = obj.cable(:, i);
                        newchannelnames{j} = obj.channelnames{i};
                        newchannelcolors{j} = obj.channelcolors{i};
                        j = j + 1;
                    end
                end
                
                % Actualiza el estado
                obj.channelnames = newchannelnames;
                obj.cable = cablenew;
                obj.cable_prev = cablenew_prev;
                obj.numchannels = newchannels;
                obj.channelcolors = newchannelcolors;
                
                % Grafica el nuevo estado
                obj.plot_reg_all();
                
            end
            
            % Pregunta cuantas ventanas se van a usar
            windows_total = inputdlg({'Cuantas ventanas se usaran?'}, 'Ventanas');
            if ~isempty(windows_total)
                obj.close_prev_plot();
                windows_total = floor(str2num(windows_total{1}));
                if windows_total < 1
                    error('El numero de ventanas debe ser mayor o igual a 1');
                end
                obj.window_fix = {}; % Almacena los puntos de cada ventana
                obj.window_time = {}; % Almacena los tiempos de cada ventana
                obj.window_color = {}; % Almacena los colores
                for i = 1:windows_total
                    obj.window_color{i} = rand(1, 3);
                    obj.window_time{i} = [];
                    
                    % Pide el primer punto
                    obj.plot_reg_all_windows(true, true);
                    title(sprintf('Ventana %d - Seleccione el primer punto', i));
                    
                    data_region = ginput(1);
                    lim1 = data_region(1, 1);
                    obj.window_time{i} = lim1;
                    obj.close_prev_plot();
                    
                    % Pide el segundo punto
                    obj.plot_reg_all_windows(true, true);
                    title(sprintf('Ventana %d - Seleccione el segundo punto', i));
                    
                    data_region = ginput(1);
                    lim2 = data_region(1, 1);
                    obj.window_time{i} = [lim1, lim2];
                    obj.close_prev_plot();
                    
                    if lim2 < lim1
                        error('Los limites de la ventana deben ser crecientes');
                    end
                    
                    % Guarda los puntos
                    lim1 = ceil(obj.fs*lim1);
                    lim2 = ceil(obj.fs*lim2);
                    obj.window_fix{i} = [lim1, lim2];
                    
                end
            else
                % Usa todo el registro
                obj.window_fix = {[1, length(obj.cable)]};
                obj.close_prev_plot();
            end
            
        end % request_windows function
        
        function tabularAnalisisIdentificacionNL(obj, tabla) %#ok<INUSL>
            % tabularAnalisisIdentificacionNL: Imprime en consola la
            % tabla de resultados de la identificacion no lineal
            
            % Obtiene numero de modos
            [n, ~] = size(tabla);
            fprintf('\t\t|\tN\t|\tFreq\t|\tBeta\t|\tTheta\t|\tRho\t\t|\n');
            fprintf('\t\t---------------------------------------------------------\n');
            for i = 1:n
                fprintf('\t\t|\t%d\t|\t%.4f\t|\t%.4f\t|\t%.4f\t|\t%.4f\t|\n', ...
                    i, tabla(i, 1)/(2 * pi), tabla(i, 2), tabla(i, 3), tabla(i, 4));
            end % for i
            fprintf('\t\t---------------------------------------------------------\n');
            
        end % tabularAnalisisIdentificacionNL function
        
        function tabularNLFreq(obj, tabla) %#ok<INUSL>
            % tabularNLFreq: Imprime la tabla de frecuencias
            
            % Obtiene numero de modos
            [n, ~] = size(tabla);
            fprintf('\t\t|\tN\t|\tPSD\t\t|\t\tPromedio\t\t|\tMediana\t|\n');
            fprintf('\t\t---------------------------------------------------------\n');
            for i = 1:n
                fprintf('\t\t|\t%d\t|\t%.4f\t|\t%.4f\t+-\t%.4f\t|\t%.4f\t|\n', ...
                    i, tabla(i, 1), tabla(i, 2), tabla(i, 3), tabla(i, 4));
            end % for i
            fprintf('\t\t---------------------------------------------------------\n');
            
        end % tabularNLFreqAmort function
        
        function tabularNLAmort(obj, tabla) %#ok<INUSL>
            % tabularNLAmort: Imprime la tabla de amortiguamientos
            
            % Obtiene numero de modos
            [n, ~] = size(tabla);
            fprintf('\t\t|\tN\t|\t\tPromedio\t\t|\tMediana\t|\n');
            fprintf('\t\t---------------------------------------------\n');
            for i = 1:n
                fprintf('\t\t|\t%d\t|\t%.4f\t+-\t%.4f\t|\t%.4f\t|\n', ...
                    i, tabla(i, 1), tabla(i, 2), tabla(i, 3));
            end % for i
            fprintf('\t\t---------------------------------------------\n');
            
        end % tabularNLAmort function
        
    end
end