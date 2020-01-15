% Con el método PSD determine las frecuencias de vibrar de cada cable.
% Para revisar su procedimiento adjunte gráficas de PSD vs frecuencias o
% log (PSD) vs frecuencias para cada set de datos entregados.

% Filter data
ensayo_filter = [650, 0, 100, 100, 200];

% En 2: Borrar registros 3 y 7
% En 4: Borrar 6

% Creamos los datos
if exist('psd', 'var')
    psd = cell(tot, 1);
end

for i=1:tot
    psd{i} = PSDAnalsys(data{i}, test_names{i}, ensayo_filter(i)); %#ok<*UNRCH,*SUSENS>
end
% i=1; psd{i} = PSDAnalsys(data{i}, test_names{i}, ensayo_filter(i));
clearvars i f;