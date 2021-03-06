% Load data forzará una eliminación total
addpath('lib');
clear all; %#ok<*CLALL>
close all;

archivos = {; ...
    'Cable_1_impacto.mat'; ...
    'Cable_2_ambiental.mat'; ...
    'Cable_3_impacto.mat'; ...
    'Cable_4_ambiental.mat'; ...
    'Cable_5_impacto.mat'; ...
    };
% archivos = {; ...
%     'Cable_1_ambiental.mat'; ...
%     'Cable_2_impacto.mat'; ...
%     'Cable_3_ambiental.mat'; ...
%     'Cable_4_impacto.mat'; ...
%     'Cable_5_ambiental.mat'; ...
%     };
tot = length(archivos);
data = cell(tot, 1);
test_names = cell(tot, 1);

%% Carga los datos de los registros
fprintf('Cargando archivos ... ');
for i = 1:tot
    fieldname = lower(archivos{i});
    fieldname = replace(fieldname, '.mat', '');
    fieldname = replace(fieldname, 'cable_', 'cable');
    fieldname = strcat('data_', fieldname);
    data{i} = load(sprintf('data/%s', archivos{i}));
    data{i} = data{i}.(fieldname);
    test_names{i} = replace(replace(archivos{i}, '.mat', ''), '_', ' ');
end
fprintf('OK\n');

%% Informacion de los sensores
sensores_loc = {; ...
    [4, 4]; ...
    [105; 105]; ...
    [138; 243]; ...
    [119; 362]; ...
    [106; 468]; ...
    [134; 602]; ...
    [124; 726]; ...
    [130; 856]; ...
    };

% Factor en A*x+B [A, B]
sensor_factor = {; ...
    [980.665/1028, 0]; ... # mV/g -> cm/s^2
    [2.7473; 2.9293]; ... # cm
    [2.7423; 3.0162]; ...
    [2.7539; 2.8831]; ...
    [2.7663; 2.5769]; ...
    [2.7502; 2.6691]; ...
    [2.7474; 2.6111]; ...
    };

%% Datos del cable
cable_largo = 8.56;
cable_densidad = 0.05; % kg/m
cable_diametro = 8.00/1000; % mm

%% Aplica los factores de conversion
for i = 1:tot
    registro = data{i};
    for j = 1:7
        registro(:, j) = registro(:, j) .* sensor_factor{j}(1) + sensor_factor{j}(2);
    end
    data{i} = registro;
end

%% Limpia las variables
clearvars fieldname i j registro;