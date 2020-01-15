% Realiza ajuste no lineal
test = 1; %#ok<*IJCL>

% Carga el ejemplo
ex = psd{test}; %#ok<*SUSENS>

nmodos = {[1, 2, 3, 4, 5], [], [1, 2, 3, 4, 5], [], [1, 2, 3]};
ex.calc_nl('modos', nmodos{test}, 'betalim', [0, 5e-3], 'load', true); % La primera ventana