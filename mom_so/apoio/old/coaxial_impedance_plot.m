% Definir os parâmetros
clc, clear;
frequencies = logspace(0, 6, 100); % Frequências de 10 Hz a 100 MHz
a = 22e-3; % Raio do condutor interno em metros
b = 39.5e-3; % Raio interno da bainha em metros
c = 44e-3; % Raio externo da bainha em metros
sigma = 5.8e7; % Condutividade do cobre em S/m

% Calcular a impedância
[R, L] = coaxial_impedance(frequencies, a, b, c, sigma);

% Plotar a indutância
figure;
subplot(1,2,1);
loglog(frequencies, L);
xlabel('Frequency (Hz)');
ylabel('Inductance p.u.l. [H/m]');
xlim([1e0, 1e6]);
ylim([1.1e-7, 1.9e-7]);
title('Inductance of Coaxial Cable');
grid on;

% Plotar a resistência
subplot(1,2,2);
loglog(frequencies, R);
xlabel('Frequency (Hz)');
ylabel('Resistance p.u.l. [\Omega/m]');
xlim([1e0, 1e6]);
% ylim([8e-5, 1e-2]);
title('Resistance of Coaxial Cable');
grid on;