%% Exercicio 5
clc,clear all,close all

% Dados 

A = [0 1; -16 -0.8];
B = [1;1];
C = [1 2];
D = 0;

sys = ss(A,B,C,D);
%% A)

% Cálculo da energia de Entrada

w = 3;
T = 0.01;
t = 0:T:10;
sinal_entrada = sin(w*t).*exp(-0.1*t);
energia_entrada = sqrt(trapz(sinal_entrada.^2)*T);  

% Cálculo da energia de Saída

sinal_saida = lsim(sys,sinal_entrada,t);  % sinal de saida dependente de u,t,sys
energia_saida = sqrt(trapz(sinal_saida.^2)*T);  

% Cáculo da relação entre Saída/Entrada
relacao =  energia_saida/energia_entrada


%% B)

infinito(A,B,C,D) % Roda a função que cálcula a LMI e informa a norma Hinf

% Verificação da simulação com a nomr Hinf
for k = 0:0.05:0.4
    vl = 20*k + 1;
    w = 3.8+k;
    T = 0.01;
    t = 0:T:10;
    sinal_entrada = sin(w*t).*exp(-0.1*t);
    energia_entrada = sqrt(trapz(sinal_entrada.^2)*T);  

    % Cálculo da energia de Saída
    sinal_saida = lsim(sys,sinal_entrada,t);  % sinal de saida dependente de u,t,sys
    energia_saida = sqrt(trapz(sinal_saida.^2)*T);  

    % Cáculo da relação entre Saída/Entrada
    disp('A relação entre a energia de saída e entrada é:')
    vl
    relacao_total =  energia_saida/energia_entrada
    w
    
end

%% C)

% Cálculo da energia de Saída

sinal_saida = impulse(sys);  % sinal de saida dependente de u,t,sys
energia_H2 = sqrt(trapz(sinal_saida).^2*T) 



%% D)
    
H2(A,B,C,D) % Roda a função que cálcula a LMI e informa a norma H2
