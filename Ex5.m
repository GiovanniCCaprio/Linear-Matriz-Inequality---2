%% Exercicio 5
clc,clear all,close all

% Dados 

A = [0 1; -16 -0.8];
B = [1;1];
C = [1 2];
D = 0;

sys = ss(A,B,C,D);
%% A)

% C�lculo da energia de Entrada

w = 3;
T = 0.01;
t = 0:T:10;
sinal_entrada = sin(w*t).*exp(-0.1*t);
energia_entrada = sqrt(trapz(sinal_entrada.^2)*T);  

% C�lculo da energia de Sa�da

sinal_saida = lsim(sys,sinal_entrada,t);  % sinal de saida dependente de u,t,sys
energia_saida = sqrt(trapz(sinal_saida.^2)*T);  

% C�culo da rela��o entre Sa�da/Entrada
relacao =  energia_saida/energia_entrada


%% B)

infinito(A,B,C,D) % Roda a fun��o que c�lcula a LMI e informa a norma Hinf

% Verifica��o da simula��o com a nomr Hinf
for k = 0:0.05:0.4
    vl = 20*k + 1;
    w = 3.8+k;
    T = 0.01;
    t = 0:T:10;
    sinal_entrada = sin(w*t).*exp(-0.1*t);
    energia_entrada = sqrt(trapz(sinal_entrada.^2)*T);  

    % C�lculo da energia de Sa�da
    sinal_saida = lsim(sys,sinal_entrada,t);  % sinal de saida dependente de u,t,sys
    energia_saida = sqrt(trapz(sinal_saida.^2)*T);  

    % C�culo da rela��o entre Sa�da/Entrada
    disp('A rela��o entre a energia de sa�da e entrada �:')
    vl
    relacao_total =  energia_saida/energia_entrada
    w
    
end

%% C)

% C�lculo da energia de Sa�da

sinal_saida = impulse(sys);  % sinal de saida dependente de u,t,sys
energia_H2 = sqrt(trapz(sinal_saida).^2*T) 



%% D)
    
H2(A,B,C,D) % Roda a fun��o que c�lcula a LMI e informa a norma H2
