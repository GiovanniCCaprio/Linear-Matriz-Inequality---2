    %% Exercicio 1 
   
    clc 
    clear all
    close all

    %% Entrada

n = 3;       % numero de estados do sistema
N = 4;       % numero de vertices
caso =  1;   % 0 = continuo ; 1 = discreto;
for i = 1:N                          % criação dos A(a) para teste
    A{i} = randn(n);
end

    %% 1000 pontos uniformemente distribuidos

q = 0; %verificação final 
maxEig = -10^10;
for i = 1:1000
    Aa = zeros(n);
    %criação dos pontos uniformes
    x(1)=1-rand^(1/(N-1));
        for k=2:N-1
            x(k)=(1-sum(x(1:k-1)))*(1-rand^(1/(N-k)));
        end
    x(N) = 1-sum(x(1:N-1));
        for t = 1:N
            Aa = Aa + (x(t)*A{t});
        end
    autoValores = eig(Aa);
        %maior valor real
        if caso == 0
            maximo = max(real(autoValores));
        else
            maximo = max(abs(autoValores));
        end
        if maximo > maxEig
            maxEig = maximo;
            q = 0;
            for s = 1:N
                g(s) = x(s);
            end
        end
    plot(real(autoValores),imag(autoValores),'*')
    xlabel('Valor Real')
    ylabel('Imaginário')
    hold on
end


    %% 100 pts igualmente espaçados em 3 vertices


for t = 1:N %varre os A
    for k = (t+1):N
        if t~=k %verifica se não é o mesmo vértice
            for i = 0:0.01:1
                Aa=i*A{t}+(1-i)*A{k};
                autoValores = eig(Aa);
                plot(real(autoValores),imag(autoValores),'*')
                hold on
                %maior valor real
                if caso == 0
                    maximo = max(real(autoValores));
                else
                    maximo = max(abs(autoValores));
                end
                if maximo > maxEig
                    maxEig = maximo; 
                    q = 1;
                    g(1) = t;
                    g(2) = k;
                    g(3) = i;
                    g(4) = 1-i;
                end
            end
        end
    end
end    

    %% 150 pts uniformemente distribuidos em 3 vertices

if N>3  %verifica se os vértices são maiores que 3
        for t = 1:N     %Varre A
            for k = (1+t):N
                    for l = (1+k):N
                        if (l~=k||l~=t||t~=k)     %só verifica se forem 3 vértices diferentes
                           for p = 1:150 
                                %criação dos pontos uniformes
                                test = 2;
                                while test > 1
                                    z(1) = (rand);
                                    z(2) = rand;
                                    test = sum(z(1:2));
                                end
                                z(3) = 1-sum(z(1:2));
                                Aa= z(1)*A{t}+z(2)*A{k}+z(3)*A{l};
                                autoValores = eig(Aa);
                                %maior valor real
                                if caso == 0
                                    maximo = max(real(autoValores));
                                else
                                    maximo = max(abs(autoValores));
                                end
                                if maximo > maxEig
                                    maxEig = maximo; 
                                    q = 2;
                                    g(1) = t;
                                    g(2) = k;
                                    g(3) = l;
                                    g(4) = z(1);
                                    g(5) = z(2);
                                    g(6) = z(3);
                                end
                                plot(real(autoValores),imag(autoValores),'*')
                                hold on
                           end
                        end
                    end
            end
        end
else
    disp('Número de vértices menor que 3')
end

%% Lógica para display dos resultados

disp('Maior valor real:' )
disp(maxEig)

disp('Este valor foi encontrado nos alphas:' )

if q == 0
   for s = 1:N
       disp(g(s));
   end
end
if q == 1
    disp(g(3))
    disp(g(4))
    disp('Entre os vértices(A) de número:' )
    disp(g(1))
    disp(g(2))
end
if q == 2
    disp(g(4))
    disp(g(5))
    disp(g(6))
    disp('Entre os vértices(A) de número:' )
    disp(g(1))
    disp(g(2))
    disp(g(3))
end
if caso == 1
    tetha=-pi:.055:pi;
    xr=cos(tetha);
    yr=sin(tetha);
    plot(xr,yr,'black');    
end









