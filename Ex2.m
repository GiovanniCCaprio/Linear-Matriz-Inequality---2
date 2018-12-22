%% Exercicio 2
clc, clear all, close all

% Entradas

n = 3;          % numero de estados
m = 2;          % numero de entradas
l = 1;          % numero de saidas
N = 4;          % numero de vértices


% Cria os A estáveis (Autovalores negativos com d < 1, diagonal princiapal)
for i = 1:N              % A = n x n
    A{i} = zeros(n);
    for x = 1:n
        A{i}(x,x) = -rand;
    end
end

% Cria B
for i = 1:N              % B = n x m
    B{i}(1,1) = rand;
    for x = 1:n
        for t = 1:m
            B{i}(x,t) = rand;
        end
    end
end

% Cria C
for i = 1:N              % B = l x n
    C{i}(1,1) = rand;
    for x = 1:l
        for t = 1:n
            C{i}(x,t) = rand;
        end
    end
end

% Cria D
for i = 1:N              % D = l x m
    D{i}(1,1) = rand;
    for x = 1:l
        for t = 1:m
            D{i}(x,t) = rand;
        end
    end
end

q = 0; %verificação final 
h2 = 0;
hinf = 0;
maxEig = -10^10;

%%   GRID


    % 1000 pontos uniformemente distribuidos

q = 0; %verificação final 
maxEig = -10^10;
for i = 1:1000
    A_ = zeros(n);
    B_ = zeros(n,m);
    C_ = zeros(l,n);
    D_ = zeros(l,m);
    %criação dos pontos uniformes
    x(1)=1-rand^(1/(N-1));
        for k=2:N-1
            x(k)=(1-sum(x(1:k-1)))*(1-rand^(1/(N-k)));
        end
    x(N) = 1-sum(x(1:N-1));
        for t = 1:N
            A_ = A_ + (x(t)*A{t});
            B_ = B_ + (x(t)*B{t});
            C_ = C_ + (x(t)*C{t});
            D_ = D_ + (x(t)*D{t});
        end
        sys = ss(A_,B_,C_,D_,-1);
        h2_ = norm(sys,2);
        hinf_ = norm(sys,inf);
        autoValores = eig(A_);
        plot(real(autoValores),imag(autoValores),'*')
        xlabel('Valor Real')
        ylabel('Imaginário')
        hold on
        maximo = max(abs(autoValores));
        if maximo > maxEig
            maxEig = maximo;
            q = 0;
        end
        if h2_ > h2
            h2 = h2_;
        end
        if hinf_ > hinf
            hinf = hinf_;
        end
end


    %% 100 pts igualmente espaçados em 3 vertices

autoValores = 0;
for t = 1:N %varre os A
    for k = (t+1):N
        if t~=k %verifica se não é o mesmo vértice
            for i = 0:0.01:1
                A_=i*A{t}+(1-i)*A{k};
                B_=i*B{t}+(1-i)*B{k};
                C_=i*C{t}+(1-i)*C{k};
                D_=i*D{t}+(1-i)*D{k};
                autoValores = eig(A_);
                plot(real(autoValores),imag(autoValores),'*')
                hold on
                sys = ss(A_,B_,C_,D_,-1);
                h2_ = norm(sys,2);
                hinf_ = norm(sys,inf);
                %maior valor real
                maximo = max(abs(autoValores));
                if maximo > maxEig
                    maxEig = maximo;
                    q = 0;
                end
                if h2_ > h2
                    h2 = h2_;
                end
                if hinf_ > hinf
                    hinf = hinf_;
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
                                A_ = z(1)*A{t}+z(2)*A{k}+z(3)*A{l};
                                B_ = z(1)*B{t}+z(2)*B{k}+z(3)*B{l};
                                C_ = z(1)*C{t}+z(2)*C{k}+z(3)*C{l};
                                D_ = z(1)*D{t}+z(2)*D{k}+z(3)*D{l};
                                autoValores = eig(A_);
                                %maior valor real
                                maximo = max(abs(autoValores));
                                if maximo > maxEig
                                    maxEig = maximo;
                                    q = 0;
                                end
                                if h2_ > h2
                                    h2 = h2_;
                                end
                                if hinf_ > hinf
                                    hinf = hinf_;
                                end
                                plot(real(autoValores),imag(autoValores),'*')
                                hold on
                           end
                        end
                    end
            end
        end
else
    disp('Número de vértices menor que 4')
end

%% Saidas

maxEig
h2
hinf
for i = 1:N
   A{i}
   B{1}
   C{1}
   D{1}
end    
