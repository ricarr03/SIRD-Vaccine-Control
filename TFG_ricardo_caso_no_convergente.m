% Este código ha sido escrito por Ricardo Rodríguez León
% estudiante de 4o de Ingeniería Matemática para el Trabajo de 
% Fin de Grado de Modelización Matemática y Control de Epidemias.

clear all
clc
clf
% constantes asociadas a nuestro problema de control 
B = 0.2491;
mu = 0.0121;
gamma = 0.1210;
omega = 0.1177;
Nu = 50000;
N = 1000000;
sigma = 0.05;
c = 10;
k = 500;

% variable a evaluar para salir o no del bucle de resolución
test = -1;
% valor mínimo para la condición de convergencia
num = 10^-6;
% definición del paso temporal h 
q = 0; d = 200;
M = 10^5; h = d/M;

% funcion PHI
Phi = @(x)[(x<=0).*(0) + (0<x & x<1/(2*N)).*(exp(-((1-2*N*x).^2)./(1-(1-2*N*x).^2))) + (x>=1/(2*N)).*(1)];
% derivada de PHI
dPhi = @(x)[(x<0).*(0) + (0<x & x<1/(2*N)).*((exp(-((1-2*N*x).^2)./(1-(1-2*N*x).^2))*(1-2*N*x))./(4*N*x^2*(N*x-1).^2))+(x>1/(2*N)).*(0)];

% funcion F
F = @(t,X,u)[mu*(X(2,1)+X(3,1))+omega*gamma*X(2,1)+sigma*X(3,1)-B*X(1,1)*X(2,1)-Nu/N*u*Phi(X(1,1));...
    B*X(1,1)*X(2,1)-mu*X(2,1)-gamma*X(2,1);...
    (1-omega)*gamma*X(2,1)-sigma*X(3,1)-mu*X(3,1)+Nu/N*u*Phi(X(1,1));...
    omega*gamma*X(2,1)];

% funcion G
G = @(t,X,P,u)[[dPhi(X(1,1))*(Nu/N*u*(P(3,1)-P(1,1))-c*(u).^2*Phi(X(1,1)))+B*X(2,1)*(P(2,1)-P(1,1));...
    omega*gamma*(P(1,1)+P(4,1)-P(3,1))+(P(1,1)-P(2,1))*mu+B*X(1,1)*(P(2,1)-P(1,1))+gamma*(P(3,1)-P(2,1));...
    (mu+sigma)*(P(1,1)-P(3,1)); -k]];

% datos iniciales
t(1) = q;
X(1,1) = 0.99; X(2,1) = 0.01;
X(3,1) = 0; X(4,1) = 0;
% datos iniciales para la Q
Q(1,1) = 0; Q(2,1) = 0;
Q(3,1) = 0; Q(4,1) = -1;
% datos finales 
t(M) = d;

% valor del control para resolver el sistema sin control
u(1:M) = 0;

% resolvemos el sistema sin control 

for r=1:(M-1)
    t(r+1) = t(r) + h;
    F1 = F(t(r),X(:,r),u(r));
    F2 = F(t(r)+h/2,X(:,r)+h*F1/2,1/2*(u(r)+u(r+1)));
    F3 = F(t(r)+h/2,X(:,r)+h*F2/2,1/2*(u(r)+u(r+1)));
    F4 = F(t(r)+h,X(:,r)+h*F3,u(r+1));
    X(:,r+1) = X(:,r) + h/6 * (F1 + 2*F2 + 2*F3 + F4);
end

% cálculo del coste sin control
coste_sin = X(4,M);
for r=1:M
    coste_sin = coste_sin + (h*k*X(4,r));
end
% representación de la evolución sin control
plot(t,X(1,:),'r',t,X(2,:),'b',t,X(3,:),'g',t,X(4,:),'k')
legend('S','I','R','D','FontSize',14)
title('Evolución del sistema sin control','FontSize',14)
xlabel('Tiempo (días)','FontSize',14)
ylabel('Porcentaje de individuos','FontSize',14)

% resolvemos el sistema con el criterio de convergencia
% valor inicial de u
u(1:M) = 0;
% contador del número de iteraciones
iter = 0;
while test < 0
  
    iter = iter + 1;
    % definimos las variables old (para comparar convergencia)
    oldu(1:M) = u(1:M);
    % definición de las variables de coste con control
    coste_control(iter) = X(4,M);
    for r=1:M
        % cálculo del coste con control para cada iteracion
        coste_control(iter) = coste_control(iter) + (h* (k*X(4,r)+c/2*(Phi(X(1,r))*u(r))^2));
    end
    % forward
    for r=1:(M-1)
        t(r+1) = t(r) + h;
        F1 = F(t(r),X(:,r),u(r));
        F2 = F(t(r)+h/2,X(:,r)+h*F1/2,1/2*(u(r)+u(r+1)));
        F3 = F(t(r)+h/2,X(:,r)+h*F2/2,1/2*(u(r)+u(r+1)));
        F4 = F(t(r)+h,X(:,r)+h*F3,u(r+1));
        X(:,r+1) = X(:,r) + h/6 * (F1 + 2*F2 + 2*F3 + F4);
    end
    % backwards

    for r=1:(M-1)
        t(r+1) = t(r) + h;
        Q1 = G(t(r),X(:,M-r+1),Q(:,r),u(M-r+1));
        Q2 = G(t(r)+h/2,1/2*(X(:,M-r+1)+X(:,M-r)),Q(:,r)+h*Q1/2,1/2*(u(M-r+1)+u(M-r)));
        Q3 = G(t(r)+h/2,1/2*(X(:,M-r+1)+X(:,M-r)),Q(:,r)+h*Q2/2,1/2*(u(M-r+1)+u(M-r)));
        Q4 = G(t(r)+h,X(:,M-r),Q(:,r)+h*Q3,u(M-r));
        Q(:,r+1) = Q(:,r) + h/6 * (Q1 + 2*Q2 + 2*Q3 + Q4);
    end
    % cálculo de P por medio de la fórmula
    for r=1:M
        P(1,r) = Q(1,M-r+1); P(2,r) = Q(2,M-r+1);
        P(3,r) = Q(3,M-r+1); P(4,r) = Q(4,M-r+1);
    end
    % cálculo del control
    for r=1:M
       u(r) = max(min(((Nu/N*(P(3,r)-P(1,r)))./(c*Phi(X(1,r)))),1),0);
    end
    
    % definición de las variables de coste con control y de la diferencia
    % de u
    aprox(iter) = 0;
    for r=1:M
        % cálculo de la diferencia de los controles para ver la
        % convergencia para cada iteracion
        aprox(iter) = aprox(iter) + (h* (u(r)-oldu(r))^2);
    end

    % criterio de convergencia
    if aprox(iter) <= num 
        test = 1;
        disp("converge");
        disp(iter);
    % para obtener la gráfica de iteraciones impares, basta cortar en una
    % iteración impar
    elseif iter==100
        test = 1;
        disp("iteracion");
        disp(iter);
    end
end

% representación de la evolución de infectados tras resolver el sistema
plot(t,X(1,:),'r',t,X(2,:),'b',t,X(3,:),'g',t,X(4,:),'k')
legend('S','I','R','D','FontSize',14)
title('Evolución del sistema con control (PAR)','FontSize',14)
title('Evolución del sistema con control (IMPAR)','FontSize',14)
xlabel('Tiempo (días)','FontSize',14)
ylabel('Porcentaje de individuos','FontSize',14)
xlim([0 200])



