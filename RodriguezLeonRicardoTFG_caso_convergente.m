% Este código ha sido escrito por Ricardo Rodríguez León
% estudiante de 4o de Ingeniería Matemática para el Trabajo de 
% Fin de Grado de Modelización Matemática y Control de Epidemias.

clear all
clc
clf
% constantes asociadas a nuestro problema de control
B = 0.2691;
mu = 0.0121;
gamma = 0.1001;
omega = 0.1252;
Nu = 50000;
N = 1000000;
sigma = 0.0798;
c = 10;
k = 500;

% variable a evaluar para salir o no del bucle de resolución
test = -1;
% valor mínimo para la condición de convergencia
num = 10^-6;
% definición del paso temporal h 
q = 0; d = 200;
M = 10^4; h = d/M;

% funcion PHI
Phi = @(x)[(x<=0).*(0) + (0<x & x<1/(2*N)).*(exp(-((1-2*N*x).^2)./(1-(1-2*N*x).^2))) + (x>=1/(2*N)).*(1)];
% Phi = @(x)[1];
% derivada de PHI
dPhi = @(x)[(x<0).*(0) + (0<x & x<1/(2*N)).*((exp(-((1-2*N*x).^2)./(1-(1-2*N*x).^2))*(1-2*N*x))./(4*N*x^2*(N*x-1).^2))+(x>1/(2*N)).*(0)];
% dPhi = @(x)[0];

% funcion F
F = @(t,X,u)[mu*(X(2,1)+X(3,1))+omega*gamma*X(2,1)+sigma*X(3,1)-B*X(1,1)*X(2,1)-Nu/N*u*Phi(X(1,1));...
    B*X(1,1)*X(2,1)-mu*X(2,1)-gamma*X(2,1);...
    (1-omega)*gamma*X(2,1)-sigma*X(3,1)-mu*X(3,1)+Nu/N*u*Phi(X(1,1));...
    omega*gamma*X(2,1)];

% funcion G
G = @(t,X,P,u)[[dPhi(X(1,1))*(Nu/N*u*(P(3,1)-P(1,1))-c*(u).^2*Phi(X(1,1)))+B*X(2,1)*(P(2,1)-P(1,1));...
    omega*gamma*(P(1,1)+P(4,1)-P(3,1))+(P(1,1)-P(2,1))*mu+B*X(1,1)*(P(2,1)-P(1,1))+gamma*(P(3,1)-P(2,1));...
    (mu+sigma)*(P(1,1)-P(3,1)); -k]];

% datos iniciales (por defecto el caso inicial)
% si se desea otras condiciones iniciales cambiarlas por esas mismas
t(1) = q;
X(1,1) = 0.99; X(2,1) = 0.01;
X(3,1) = 0; X(4,1) = 0;
% datos iniciales para la Q
Q(1,1) = 0; Q(2,1) = 0;
Q(3,1) = 0; Q(4,1) = -1;
% datos finales 
t(M) = d;

% resolvemos el sistema con el criterio de convergencia
% valor inicial de u
% u(1:M) = 1;
% u(1:M) = sin(pi*(1:M)/M);
u(1:M) = 0;
% contador del número de iteraciones, se inicia en 0
iter = 0;
% valor de coste_min muy alto para poder extraer el coste mínimo.
coste_min = 10^7;
while test < 0
    % La iteración 1, corresponde con i=0, siguiendo la notación de la
    % memoria
    iter = iter + 1;
    
    % definimos las variables old (para comparar convergencia)
    oldu(1:M) = u(1:M);
    
    % forward
    for r=1:(M-1)
        t(r+1) = t(r) + h;
        F1 = F(t(r),X(:,r),u(r));
        F2 = F(t(r)+h/2,X(:,r)+h*F1/2,1/2*(u(r)+u(r+1)));
        F3 = F(t(r)+h/2,X(:,r)+h*F2/2,1/2*(u(r)+u(r+1)));
        F4 = F(t(r)+h,X(:,r)+h*F3,u(r+1));
        X(:,r+1) = X(:,r) + h/6 * (F1 + 2*F2 + 2*F3 + F4);
    end
    
    % definición de las variables de coste con control
    J(iter) = X(4,M);
    for r=1:M
        % cálculo del coste con control para cada iteracion
        J(iter) = J(iter) + (h* (k*X(4,r)+c/2*(Phi(X(1,r))*u(r))^2));
    end
    
    if iter == 1
        % comentar y descomentar para obtener las gráficas deseadas
        % representación de la evolución de fallecidos sin control
        plot(t,X(1,:),'r')
        % representación de la evolución sin control
        % plot(t,X(1,:),'r',t,X(2,:),'b',t,X(3,:),'g',t,X(4,:),'k')
        % title('Evolución sistema sin control','FontSize',14)
        % legend('S','I','R','D','FontSize',14)
        % xlabel('Tiempo (días)','FontSize',14)
        % ylabel('Porcentaje de individuos','FontSize',14)
        hold on
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
    
    aprox(iter) = 0;
    for r=1:M
        % cálculo de la diferencia de los controles para ver la
        % convergencia para cada iteracion
        aprox(iter) = aprox(iter) + (h* (u(r)-oldu(r))^2);
    end
    
    % guardamos el valor del sistema y del control para el valor del mínimo
    % coste
    if iter>1
        if J(iter)<coste_min
            xmin((1:4),(1:M)) = X((1:4),(1:M));
            umin(1:M) = oldu(1:M);
            coste_min = J(iter);
            iter_min = iter;
        end
    end    
    % criterio de convergencia
    if aprox(iter) <= num 
        test = 1;
        disp("converge");
        disp(iter);

    elseif iter==100
        test = 1;
        disp("iteracion");
        disp(iter);
    end

end


% comentar y descomentar para obtener las gráficas deseadas
% representación de la evolución de fallecidos tras resolver el sistema
plot(t,X(1,:),'b',t,xmin(1,:),'g')
% plot(t,X(1,:),'r',t,X(2,:),'b',t,X(3,:),'g',t,X(4,:),'k')
% title('Evolución del sistema con control','FontSize',14)
% legend('S','I','R','D','FontSize',14)
% xlabel('Tiempo (días)','FontSize',14)
% ylabel('Porcentaje de individuos','FontSize',14)
title("Evolución de (S)",'FontSize',14)
xlabel("Tiempo (días)",'FontSize',14)
ylabel("Susceptibles",'FontSize',14)
ultimovalor = strcat('con último control (u_{',num2str(iter-1),'})');
minimo = strcat('mínimo coste (u_{',num2str(iter_min-1),'})');
legend("sin control (u=0)", ultimovalor, minimo,'FontSize',14)
xlim([0 200])
