% dTheta = 1e-5:1e-1:1; Gaussian Model mit fitrgp(X,y); X ist Matrix mit
% input variablen, y output variable (Viskosität)
% Variablen noch so weit wie möglich reduzieren -> nur für eine
% Extrudergeometrie, d.h. R_e, R_i vorgeben
% t = 0.05;
% dTheta = 1;
% dP = 2e5;
% N = 5;
% phi = 0.5;
% T = 150+273.15;
% w_p = 1;
% M_w = 75000;
% 
% R_e = 5.5e-3;
% 
% eta0 = calcViscosity0(M_w,w_p,T);
% [eta, gamma] = calcViscosityPLA(dTheta, dP, t, N, phi, eta0);

%% Input Matrix erstellen
% sinnvolle Grenzen festlegen
R_e = 5.5e-3;
t = [R_e; 2*R_e; 3*R_e]';
dTheta = linspace(0.1,1,10);
dP = logspace(4,7,10);
N = linspace(1/60,500/60,10);
phi = linspace(0.1,1,10);
eta0 = logspace(-2,5,10);

X = combinations(dTheta,dP,t,N,phi,eta0);
X_mat = table2array(X);
Y = zeros(length(X_mat(:,1)),1);
for i = 1 : length(X_mat(:,1))
    Y(i) = calcViscosityPLA(X_mat(i,1),X_mat(i,2),X_mat(i,3),X_mat(i,4),X_mat(i,5),X_mat(i,6));
end

gprMdl = fitrgp(X_mat,Y,'KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','BasisFunction','pureQuadratic');

ypred = resubPredict(gprMdl);
x_lin = [0 1e5];
err = 0.1;
x_lin_ub = x_lin.*(1+err);
x_lin_lb = x_lin.*(1-err);

figure
plot(Y,ypred,'x',x_lin,x_lin_lb,x_lin,x_lin_ub,x_lin,x_lin)
xlabel('Real model')
ylabel('Gaussian model')

figure
plot(X_mat(1:5,6),Y(1:5),X_mat(1:5,6),ypred(1:5))
legend('physikalisches Modell','Gauss Modell')
xlabel('eta_0')
ylabel('eta')

function eta_0 = calcViscosity0(M_w,w_p,T)
    E_a = 104000;
    R = 8.314;
    T_ref = 473.15;
    a = 2.24;
    eta_0_ref = 97;    
    eta_0 = eta_0_ref * power(w_p,a) * power(M_w/100000, a) * exp(E_a/R*(1/T-1/T_ref));
end

function [eta, gamma] = calcViscosityPLA(dTheta, dP, t, N, phi, eta_0)
    B = 1.67e-5;
    R_e = 5.5e-3;
    R_i = 3.5e-3;
        
    %Declaration of necessary variables
    dz = t/(2*pi)*dTheta;
    v_e_ax = N*2*pi*R_e*cos(phi);
    v_e_u = N*2*pi*R_e*sin(phi);
    
    %eta_0 = 0.1;
    
    %Initialization for iterative procedure
    eta = eta_0;
    dE = 100;
    %Iterative Loop with 10^-4 tolerance
    while dE > 0.001
    
        %Berechung der Schergeschwindigkeit nach Vergnes 1998
        fun = @(r) r*sqrt((1/(2*eta)*(dP)/dTheta*(1-2/r^2*R_e^2* ...
            R_i^2/(R_e^2-R_i^2)*log(R_e/R_i))+2*v_e_ax*R_e/(R_e^2-R_i^2)* ...
            (1-(r^2-R_i^2)/(r^2)))^2+(1/(4*eta)*(dP)/dz*(2*r-1/r* ...
            (R_e^2-R_i^2)/(log(R_e/R_i)))+v_e_u/(r*log(R_e/R_i)))^2);
    
        gamma = abs(2/(R_e^2-R_i^2)*integral(fun, R_e,R_i, "ArrayValued",1));
    
        %Viskositätsberechnung nach Witzke 1999
        eta1 = asinh(eta_0*gamma*B)/(gamma*B);
    
    
        dE = abs(eta1-eta);
        eta = eta1;
    end
end