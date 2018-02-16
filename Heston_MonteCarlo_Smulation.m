clear;
clc;
close all;

S = 100; 
K = 100; 
T = 5;  
r = 0.05;
v = 0.09;  % initial var
theta = 0.09;  
kappa = 2;  
omega = 0.3;  
rho = -0.3; 
x = 34.9998;  
 
path = 10000;  %sample paths
iterations = 50;   
N = 20;   
 
dt = 1/N;  

% Absorption
S_A = zeros(path,T*N+1);
V_At = zeros(path,T*N+1);
S_A(:,1) = S;
V_At(:,1) = v;
V_A(:,1) = v;
V_A = zeros(path,T*N+1);
C_A = zeros(1,iterations);
LS_A=zeros(path,T*N+1);


tic;
for i = 1:iterations
    dW1 = sqrt(dt)*randn(path,T*N);
    Z = randn(path,T*N);
    dW2 = rho*dW1+sqrt(1-rho^2)*Z*sqrt(dt);
    for j = 1:T*N
        V_At(:,j+1) = max(V_At(:,j),0) - kappa*(max(V_At(:,j),0)-theta)*dt + omega*sqrt(max(V_At(:,j),0)).*dW1(:,j);
        V_A(:,j+1) = max(V_At(:,j+1),0);
        LS_A(:,j+1) = log(S_A(:,j))+(r-(1/2).*V_A(:,j))*dt+sqrt(V_A(:,j)).*dW2(:,j);
        S_A(:,j+1)=exp(LS_A(:,j+1));
    end
    
    C_A(i) = exp(-r*T)*(mean(max(S_A(:,end)-K,0)));
    
end
Price_Absorption = mean(C_A)
SE_Absorption = sqrt(mean((C_A-mean(C_A)).^2))
Bias_Absorption = abs(Price_Absorption-x)
RMSE_Absorption = sqrt(Bias_Absorption^2+SE_Absorption^2)
toc;
 
% Reflection 
S_R = zeros(path,T*N+1);
V_Rt = zeros(path,T*N+1);
S_R(:,1) = S;
V_Rt(:,1) = v;
V_R(:,1) = v;
V_R = zeros(path,T*N+1);
C_R = zeros(1,iterations);
LS_R=zeros(path,T*N+1);
tic;
for k = 1:iterations
    dW1 = sqrt(dt)*randn(path,T*N);
    Z = randn(path,T*N);
    dW2 = rho*dW1+sqrt(1-rho^2)*Z*sqrt(dt);
    for j = 1:T*N
        V_Rt(:,j+1) = abs(V_Rt(:,j)) - kappa*(abs(V_Rt(:,j))-theta)*dt + omega*sqrt(abs(V_Rt(:,j))).*dW1(:,j);
        V_R(:,j+1) = abs(V_Rt(:,j+1));
        LS_R(:,j+1) = log(S_R(:,j))+(r-(1/2).*V_R(:,j))*dt+sqrt(V_R(:,j)).*dW2(:,j);
        S_R(:,j+1)=exp(LS_R(:,j+1));
    end
    
    C_R(k) = exp(-r*T)*(mean(max(S_R(:,end)-K,0)));
end
 
Price_Reflection = mean(C_R)
SE_Reflection = sqrt(mean((C_R-mean(C_R)).^2))
Bias_Reflection = abs(Price_Reflection-x)
RMSE_Reflection = sqrt(Bias_Reflection^2+SE_Reflection^2)
toc;

% Partial Truncation
S_P = zeros(path,T*N+1);
V_Pt = zeros(path,T*N+1);
S_P(:,1) = S;
V_Pt(:,1) = v;
V_P(:,1) = v;
V_P = zeros(path,T*N+1);
C_P = zeros(1,iterations);
LS_P=zeros(path,T*N+1);
 

tic;
for k = 1:iterations
    dW1 = sqrt(dt)*randn(path,T*N);
    Z = randn(path,T*N);
    dW2 = rho*dW1+sqrt(1-rho^2)*Z*sqrt(dt);
    for j = 1:T*N
        V_Pt(:,j+1) = V_Pt(:,j) - kappa*(V_Pt(:,j)-theta)*dt + omega*sqrt(max(V_Pt(:,j),0)).*dW1(:,j);
        V_P(:,j+1) = max(V_Pt(:,j+1),0);
        LS_P(:,j+1) = log(S_P(:,j))+(r-(1/2).*V_P(:,j))*dt+sqrt(V_P(:,j)).*dW2(:,j);
        S_P(:,j+1)=exp(LS_P(:,j+1));
    end
    
    C_P(k) = exp(-r*T)*(mean(max(S_P(:,end)-K,0)));
end
 
Price_Partial = mean(C_P)
SE_Partial = sqrt(mean((C_P-mean(C_P)).^2))
Bias_Partial = abs(Price_Partial-x)
RMSE_Partial = sqrt(Bias_Partial^2+SE_Partial^2)
toc;


%FULL Trun
S_F = zeros(path,T*N+1);
V_Ft = zeros(path,T*N+1);
S_F(:,1) = S;
V_Ft(:,1) = v;
V_F(:,1) = v;
V_F = zeros(path,T*N+1);
C_F = zeros(1,iterations);
LS_F=zeros(path,T*N+1);


tic;
for l = 1:iterations
    dW1 = sqrt(dt)*randn(path,T*N);
    Z = randn(path,T*N);
    dW2 = rho*dW1+sqrt(1-rho^2)*Z*sqrt(dt);
    for j = 1:T*N
        V_Ft(:,j+1) = V_Ft(:,j) - kappa*(max(V_Ft(:,j),0)-theta)*dt + omega*sqrt(max(V_Ft(:,j),0)).*dW1(:,j);
        V_F(:,j+1) = max(V_Ft(:,j+1),0);
        LS_F(:,j+1) = log(S_F(:,j))+(r-(1/2).*V_F(:,j))*dt+sqrt(V_F(:,j)).*dW2(:,j);
        S_F(:,j+1)=exp(LS_F(:,j+1));
    end
    C_F(l) = exp(-r*T)*(mean(max(S_F(:,end)-K,0)));
        
end
 
Price_Full = mean(C_F)
SE_Full = sqrt(mean((C_F-mean(C_F)).^2))
Bias_Full = abs(Price_Full-x)
RMSE_Full = sqrt(Bias_Full^2+SE_Full^2)
toc;

