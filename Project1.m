%% Sara Huang & Daniel Nakhimovich
% ECE302: MATLAB Project 1
clear all; close all; clc

%% 1 (a) (b) 

h = 0.5;
mu = 2;
sigma0 = 0.5;
sigma2 = 1;

msemap = @(x,n) (h*sigma0^2*sum(x)+mu*sigma2^2)./(h^2*[1:n]*sigma0^2+sigma2^2);

ml = @(x,n) sum(x)./[1:n];

%% 1 (c)

n = 222;
m = 666;

pmserrs = zeros(m,n);
pmlerrs = zeros(m,n);

for k = 1:m
    theta = triu(normrnd(mu,sigma0,n));
    x = triu(h*theta+normrnd(0,sigma2,n));

    estmse = triu(msemap(x,n));
    mserrs = sum((theta-estmse).^2);
    pmserrs(k,:) = mserrs./[1:n];

    estml = triu(ml(x,n));
    mlerrs = sum((theta-estml).^2);
    pmlerrs(k,:) = mlerrs./[1:n];
end

pmserrs = mean(pmserrs);
pmlerrs = mean(pmlerrs);

plot([1:n],pmserrs,'m')
hold on
plot([1:n],pmlerrs)
axis([0 n 0 100])
xlabel('Number of Measurements Taken')
ylabel('Square Error')
title('Minimum Mean-Squared Error of \theta estimate')
[~,~,~,leg] = legend('MMSE (normal)','ML');

%% 1 (d)

pmserrs2 = zeros(m,n);axis([0 n 0 100])

pmserrs3 = zeros(m,n);

for k = 1:m
    theta2 = triu(unifrnd(-5,5,n));
    x = triu(h*theta2+normrnd(0,sigma2,n));

    estmse2 = triu(msemap(x,n));
    mserrs2 = sum((theta2-estmse2).^2);
    pmserrs2(k,:) = mserrs2./[1:n];

    theta3 = triu(exprnd(mu,n));
    x = triu(h*theta3+normrnd(0,sigma2,n));

    estmse3 = triu(msemap(x,n));
    mserrs3 = sum((theta3-estmse3).^2);
    pmserrs3(k,:) = mserrs3./[1:n];
end

pmserrs2 = mean(pmserrs2);
pmserrs3 = mean(pmserrs3);

plot([1:n],pmserrs2)
hold on
plot([1:n],pmserrs3)
axis([0 n 0 100])
xlabel('Number of Measurements Taken')
ylabel('Square Error')
title('Minimum Mean-Squared Error of \theta estimate')
if length(size(leg)) == 2
    legend(leg{1},leg{2},'MMSE (uniform)', 'MMSE (exponential)')
else
    legend('MMSE (uniform)', 'MMSE (exponential)')
end

%% 2 (a)

% Mathematical Model:
fx = @(x,t,t1,t2,sigma) ...
     (t<t1).*(normpdf(x,-1,sigma)/2+normpdf(x,1,sigma)/2) + ...
     (t>=t1).*(t<=t2).*(normpdf(x,-2,sigma)/4+normpdf(x,0,sigma)/2+normpdf(x,2,sigma)/4) + ...
     (t>t2).*(normpdf(x,-1,sigma)/2+normpdf(x,1,sigma)/2);
     
% Likelihood Function:
L = @(x,t,t1,t2,sigma) sum(log(fx(x,t,t1,t2,sigma)));

%% 2 (b)

tt1 = 10;
tt2 = 80;
tt3 = 100;
sigma = 0.11;

vals = [-1 1];
sig1 = vals(randi([1,2],1,tt1));
interf = vals(randi([1,2],1,tt2-tt1))+vals(randi([1,2],1,tt2-tt1));
sig2 = vals(randi([1,2],1,tt3-tt2));

X = [sig1,interf,sig2];
Xobs = X+sigma*randn(1,tt3);
Xind = [1:tt3];

max = -Inf;
pt1 = 0;
pt2 = 0;
for k1 = [1:tt3]
   for k2 = [k1:tt3]
       m = L(Xobs,Xind,k1,k2,sigma);
       if m > max
           max = m;
           pt1 = k1-1;
           pt2 = k2;
       end
   end
end
pt1
pt2

%% 2 (c)
close all;
figure
plot(Xind,Xobs,'m')
hold on
plot([pt1,pt1],[-3.5,-2.5],'k')
plot([pt2,pt2],[-3.5,-2.5],'k')
plot([pt1,pt2],[-3,-3],'k')
title('Received Signal')
xlabel('time')
ylabel('Voltage (V)')
legend('Signal','+Interference')
