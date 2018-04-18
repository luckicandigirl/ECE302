%% Daniel Nakhimovich & Sara Huang
% ECE302: MATLAB Project 2
clear all; clc
colors = get(gca,'colororder');
close all;

%% (1) (a)
A = 2;
sigmaV2 = 1;
C00=0;
C11=0;
C10=1;
C01=1;
P0=0.8;
P1=0.2;
eta = (C10-C00)/(C01-C11)*(P0/P1);
gamma1 = (A/2) + sigmaV2*log(eta)/A;

s1 = randi(10,1,1000)>8;
y1 = normrnd(A*s1,sqrt(sigmaV2));
rec1 = y1>gamma1;
perr1 = sum(s1~=rec1)/1000
perrt1 = 0.8*(1-normcdf(gamma1,0,sqrt(sigmaV2)))+0.2*normcdf(gamma1,A,sqrt(sigmaV2))

%% (1) (b) (c)
n = 1024;
m = 64;

C00=0;
C11=0;
C10=1;
C01=10;
eta = (C10-C00)/(C01-C11)*(P0/P1);

snrs = [5 4 3 2
       ;4 3 2 1];
      
fs = length(snrs);
figure
for l = 1:fs
    snr = snrs(:,l)';
    N = logspace(-10,10,2048);
    gamma = (snr(1)/2) + snr(2)*log(N)/snr(1);
    ggamma = (snr(1)/2) + snr(2)*log(eta)/snr(1);

    for k = 1:m
        s = randi(10,1,n)>8;
        y = normrnd(snr(1)*s,sqrt(snr(2)),1,n);
        rec2 = y>gamma';
        pdet1(k,:) = sum((s & rec2)')./sum(s');
        pfa1(k,:) = sum((~s & rec2)')./sum((~s)');
        
        %prec = y>ggamma;
        %ppdet(k,:) = sum((s & prec)')./sum(s');
        %ppfa(k,:) = sum((~s & prec)')./sum((~s)');
    end
    pdet1 = mean(pdet1);
    pfa1 = mean(pfa1);
    Pdet1 = qfunc((gamma-snr(1))/sqrt(snr(2)));
    Pfa1 = qfunc(gamma/sqrt(snr(2)));
    ppdet = qfunc((ggamma-snr(1))/sqrt(snr(2)));
    ppfa = qfunc(ggamma/sqrt(snr(2)));
    plot(pfa1,pdet1,'--','color',colors(l,:))
    hold on
    plot(Pfa1,Pdet1,ppfa,ppdet,'o','color',colors(l,:))
end
title('ROC Curves')
xlabel('P_F')
ylabel('P_D')
legend('SNR=5/4 (exp)','SNR=5/4','C01=10*C10','SNR=4/3 (exp)','SNR=4/3','C01=10*C10','SNR=3/2 (exp)','SNR=3/2','C01=10*C10','SNR=2/1 (exp)','SNR=2/1','C01=10*C10','Location','East')

%% (1) (d)
p = linspace(0,1,n);
p = p(1:end-1);
cost1 = ((C01-C00)-(C10-C00)*ppfa-(C01-C11)*ppdet)*p+C00+(C10-C00)*ppfa;
figure
plot(p,cost1,'m')
hold on
title('Expected Cost')
xlabel('p*')
ylabel('E[C]')

%% (1) (e)
Dif = abs(Pdet1-((C01-C00)/(C01-C11)-(C10-C00)/(C01-C11)*Pfa1));
index = find(Dif == min(Dif(:)));
PD = Pdet1(index)
PF = Pfa1(index)
cost2 = ((C01-C00)-(C10-C00)*PF-(C01-C11)*PD)*p+C00+(C10-C00)*PF;
plot(p,cost2)
legend('C01=10*C10','Minimax','Location','East')

figure
plot(Pfa1,Pdet1)
hold on
plot(p,((C01-C00)/(C01-C11)-(C10-C00)/(C01-C11)*p),'m')
scatter(PF,PD)
title('ROC Curve + Minimax')
xlabel('P_F')
ylabel('P_D')
legend('SNR=2/1','Minimax Equation','Location','East')

%% (1) (f)
C00=0;
C11=0;
C10=1;
C01=1;
eta = (C10-C00)/(C01-C11)*(P0/P1);
sigmaV2=1;
sigmazV2=4;
gamma2 = abs(2*(sigmaV2*sigmazV2)/(sigmaV2-sigmazV2)*log(sqrt(sigmaV2/sigmazV2)*eta));

sigma = [sigmazV2 sigmaV2];
s2 = randi(10,1,1000)>8;
y2 = normrnd(0,sqrt(sigma(2-s2)));
rec2 = y2.^2>gamma2;
perr2 = sum(s2~=rec2)/1000
perrt2 = 0.8*(normcdf(sqrt(abs(gamma2)),0,sqrt(sigmazV2))-0.5)+0.2*(normcdf(sqrt(abs(gamma2)),0,sqrt(sigmaV2))-0.5)

C00=0;
C11=0;
C10=1;
C01=10;
eta = (C10-C00)/(C01-C11)*(P0/P1);

snrs = [32 16 8 4
       ;1  1  1 1];
      
fs = length(snrs);
figure
for l = 1:fs
    snr = snrs(:,l)';
    N = logspace(log10(sqrt(snr(1)/snr(2))),-10,2048);
    gamma = abs(2*(snr(2)*snr(1))/(snr(2)-snr(1))*log(sqrt(snr(2)/snr(1))*N));
    
    for k = 1:m
        s = randi(10,1,n)>8;
        y = normrnd(0,sqrt(snr(2-s)),1,n);
        rec3 = y.^2>gamma';
        pdet3(k,:) = sum((s & rec3)')./sum(s');
        pfa3(k,:) = sum((~s & rec3)')./sum((~s)');
    end
    pdet3 = mean(pdet3);
    pfa3 = mean(pfa3);
    Pdet3 = 2*qfunc(sqrt(gamma)/sqrt(snr(1)));
    Pfa3 = 2*qfunc(sqrt(gamma)/sqrt(snr(2)));
    plot(pfa3,pdet3,'--','color',colors(l,:))
    hold on
    plot(Pfa3,Pdet3,'color',colors(l,:))
end
title('ROC Curve for Different Ratios of \sigma_z^2 to \sigma^2')
xlabel('P_F')
ylabel('P_D')
legend('32:1 (exp)','32:1','16:1 (exp)','16:1','8:1 (exp)','8:1','4:1 (exp)','4:1','Location','East')

%% (2)
n = 1024;
m = 64;

lambdas = [1 1 1  1
          ;2 8 32 128];
      
fs = length(lambdas);
figure
for l = 1:fs
    lambda = lambdas(:,l)';
    N = logspace(log10(lambda(2)/lambda(1)),-10,2048);
    gamma = log(N*lambda(1)/lambda(2))/(lambda(1)-lambda(2));

    for k = 1:m
        sent = randi([1 2],1,n);
        obs = exprnd(1./lambda(sent));
        rec = obs<gamma';
        pdet(k,:) = sum((sent==2 & rec)')./sum((sent==2)');
        pfa(k,:) = sum((sent==1 & rec)')./sum((sent==1)');
    end
    pdet = mean(pdet);
    pfa = mean(pfa);
    Pdet = expcdf(gamma,1./lambda(2));
    Pfa = expcdf(gamma,1./lambda(1));
    plot(pfa,pdet,'--','color',colors(l,:))
    hold on
    plot(Pfa,Pdet,'color',colors(l,:))
end
title('ROC Curve for Different Ratios of Rates')
xlabel('P_F')
ylabel('P_D')
legend('1:2 (exp)','1:2','1:8 (exp)','1:8','1:32 (exp)','1:32','1:128 (exp)','1:128','Location','East')

%% (3)
load('Iris.mat')

shuf = randperm(length(labels));
features = features(shuf,:);
labels = labels(shuf);

trainFeatures = features(1:2:end,:);
trainLabels = labels(1:2:end);
testFeatures = features(2:2:end,:);
testLabels = labels(2:2:end);
testLabels = [1:length(testLabels);testLabels';ones(1,length(testLabels))]';
testLabels = full(spconvert(testLabels));

ulabs = unique(labels);
for k = 1:length(ulabs)
    MU(k,:) = mean(trainFeatures(trainLabels==k,:));
    SIGMA(:,:,k) = cov(trainFeatures(trainLabels==k,:));
    likelihoods(:,k) = mvnpdf(testFeatures,MU(k,:),SIGMA(:,:,k));
end
[perr,cmat,~,~] = confusion(testLabels',likelihoods')
