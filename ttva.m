% ttva.m

%% inputs
soa = 250;
w = [.75 .25];
targetDur = 30;

%% params
tau = targetDur;
muLock = 80;
muDwell = 400;
muDecay = 200;
C = 40;
tauStar = [NaN tau+muDecay]; % leave tauStar(1) empty for now

lambdaLock = 1/muLock;

%% calculate tauStar
if soa <= muDecay
    tauStar(1) = tau + soa;
else
    tauStar(1) = tau + muDecay;
end
    


%% vT2
pFree = exp(C*soa) + C/(lambdaLock-C)*(exp(C*soa)-exp(lambdaLock*soa));

v(2) = C*(w(2)/sum(w).*pFree + pReleased);

%% just approximating based on Petersen 2012 Fig 6
C = 1;
w = [.5 .5];
x = 0:800;
tau = 2;
pFree = exp(-.03*x);
pReleased = 1-exp(-.0015*x);

vT2 = C*(w(2)/sum(w).*pFree + pReleased);
vT2_valid = C*(.75/sum(w).*pFree + pReleased);
vT2_invalid = C*(.25/sum(w).*pFree + pReleased);

pT2 = 1-exp(-tau*vT2);
pT2_valid = 1-exp(-tau*vT2_valid);
pT2_invalid = 1-exp(-tau*vT2_invalid);

figure
hold on
plot(x, [pFree' pReleased' pT2'])
plot(x, pT2_valid, 'r--')
plot(x, pT2_invalid, 'r-.')
legend('p(free)','p(released)','p(T2)','p(T2) valid','p(T2) invalid')


    