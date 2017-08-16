% Script to generate sample for empirical example from the manual
n = 25;
T = 10;
nT = n*T;

ref_sig_u = 0.2;
ref_sig_v = 0.1;
beta_0 = 2;
beta_1 = 0.5;
beta_2 = 0.5;

x = randn(nT,1);
lnx1 = (x -ones(nT,1)*mean(x,1))./(ones(nT,1)*std(x,1,1));

x = randn(nT,1);
lnx2 = (x -ones(nT,1)*mean(x,1))./(ones(nT,1)*std(x,1,1));

ref_u = randg(1,nT,1).*ref_sig_u; %inefficiencies from exponential distribution

ref_v = normrnd(0,ref_sig_v,[nT 1]);
ref_v = ref_sig_v.*((ref_v-mean(ref_v))./std(ref_v,1)); %disturbances from normal distribution

lny = beta_0.*ones(nT,1) + beta_1.*lnx1 + beta_2.*lnx2 + ref_v - ref_u;

data_y  = exp(reshape(lny,25,10));
data_x1 = exp(reshape(lnx1,25,10));
data_x2 = exp(reshape(lnx2,25,10));