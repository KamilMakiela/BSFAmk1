function wynik = lindley_test(cov,est,vect)

n = size(vect,1);
l = size(est,1);
df = l-n;
for i=1:n
    cov(vect(i),:) = [];
    cov(:,vect(i)) = [];
    est(vect(i)) = [];
    vect = vect-1;
end

chikw = est'*inv(cov)*est;

wynik.cov = cov;
wynik.est = est;
wynik.chikw = chikw;
wynik.pvalue =  1-chi2cdf(chikw,df);