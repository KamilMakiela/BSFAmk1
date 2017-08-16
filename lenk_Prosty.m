%wywo³anie dla Cobb-Douglasa
%w = lenk_Prosty(s_B(1,100001:1000000),b_B(:,100001:1000000),user.n,user.t,user.priors.ex,user.priors.cov,user.vecs.all(:,1:2),user.model);
%Wywo³anie dla transloga
%w = lenk_Prosty(s_B(1,100001:1000000),b_B(:,100001:1000000),user.n,user.t,user.priors.ex,user.priors.cov,[user.vecs.all(:,:,1)' user.vecs.all(:,:,2)'],user.model);

function korekta = lenk_Prosty(s, b, n, t, be, C, vecs, model)

s = 1./(s.^2);
korekta.z_prec(1,1) = min(s);
korekta.z_prec(1,2) = max(s);

dim = size(b,1);

for i = 1:dim
    korekta.z_b(i,1) = min(b(i,:));
    korekta.z_b(i,2) = max(b(i,:));
end

%DIRECT SAMPLING DLA bet z ucieciem
%Wykombinowaæ ¿eby siê coœ macierzowo liczy³o

switch model
    case {1,2, 3, 4}
        LOSOWAN = 200;
        Lenk_CD;
    case {5,6}
        %translog
        LOSOWAN = 10000;
        Lenk_TR;
    case {7,8}
        %translog t
        LOSOWAN = 20;
        Lenk_TRt;
    case {9,10}
        %translog LT
        LOSOWAN = 20;
        Lenk_TRLT
        
    case {11,12}
        %translog QT
        LOSOWAN = 10;
        Lenk_TRQT
        %korekta.bet(1,1) = -7;
        %sigma = sqrt(diag(C));
        %korekta.bet = log10(prod(normcdf(korekta.z_b(:,2),be,sigma,C) - normcdf(korekta.z_b(:,1),be,sigma,C)));
end

korekta.prec_alt =  log10(gamcdf(korekta.z_prec(1,2), 10^-6,10^6) - gamcdf(korekta.z_prec(1,1), 10^-6,10^6));

korekta.lenk2 = korekta.prec_alt+korekta.bet(1,1);
%korekta.lenk2 = korekta.prec(1,3)*korekta.fi_u(1,3)*korekta.bet(1,3);

