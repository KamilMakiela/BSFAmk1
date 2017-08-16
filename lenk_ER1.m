%wywo³anie dla Cobb-Douglasa
%w = lenk_ER1(s_B(1,poczatek:koniec),fi_B(1,poczatek:koniec),u_B(:,poczatek:koniec),b_B(:,poczatek:koniec),user.n,user.t,user.med_apr,user.priors.ex,user.priors.cov,user.vecs.all(:,1:2), user.model);
%Wywya³anie dla transloga
%w = lenk_ER1(s_B(1,poczatek:koniec),fi_B(1,poczatek:koniec),u_B(:,poczatek:koniec),b_B(:,poczatek:koniec),user.n,user.t,user.med_apr,user.priors.ex,user.priors.cov,[user.vecs.all(:,:,1)' user.vecs.all(:,:,2)'], user.model);

function korekta = lenk_ER1(s, fi, u, b, n, t, r0, be,C,vecs, model)
s = 1./(s.^2);
korekta.z_prec(1,1) = min(s);
korekta.z_prec(1,2) = max(s);

korekta.z_fi(1,1) = min(fi);
korekta.z_fi(1,2) = max(fi);

for i = 1:n*t
    korekta.z_u(i,1) = min(u(i,:));
    korekta.z_u(i,2) = max(u(i,:));
end
dim = size(b,1);

for i = 1:dim
    korekta.z_b(i,1) = min(b(i,:));
    korekta.z_b(i,2) = max(b(i,:));
end

%DIRECT SAMPLING DLA bet z ucieciem
%Wykombinowaæ ¿eby siê coœ macierzowo liczy³o


switch model
    case {1,3,4}
        LOSOWAN = 2000;
        Lenk_CD;
    case {5,6}
        %translog
        LOSOWAN = 1000;
        Lenk_TR;
    case {2, 7,8}
        %translog t
        LOSOWAN = 40;
        Lenk_TRt;
    case {9,10}
        %translog LT
        LOSOWAN = 1000;
        Lenk_TRLT
        
    case {11,12}
        %translog QT
        LOSOWAN = 1000;
        Lenk_TRQT
        %korekta.bet(1,1) = -7;
        %sigma = sqrt(diag(C));
        %korekta.bet = log10(prod(normcdf(korekta.z_b(:,2),be,sigma,C) - normcdf(korekta.z_b(:,1),be,sigma,C)));
end

% LOSOWAN = 1000;
% korekta.bet = zeros(1,3);
% h = waitbar(0,'bety...');
% while korekta.bet(1,2) < LOSOWAN
%     %for i = 1:LOSOWAN
%     draw = be + (chol(C)')*randn(dim,1);
%     %disp(size(draw)); disp(size(vecs));
%     elast = draw'*vecs;
%     if isempty(find(elast<0, 1))
%         korekta.bet(1,1) = korekta.bet(1,1)+1;
%         if size(draw(korekta.z_b(:,1) <= draw & draw <= korekta.z_b(:,2)),1) == dim
%             korekta.bet(1,2) = korekta.bet(1,2) + 1;
%             waitbar(korekta.bet(1,2)/LOSOWAN,h);
%         end
%     end
%     
% end
% korekta.bet(1,3) = korekta.bet(1,2)./korekta.bet(1,1);
% close(h);

%DIRECT SAMPLING DLA FI oraz u
LOSOWAN = 1;
LICZBA_FI = 100000;
h = waitbar(0,'fi oraz uit...');
korekta.fi_u = zeros(1,3);
while korekta.fi_u(1,2) < LOSOWAN
    draw = randg(1,[LICZBA_FI,1])/(-log(r0));
    korekta.fi_u(1,1) = korekta.fi_u(1,1)+LICZBA_FI;
    trafione_fi = draw(korekta.z_fi(1,1) <= draw & draw <= korekta.z_fi(1,2));
    l_trafien_fi = size(trafione_fi,1);
    %disp(l_trafien_fi);
    if l_trafien_fi > 0
        %for a = 1:l_trafien_fi
        %    draw = randg(1,[n*t,1])./trafione_fi(a);
        %    rozm = size(draw(korekta.z_u(:,1) <= draw & draw <= korekta.z_u(:,2)),1);
        %    if rozm == n*t
        %        korekta.fi_u(1,2) = korekta.fi_u(1,2) + 1;
        %        %disp(rozm);
        %    end
        %    waitbar(a/l_trafien_fi,h);
        %end
        korekta.fi_u(1,2) = 1;
        for i = 1:size(trafione_fi,1)
            cos(:,i) = gamcdf(korekta.z_u(:,2), 1, 1./trafione_fi(i)) - gamcdf(korekta.z_u(:,1), 1, 1./trafione_fi(i));
            waitbar(i/l_trafien_fi,h);
        end
        x = prod(cos,1);
        y = mean(x);
        korekta.fi_u_alt = log10(y*l_trafien_fi/LICZBA_FI); % iloczyn prawdopodobienstw p(fi) oraz p(u wylosowanych na podstawie przyjetych fi)
    end
end
korekta.fi_u(1,3) = korekta.fi_u(1,2)/korekta.fi_u(1,1);
clear draw;
close(h);


%DIRECT SAMPLING DLA PRECYZJI

% LOSOWAN = 5;
% LICZBA_PREC = 10000000;
% a=0;
% korekta.prec = zeros(1,3);
% h = waitbar(0,'precyzja...');
% while korekta.prec(1,2) < LOSOWAN
%     a = a+1;
%     draw = randg(0.5*10^(-6),[LICZBA_PREC,1])/(0.5*10^(-6));
%     korekta.prec(1,1) = korekta.prec(1,1)+LICZBA_PREC;
%     korekta.prec(1,2) = korekta.prec(1,2) + size(draw(korekta.z_prec(1,1) <= draw & draw <= korekta.z_prec(1,2)),1);
%     waitbar(a/LOSOWAN,h);
% end
% close(h)
% korekta.prec(1,3) = korekta.prec(1,2)/korekta.prec(1,1);
korekta.prec_alt =  log10(gamcdf(korekta.z_prec(1,2), 10^-6,10^6) - gamcdf(korekta.z_prec(1,1), 10^-6,10^6));

% DIRECT SAMPLING DLA u
% h = waitbar(0,'uit...');
% for i = 1:n*t
%     korekta.u(i,1:3) = zeros(1,3);
%     fi_draws = randg(1,[10000,1])/(-log(r0));
%     while korekta.u(i,2) < 10
%         draw = randg(1,[10000,1])./(fi_draws);
%         korekta.u(i,1) = korekta.u(i,1)+10000;
%         korekta.u(i,2) = korekta.u(i,2) + size(draw(korekta.z_u(i,1) <= draw & draw <= korekta.z_u(i,2)),1);
%     end
%     waitbar(i/10,h)
%     korekta.u(i,3) = korekta.u(i,2)/korekta.u(i,1);
% end
% close(h);

%korekta.bet_alt = log10(korekta.bet(1,3));
korekta.bet_alt = korekta.bet(1,1);
%korekta.lenk1 = log10(korekta.prec(1,3))+log10(korekta.fi_u(1,3))+log10(korekta.bet(1,3));
korekta.lenk2 = korekta.prec_alt+korekta.fi_u_alt+korekta.bet_alt;

%korekta.lenk2 = korekta.prec(1,3)*korekta.fi_u(1,3)*korekta.bet(1,3);

