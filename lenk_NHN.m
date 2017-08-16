%wywo³anie dla Cobb-Douglasa
%w = lenk_NHN(s_B(1,poczatek:koniec),fi_B(1,poczatek:koniec),u_B(:,poczatek:koniec),b_B(:,poczatek:koniec),user.n,user.t,user.med_apr,user.priors.ex,user.priors.cov,user.vecs.all(:,1:2),user.model);
%Wywo³anie dla transloga
%w = lenk_NHN(s_B(1,poczatek:koniec),fi_B(1,poczatek:koniec),u_B(:,poczatek:koniec),b_B(:,poczatek:koniec),user.n,user.t,user.med_apr,user.priors.ex,user.priors.cov,[user.vecs.all(:,:,1)' user.vecs.all(:,:,2)'],user.model);

function korekta = lenk_NHN(s, fi, u, b, n, t, r0, be,C,vecs, model)

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
    case {1,3, 4}
        LOSOWAN = 1000;
        Lenk_CD;
    case {5,6}
        %translog
        LOSOWAN = 10000;
        Lenk_TR;
    case {2,7,8}
        %translog t
        LOSOWAN = 1000;
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

%Ewaluacja omegi na podstawie ró¿nicy wartoœci dystrybuant w dwóch punktach
korekta.fi_alt =  log10(gamcdf(korekta.z_fi(1,2), 1, 1/(10*((log(r0)^2)))) - gamcdf(korekta.z_fi(1,1),  1, 1/(10*((log(r0)^2)))));
%disp(korekta.fi_alt);

if korekta.fi_alt > -7
    %Direct sampling fi w celu ewaluacji ewaluacji p(u)
    LOSOWAN = 1;
    LICZBA_FI = 2000*round(1/power(10,korekta.fi_alt));
    h = waitbar(0,'fi oraz uit...');
    korekta.fi_u = zeros(1,3);
    while korekta.fi_u(1,2) < LOSOWAN
        draw = randg(5,[LICZBA_FI,1])/(10*((log(r0)^2)));
        korekta.fi_u(1,1) = korekta.fi_u(1,1)+LICZBA_FI;
        trafione_fi = draw(korekta.z_fi(1,1) <= draw & draw <= korekta.z_fi(1,2));
        l_trafien_fi = size(trafione_fi,1);
        if l_trafien_fi > 0
            %disp(l_trafien_fi);
            for i = 1:size(trafione_fi,1)
                korekta.u.wynik(:,i) = 2*(normcdf(korekta.z_u(:,2), 0, 1./trafione_fi(i)) - normcdf(korekta.z_u(:,1), 0, 1./trafione_fi(i)));
                korekta.u.log10prob(i) = sum(log10(korekta.u.wynik(:,i)),1);
                korekta.u.prob(i) = power(10,korekta.u.log10prob(i));
                waitbar(i/l_trafien_fi,h)
            end
            korekta.u_alt = log10(mean(korekta.u.prob));
            korekta.fi_u(1,2) = korekta.fi_u(1,2) + 1;
        end
    end
    %korekta.fi_u(1,3) = korekta.fi_u(1,2)/korekta.fi_u(1,1);
    clear draw;
    close(h)
else
    %min
    korekta.u.wynik(:,1) = 2*(normcdf(korekta.z_u(:,2), 0, 1./korekta.z_fi(1,2)) - normcdf(korekta.z_u(:,1), 0, 1./korekta.z_fi(1,2)));
    korekta.u.log10prob(1) = sum(log10(korekta.u.wynik(:,1)),1);
    korekta.u.prob(1) = prod(korekta.u.wynik(:,1),1);
    %max
    korekta.u.wynik(:,2) = 2*(normcdf(korekta.z_u(:,2), 0, 1./korekta.z_fi(1,1)) - normcdf(korekta.z_u(:,1), 0, 1./korekta.z_fi(1,1)));
    korekta.u.log10prob(2) = sum(log10(korekta.u.wynik(:,2)),1);
    korekta.u.prob(1) = prod(korekta.u.wynik(:,2),1);
    %mean
    srodek = 0.5*(korekta.z_fi(1,1)+korekta.z_fi(1,2));
    korekta.u.wynik(:,3) = 2*(normcdf(korekta.z_u(:,2), 0, 1./srodek) - normcdf(korekta.z_u(:,1), 0, 1./srodek));
    korekta.u_alt = sum(log10(korekta.u.wynik(:,3)),1);
    %blad
    korekta.u.odchyl = 0.5*(korekta.u.log10prob(2)-korekta.u.log10prob(1));
end
% %----STARE----
% %DIRECT SAMPLING DLA FI oraz u
% LOSOWAN = 1;
% LICZBA_FI = 1000000;
% h = waitbar(0,'fi oraz uit...');
% korekta.fi_u = zeros(1,3);
% while korekta.fi_u(1,2) < LOSOWAN
%     draw = randg(5,[LICZBA_FI,1])/(10*((log(r0)^2)));
%     korekta.fi_u(1,1) = korekta.fi_u(1,1)+LICZBA_FI;
%     trafione_fi = draw(korekta.z_fi(1,1) <= draw & draw <= korekta.z_fi(1,2));
%     l_trafien_fi = size(trafione_fi,1);
%     if l_trafien_fi > 0
%         %disp(l_trafien_fi);
%         for a = 1:l_trafien_fi
%             %draw = randn(n*t,1)*sqrt(1/trafione_fi(a)); 
%             draw = kmdraw(0,sqrt(1/trafione_fi(a)),n*t);
%             %draw = randn([n*t,1])./trafione_fi(a);
%             rozm = size(draw(korekta.z_u(:,1) <= draw & draw <= korekta.z_u(:,2)),1);
%             if rozm == n*t
%                 korekta.fi_u(1,2) = korekta.fi_u(1,2) + 1;
%                 disp('SUKCES');
%             end
%             waitbar(a/l_trafien_fi,h)
%         end
%         
%     end
% end
% korekta.fi_u(1,3) = korekta.fi_u(1,2)/korekta.fi_u(1,1);
% clear draw;
% close(h);


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


%korekta.lenk1 = log10(korekta.prec(1,3))+log10(korekta.fi_u(1,3))+log10(korekta.bet(1,3));
%korekta.lenk2 = log10(korekta.prec_alt)+log10(korekta.fi_u(1,3))+log10(korekta.bet(1,3));
korekta.lenk2 = korekta.prec_alt+korekta.fi_alt+korekta.u_alt+korekta.bet(1,1);
%korekta.lenk2 = korekta.prec(1,3)*korekta.fi_u(1,3)*korekta.bet(1,3);

