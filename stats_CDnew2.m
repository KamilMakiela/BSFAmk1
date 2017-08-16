

function [beta, elastycznosci, skala, eff, lamb, SS, mdd, ile_danych] = stats_CDnew2(folder, burnin, hand, vectors, pliki, n, t)

% Nowa (oszczedniejsza) funkcja do liczenia statystyk po
if ~isdir(folder)
    msgbox('There is no such results folder','Getting simulations data','error');
    return
end

l_plikow = size(pliki,2);

%TWORZE WAITBAR
h = waitbar(0,'Marginal data density',...
    'Name','Calculating...');
assignin('base', 'Waitbar_handle2', h);

%----------------petla pobierajaca PDF - skladniki do liczenie gestosci
%aposteriori
var(1:l_plikow) = struct('mdd_d', 0);
for b = 1:l_plikow
    var(b) = load(pliki(b).mdd,'mdd_d');
    
end
%assignin('base','wwr',var);
%tworze wektor losowan pdf
disp(size(var(1).mdd_d));
mddd = [var(:).mdd_d]';
disp(size(mddd));
clear var;
rlanc = size(mddd,1); %TA ZMIENNA SI? PRZYDA!!
%save_mdd = strcat(folder,'\lg10mdd');
%disp(save_mdd);
%save(save_mdd,'mddd');
mddd = mddd((burnin+1):rlanc,:);
mdd = rn2(mddd);

%mdd = 1/mean(exp(mddd((burnin+1):rlanc,:)));
%n = hist(zmienna(burnin:rlanc,1),40);
%plot(n);
%h = figure;
%hist(INEF(burnin:rlanc,1),40);
%czyszcze pamiec
%pomoc = rlanc - burnin+1;
%xxhat = (cumsum(mddd(burnin:rlanc,1).^2)')./(1:pomoc);
%plot(hand, xxhat);
clear mddd;

%---- NOWE DLA SIGMA
counter = 2/8;
msg = sprintf('Standard deviation of the symmetric error');
waitbar(counter,h,msg);

%--- Ta czesc kodu jest potrzebna do pozniejszych obliczen
pomoc = rlanc - burnin;
ile_danych = zeros(l_plikow,2);
var(1:l_plikow) = struct('s_B', 0);
%srednia = zeros(l_plikow,1);
for b = 0:(l_plikow-1)
    c = l_plikow-b;
    var(c) = load(pliki(c).s,'s_B');
    ile_danych(c) = size(var(c).s_B,2);
    if sum(ile_danych) <= pomoc
        %srednia(c) = mean(var(c).mdd_d);
    else
        index = pomoc - sum(ile_danych(c+1:l_plikow));
        if index > 0
            ile_danych(c) = index;
            %disp('tu');
        else
            ile_danych(c) = 0;
            %disp('ale');
        end
    end
end
ile_danych(:,2) = ile_danych(:,1)./sum(ile_danych(:,1));
% Koniec tej czesci kodu


SS = zeros(1,2);
flaga = 0;
%var(1:l_plikow) = struct('s_B', 0);
for b = 1:l_plikow
    if ile_danych(b,1) ~= 0
        %var(b) = load(pliki(b).s,'s_B');
        rozmiar = size(var(b).s_B,2);
        zacznij = rozmiar - ile_danych(b,1)+1;
        if flaga == 0
            param = size(var(b).s_B,1);
            SS = zeros(param,2);
            flaga = 1;
        end
        for a = 1:param
            SS(a,1) = SS(a,1) + mean(var(b).s_B(a,zacznij:rozmiar))*ile_danych(b,2);
            SS(a,2) = SS(a,2) + std(var(b).s_B(a,zacznij:rozmiar))*ile_danych(b,2);
        end
    end
end
clear var;
flaga = 0;


%----------------petla pobierajaca fi
counter = 3/8;
msg = sprintf('Standard deviation of inefficiency');
waitbar(counter,h,msg);
try
    lamb = zeros(1,2);
    var = struct('fi_B', 0);
    for b = 1:l_plikow
        if ile_danych(b,1) ~= 0
            var = load(pliki(b).fi,'fi_B');
            rozmiar = size(var.fi_B,2);
            zacznij = rozmiar - ile_danych(b,1)+1;
            if flaga == 0
                param = size(var.fi_B,1);
                lamb = zeros(param,2);
                flaga = 1;
            end
            for a = 1:param
                lamb(a,1) = lamb(a,1) + mean(1./var.fi_B(a,zacznij:rozmiar))*ile_danych(b,2);
                lamb(a,2) = lamb(a,2) + std(1./var.fi_B(a,zacznij:rozmiar))*ile_danych(b,2);
            end
        end
    end
    clear var;
    flaga = 0;
catch exeption
    disp('There is no inefficiency');
    lamb(1,1) = 0;
    lamb(1,2) =0;
end

%-------BETY, ELASTYCZNOSCI I EFEKT SKALI
counter = 4/8;
msg = sprintf('Elasticities and scale effect');
waitbar(counter,h,msg);

zmiennych = size(vectors.all,3);
obserwacji = size(vectors.all,1);
var = struct('b_B', 0);
if zmiennych == 1
    disp('Executing code for Cobb-Douglas');
    beta = zeros(1,2);
    for b = 1:l_plikow
        %disp('HALO petla');
        if ile_danych(b,1) ~= 0
            var = load(pliki(b).bety,'b_B');
            %disp(size(vectors.all));
            %disp(size(var.b_B));
            r = size(var.b_B,1);
            ELAST = (vectors.all')*var.b_B;
            %ELAST = var.b_B';
            %wyznaczam koniec i poczatek sczytywania danych
            rozmiar = size(var.b_B,2);
            zacznij = rozmiar - ile_danych(b,1)+1;
            %disp('HALO if');
            if flaga == 0
                param = size(ELAST,1);
                beta = zeros(param,2);
                flaga =1;
            end
            for a =1:param
                beta(a,1) = beta(a,1) + mean(ELAST(a,zacznij:rozmiar))*ile_danych(b,2);
                beta(a,2) = beta(a,2) + std(ELAST(a,zacznij:rozmiar))*ile_danych(b,2);
            end
        end
    end
    skala = beta(param,:);
    %elastycznosci = beta(1:(param-2),:);
    elastycznosci = 0;
    beta(param,:) = [];
    clear var;
else
    disp('Executing code for translog');
    beta = zeros(1,2);
    elastycznosci.all = zeros(obserwacji,2,zmiennych); %'2' - bo licze tylko srednia i odchylenie
    elastycznosci.obj = zeros(n,2,zmiennych);
    elastycznosci.time = zeros(t,2,zmiennych);
    for b = 1:l_plikow
        %disp('HALO petla');
        if ile_danych(b,1) ~= 0
            var = load(pliki(b).bety,'b_B');
            %disp(size(var.b_B));
            r = size(var.b_B,1);
            %disp(size(vectors));
            %ELAST = (vectors')*var.b_B;
            %wyznaczam koniec i poczatek sczytywania danych
            rozmiar = size(var.b_B,2);
            zacznij = rozmiar - ile_danych(b,1)+1;
            %disp('HALO if');
            if flaga == 0
                param = size(var.b_B,1);
                beta = zeros(param,2);
                flaga =1;
            end
            for a =1:param
                beta(a,1) = beta(a,1) + mean(var.b_B(a,zacznij:rozmiar))*ile_danych(b,2);
                beta(a,2) = beta(a,2) + std(var.b_B(a,zacznij:rozmiar))*ile_danych(b,2);
            end
            %wyznaczam wszystkie elastycznosci
            for a = 1:zmiennych
                ELAST = ((vectors.all(:,:,a))*(var.b_B(:,zacznij:rozmiar)))';
                %disp(size(ELAST));
                %disp(size(elastycznosci(:,1,a)));
                elastycznosci.all(:,1,a) = elastycznosci.all(:,1,a) + mean(ELAST,1)'*ile_danych(b,2);
                elastycznosci.all(:,2,a) = elastycznosci.all(:,2,a) + std(ELAST,1)'*ile_danych(b,2);
                %save(strcat(folder,'\',num2str(a),'Elast'),'ELAST');
            end
            %wyznaczam elastycznosci srednie po obserwacjach
            for a = 1:zmiennych
                ELAST = ((vectors.obj(:,:,a))*(var.b_B(:,zacznij:rozmiar)))';
                %disp(size(ELAST));
                %disp(size(elastycznosci(:,1,a)));
                elastycznosci.obj(:,1,a) = elastycznosci.obj(:,1,a) + mean(ELAST,1)'*ile_danych(b,2);
                elastycznosci.obj(:,2,a) = elastycznosci.obj(:,2,a) + std(ELAST,1)'*ile_danych(b,2);
                %save(strcat(folder,'\',num2str(a),'Elast'),'ELAST');
            end
            for a = 1:zmiennych
                ELAST = ((vectors.time(:,:,a))*(var.b_B(:,zacznij:rozmiar)))';
                %disp(size(ELAST));
                %disp(size(elastycznosci(:,1,a)));
                elastycznosci.time(:,1,a) = elastycznosci.time(:,1,a) + mean(ELAST,1)'*ile_danych(b,2);
                elastycznosci.time(:,2,a) = elastycznosci.time(:,2,a) + std(ELAST,1)'*ile_danych(b,2);
                %save(strcat(folder,'\',num2str(a),'Elast'),'ELAST');
            end
            
        end
    end
    clear var ELAST;
    skala.all = elastycznosci.all(:,:,zmiennych);
    skala.obj = elastycznosci.obj(:,:,zmiennych);
    skala.time = elastycznosci.time(:,:,zmiennych);
    elastycznosci.all(:,:,zmiennych) = [];
    elastycznosci.obj(:,:,zmiennych) = [];
    elastycznosci.time(:,:,zmiennych) = [];
end
flaga =0;



%------Wydajnosci (+odchylenia), nieefektywnoesci (+odchylenia)
%--Nie zapisuje zmiennej "INEF" do pliku - nie ma takiej potrzeby
counter = 6/8;
msg = sprintf('Efficiency scores');
waitbar(counter,h,msg);

try
    var = struct('u_B', 0);
    eff.all = zeros(1,4);
    for b = 1:l_plikow
        %disp('HALO petla');
        if ile_danych(b,1) ~= 0
            var = load(pliki(b).u,'u_B');
            %disp(size(var.b_B));
            %disp(size(vectors));
            %wyznaczam koniec i poczatek sczytywania danych
            rozmiar = size(var.u_B,2);
            zacznij = rozmiar - ile_danych(b,1)+1;
            %disp('HALO if');
            if flaga == 0
                param = size(var.u_B,1);
                eff.all = zeros(param,4);
                eff.time = zeros(t,2);
                eff.obj = zeros(n,2);
                flaga =1;
            end
            %disp(size(eff));
            %disp(size(mean(exp(-var.u_B(:,zacznij:rozmiar)),2)));
            %disp(param);
            eff.all(:,1) = eff.all(:,1) + mean(exp(-var.u_B(:,zacznij:rozmiar)),2).*ile_danych(b,2);
            eff.all(:,2) = eff.all(:,2) + std(exp(-var.u_B(:,zacznij:rozmiar)'))'.*ile_danych(b,2);
            eff.all(:,3) = eff.all(:,3) + mean(var.u_B(:,zacznij:rozmiar),2).*ile_danych(b,2);
            eff.all(:,4) = eff.all(:,4) + (std(var.u_B(:,zacznij:rozmiar)').*ile_danych(b,2))';
            g1 = 1;
            for g = 1:n:(n*t)
                pomoc = mean(exp(-var.u_B(g:(g+n-1),zacznij:rozmiar)),2).*ile_danych(b,2);
                pomoc2 = std(exp(-var.u_B(g:(g+n-1),zacznij:rozmiar))')'.*ile_danych(b,2);
                %disp(size(pomoc));
                %disp(size(pomoc2));
                eff.obj(:,1) = eff.obj(:,1) + pomoc./t;
                eff.obj(:,2) = eff.obj(:,2) + pomoc2./t;
                eff.time(g1,1) = eff.time(g1,1)+ mean(pomoc,1);
                eff.time(g1,2) = eff.time(g1,2)+ mean(pomoc2,1);
                g1 = g1 +1;
            end
        end
    end
    clear var pomoc pomoc2;
    
catch exception
    disp(exception);
    disp('It seems there are no inefficiencies in this model');
    %INEF = 0;
    eff = 0;
end


%Wyznaczam wykresy
%var(1) = struct('b_B', 0);
try
    var(1:l_plikow) = struct('b_B', 0);
    for b = 1:l_plikow
        var(b) = load(pliki(b).bety,'b_B');
        %SE(b,:) = var.b_B(r,:);
    end
    %assignin('base','zmienna',var);
    %tworze wektor losowan s
    SE = [var(:).b_B]';
    SE = SE(:,r);
    disp(size(SE));
    %SE = SE(:);
    %disp(size(SE));
    clear var;
    save(strcat(folder,'\', 'sequential'),'SE');
    %tu licze interesujace mnie statystyki i generuje wykresy
    SSS(:,1) = mean(SE((burnin+1):rlanc,:));
    SSS(:,2) = std(SE((burnin+1):rlanc,:));
    %try
    pomoc = rlanc - burnin;
    pomoc2 = (1:pomoc);
    %wyznaczam cusumy
    kusum = cumsum(SE((burnin+1):rlanc,1)')./pomoc2;
    clear SE;
    %standaryzuje
    CUSUM = (kusum-SSS(1,1))./SSS(1,2);
    clear kusum;
    %wyznaczam benchmark path
    %bch_var = SSS(1,2)*randn(1,pomoc)+SSS(1,1);
    bch_var = randn(1,pomoc);
    bch_var = (bch_var-mean(bch_var))./std(bch_var);
    bch_var = SSS(1,2).*bch_var+SSS(1,1);
    bch_kusum = cumsum(bch_var(1,1:pomoc))./pomoc2;
    %SSS(1,1) = mean(bch_var);
    %SSS(1,2) = std(bch_var);
    bch_pth = (bch_kusum-SSS(1,1))./SSS(1,2);
    %rysuje wykres
    plot(hand, pomoc2, CUSUM, '-', pomoc2, bch_pth, ':');
    save(strcat(folder,'\', 'cusum'),'CUSUM');
    a = min(CUSUM(1,round(pomoc/8):pomoc));
    a1 = min(bch_pth(1,round(pomoc/8):pomoc));
    if a1<a
        a = a1;
    end
    a = 2*a;
    b = max(CUSUM(1,round(pomoc/8):pomoc));
    b1 = max(bch_pth(1,round(pomoc/8):pomoc));
    if b1>b
        b = b1;
    end
    %a = -0.5;
    %b = 0.5;
    b = 2*b;
    ylim(hand, [a b]);
    xlim(hand, [0, pomoc]);
    title(hand, 'CUSUM path plot for the intercept');
    legend(hand, 'cusum path', 'benchmark path');
    legend(hand, 'boxoff');
    clear CUSUM sred_skum pomoc pomoc2 a b a1 b1;
catch e
    disp(e)
end
delete(h);


end

