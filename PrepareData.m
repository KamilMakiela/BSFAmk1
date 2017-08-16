%Tutaj przygotowuje odpowiednia postac macierzy x oraz 
%wyliczam parametry startowe dla danego modelu
function [elast_C, b_C, u_C, s_C, x, xlabel, vectors, vectors_obs, vectors_time, restrykcje, OLS_resid, priors] = PrepareData(y, x, n, t, xlabel, model, war)
%{
  1: prosty CD
  2: CD z czasem
  3: Tom-EPI
  4: CD-LT //tu trza zadzialac!!
  5: trsnslog bez trendu
%}
dim = size(x); %liczba zmiennych objasniajacych!
if ~(dim(1) == n*t)
    msgbox('Cos nie tak...','Blad szwagier','error');
    return
end
%disp(dim);
switch model
    case 1 %simple CD
        restrykcje.restrict = NaN;
        k = dim(2)+1;
        x(:, k) = ones(dim(1),1);
        xlabel(k) = cellstr('Const');
        [elast_C, b_C, u_C, s_C, OLS_resid] = licz_COLS(y, x,0);
        vectors = zeros(k);
        for i = 1:k 
            vectors(i,i) = 1;
        end
        vectors(:,k+1) = ones(k,1); %dodaje kolumne dla skali
        vectors(k,k+1) = 0;
        vectors_obs =0;
        vectors_time=0;
        
        
        % PRIORS SETTING FOR BETA
        beta_length = size(b_C,1);
        priors.ex = 1/(beta_length-1) + zeros(beta_length,1); % expected value
        priors.ex(beta_length) = 0;
        priors.cov = zeros(beta_length,beta_length);
        priors.diag = diag(ones(beta_length,1)); %for fast matrix inverse
        for a = 1:(beta_length-1)
            for b = a:(beta_length-1)
                if a==b
                    priors.cov(a,b) = 100; %covariance matrix
                end
            end
        end
        priors.cov(beta_length,beta_length) = 100; %covariance matrix - intercept
    case 2 %CD with trend
        restrykcje.restrict = NaN;
        time = zeros(dim(1),1);
        k = dim(2)+2;
        a0 = 0;
        for a1=1:t
            for a2 = 1:n
                a0 = a0 +1;
                time(a0,1) = a1;
            end
        end
        x(:, k-1) = time;
        xlabel(k-1) = cellstr('Time');
        x(:, k) = ones(dim(1),1);
        xlabel(k) = cellstr('Const');
        [elast_C, b_C, u_C, s_C, OLS_resid] = licz_COLS(y, x,0);
        
        vectors = zeros(k);
        for i = 1:k 
            vectors(i,i) = 1;
        end
        vectors(:,k+1) = ones(k,1);
        vectors(k-1,k+1) = 0; %zero for trend
        vectors(k,k+1) = 0; %zero for the wolny
        vectors_obs =0;
        vectors_time=0;
        
        % PRIORS SETTING FOR BETA
        beta_length = size(b_C,1);
        priors.ex = 1/(beta_length-2) + zeros(beta_length,1);
        priors.ex(beta_length-1) = 0.02;
        priors.ex(beta_length) = 0;
        priors.cov = zeros(beta_length,beta_length);
        priors.diag = diag(ones(beta_length,1));
        for a = 1:(beta_length-2)
            for b = a:(beta_length-2)
                if a==b
                    priors.cov(a,b) = 100;
                end
            end
        end
        priors.cov(beta_length-1,beta_length-1) = 0.01;
        priors.cov(beta_length,beta_length) = 100;
        
    case 3 %CD Tom-EPI 
        restrykcje.restrict = NaN;
        time = zeros(dim(1),1);
        k = dim(2)+2;
        a0 = 0;
        for a1=1:t
            for a2 = 1:n
                a0 = a0 +1;
                time(a0,1) = a1;
            end
        end
        x(:, k-1) = time;
        xlabel(k-1) = cellstr('Time');
        x(:, k) = ones(dim(1),1);
        xlabel(k) = cellstr('Const');
        [elast_C, b_C, u_C, s_C, OLS_resid]= licz_COLS(y, x,0);
        
        vectors = zeros(k);
        for i = 1:k 
            vectors(i,i) = 1;
        end
        vectors(:,k+1) = ones(k,1);
        vectors(k-1,k+1) = 0; %zero tam gdzie trend
        vectors(k,k+1) = 0; %zero tam gdzie wyraz wolny
        vectors_obs =0;
        vectors_time=0;
                
    case 4 %CD-LT
        restrykcje.restrict = NaN;
        %dodaje oznaczenia
        k = 2*(dim(2)+1);
        xlabel(0.5*k) = cellstr('Const');
        xlabel(0.5*k+1) = cellstr('TimeVars');
        %dodaje jedynki
        x(:, 0.5*k) = ones(dim(1),1);
        %wyznaczam wektor "czasu"
        time = zeros(dim(1),dim(2)+1);
        a0 = 0;
        for a1=1:t
            for a2 = 1:n
                a0 = a0 +1;
                time(a0,:) = a1;
            end
        end
        xt = x.*time;
        x = [xt,x];
        %x(:, (0.5*k+1):k) = xt;
        [elast_C, b_C, u_C, s_C, OLS_resid]= licz_COLS(y, x,0);
        
        vectors = zeros(k);
        for i = 1:k 
            vectors(i,i) = 1;
        end
        vectors(:,k+1) = ones(k,1);
        vectors(0.5*k,k+1) = 0;%zero tam gdzie wyraz wolny
        vectors(k,k+1) = 0;%zero tam gdzie trend
        vectors((0.5*k+1):(k-1),k+1) = 0.5*vectors((0.5*k+1):(k-1),k+1)*t;
        vectors_obs =0;
        vectors_time=0;
        
        
        beta_length = size(b_C,1);
        % PRIORS SETTING FOR BETA
        priors.ex = 1/dim(2) + zeros(beta_length,1);
        priors.cov = diag(ones(beta_length,1)).*100;
        priors.diag = diag(ones(beta_length,1));
        
    case {5,6} %Translog without trend
        [x, vectors, vectors_obs, vectors_time, v_restrict, v_indeks, elast_C, b_C, u_C, s_C, OLS_resid] = wektory_TR(x,y,0,n,t,war);

        %restrykcje.restrict = vectors(:,:,2);
        restrykcje.restrict = v_restrict;
        restrykcje.indeksy = v_indeks;
        
        beta_length = size(b_C,1);
        priors.diag = diag(ones(beta_length,1));
        
        % PRIORS SETTING FOR BETA
        translog_priors;
        
    case {7, 8} %translog with trend
        [x, vectors, vectors_obs, vectors_time, v_restrict, v_indeks, ~, ~, ~, ~,~] = wektory_TR(x,y,0,n,t,war);

        %dokladam trend do macierzy x
        time = zeros(dim(1),1);
        a0 = 0;
        for a1=1:t
            for a2 = 1:n
                a0 = a0 +1;
                time(a0,1) = a1;
            end
        end
        x(:,size(x,2)) = [];
        x = [x,time,ones(dim(1),1)];
        
        % dokladam po "0" do wektorow tak by zdadzalo sie z wymiarem x
        place = size(vectors,2)+1;
        for a=1:size(vectors,3)
            vectors(:,place,a) = zeros(dim(1),1);
            vectors_obs(:,place,a) = zeros(n,1);
            vectors_time(:,place,a) = zeros(t,1);
        end
        for a=1:size(v_restrict,3)
            %disp(size(v_restrict(:,1,a)));
            %disp(v_restrict(:,1,a));
            v_restrict(:,place,a) = zeros(size(v_restrict(:,1,1)));
        end
        [elast_C, b_C, u_C, s_C, OLS_resid] = licz_COLS(y, x, vectors);
        restrykcje.restrict = v_restrict;
        restrykcje.indeksy = v_indeks;
        
        %Parametry a piori dla bet
        beta_length = size(b_C,1);
        priors.diag = diag(ones(beta_length,1));
        
        % PRIORS SETTING FOR BETA
        translog_priors;
        
    case {9, 10} %translogo LT
        [x, vectors, vectors_obs, vectors_time, v_restrict, v_indeks, elast_C, b_C, u_C, s_C, OLS_resid] = wektory_TR(x,y,1,n,t,war);
        restrykcje.restrict = v_restrict;
        restrykcje.indeksy = v_indeks; 
        % Potem tylko dodac opcje z restrykcjami i bez!!
        
        %Parametry a piori dla bet - DO ZMIANY
        beta_length = size(b_C,1);
        priors.diag = diag(ones(beta_length,1));
        
        % PRIORS SETTING FOR BETA
        translog_priors; 
        
    case {11, 12} %translogo QT
        % dodac parametr do wektory TR! 
        % na tej podstawie wprowadzic funkcje dodajaca zmienne do
        % popdstawowej funkcji translog bez trendu. Powinno szybko sie dac
        % rozbudowac o liczenie TR-LT i TR-QT. 
        [x, vectors, vectors_obs, vectors_time, v_restrict, v_indeks, elast_C, b_C, u_C, s_C, OLS_resid] = wektory_TR(x,y,2,n,t,war);
        %[elast_C, b_C, u_C, s_C] = licz_COLS(y, x, vectors);
        %restrykcje.restrict = vectors(:,:,2);
        restrykcje.restrict = v_restrict;
        restrykcje.indeksy = v_indeks;     
        

        beta_length = size(b_C,1);
        priors.diag = diag(ones(beta_length,1));
        
        % PRIORS SETTING FOR BETA
        translog_priors;
        
    otherwise
        msgbox('Unknown model','Well this is embarrassing...','error');
end
     
return


function [x, vectors, vectors_obs, vectors_time, v_restrict, v_indeks, elast_C, b_C, u_C, s_C, OLS_resid] = wektory_TR(x,y,opcja,n,t,war)
        dim = size(x);
        %wektory do liczenia

        pomoc = [x, ones(dim(1),1)];
        poczatek = 1;
        koniec = size(pomoc,2);
        suwak = koniec;
        wynik = [];
        for a = 1:dim(2)
            for b = 1:dim(1)
                wynik(b,poczatek:koniec) = pomoc(b,:).*x(b,a);
            end
            pomoc(:,1) = [];
            suwak = suwak - 1;
            poczatek = koniec +1;
            koniec = koniec+suwak;
        end
        vector = zeros(dim(2), size(wynik,2));
        poczatek = 0;
        skacz_po = 2;
        suwak = 0;
        help = 0;
        for a = 1:dim(2)
            for b = 1:(dim(2)+1)
                if b > skacz_po
                    suwak = suwak + help + 1;
                    help = help +1;
                end
                vector(a,b+poczatek+suwak) = 1;
                %fprintf('iteracja: %d indeks: %d suwak: %d poczatek: %d\n', b, b+poczatek+suwak,suwak, poczatek);
            end
            help=a;
            suwak =0;
            poczatek = poczatek + skacz_po;
            skacz_po = skacz_po +1;
            vector(a,poczatek) = 2;
        end

        vectorz = zeros(size(vector));
        vdim = size(vectorz);
        for a = 1:vdim(1);
            for b = 1:vdim(2);
                vectorz(a,b) = vector(vdim(1)-a+1, vdim(2)-b+1);
            end
        end
        wynik = [wynik, ones(dim(1),1)];
        vectorz = [vectorz, zeros(dim(2),1)];

        %Wektory dla zmiennych
        vectors2 = zeros(dim(1),size(wynik,2),dim(2));
        vectors = zeros(dim(1),size(wynik,2),dim(2)+1);

        for a = 1:dim(2)
            for b = 1:dim(1)
                vectors2(b,:,a) = vectorz(a,:) ./ x(b,a);
            end
        end

        for a = 1:dim(2)
            vectors(:,:,a) = vectors2(:,:,a) .* wynik ;
        end
        x = wynik;
        vectors(:,:,dim(2)+1) = sum(vectors,3); %KONIEC
%--to sie odpala w przypadku gdy mam do czynienia z TR-LT lub TR-QT
        switch opcja
            case 1 %LT
                time = zeros(dim(1),1);
                a0 = 0;
                for a1=1:t
                    for a2 = 1:n
                        a0 = a0 +1;
                        time(a0,1) = a1;
                    end
                end
                x_t = zeros(size(x));
                vectors_t = zeros(size(vectors));
                for a=1:(n*t)
                    x_t(a,:) = x(a,:)*time(a);
                    vectors_t(a,:) = vectors(a,:).*time(a);
                end
                x = [x_t, x];
                vectors = [vectors_t, vectors];
        case 2 %QT
                time = zeros(dim(1),1);
                a0 = 0;
                for a1=1:t
                    for a2 = 1:n
                        a0 = a0 +1;
                        time(a0,1) = a1;
                    end
                end
                x_t = zeros(size(x));
                x_t2 = zeros(size(x));
                vectors_t = zeros(size(vectors));
                vectors_t2 = zeros(size(vectors));
                for a=1:(n*t)
                    x_t(a,:) = x(a,:)*time(a);
                    x_t2(a,:) = x(a,:)*(time(a)^2);
                    vectors_t(a,:) = vectors(a,:)*time(a);
                    vectors_t2(a,:) = vectors(a,:)*(time(a)^2);
                end
                x = [x_t2, x_t, x];
                vectors = [vectors_t2, vectors_t, vectors];
            otherwise
        end
%--koniec        
        
%wyznaczam wektory dla elastycznosci i skali po czasie i po obserwacjach
v_dim = size(vectors);
vectors_obs = zeros(n,v_dim(2),v_dim(3));
vectors_time = zeros(t,v_dim(2),v_dim(3));

c1 = 1;
for c = 1:n:(n*t)
    %disp(size(vectors_obs));
    %disp(size(vectors_time));
    %disp(size(vectors(c:(c+n-1),:,:)));
    vectors_obs = vectors_obs + vectors(c:(c+n-1),:,:)./t;
    vectors_time(c1,:,:) = vectors_time(c1,:,:) + mean(vectors(c:(c+n-1),:,:),1);
    c1 = c1+1;
end
%koniec
        
        [elast_C, b_C, u_C, s_C, OLS_resid] = licz_COLS(y, x, vectors);
        
        %----- to powinno zosta? wydzielone i dopuszcza? kilka mo?liwo?ci
        %----- wprowadzenia restrykcji, np na ?rednie elastyczno?ci :)
        
        %NOWE
%case 1: po dwie restrykcje na zmienna 
%case 2: po dwie restrykcje na zmienna, dla 'm' pierwszych obiektow
%case 3: restrykcja w calosci na jedna ze zmiennych
%case 4: restrykcja w calosci na jedna ze zmiennych
%        dla 'm' pierwszych obiektow
%case 5: restr. dla srednich po obserwacjach, po 1 na zmienna
%case 6: restr. dla srednich po obserwacjach, po 1 na zmienna
%        dla 'm' pierwszych obiektow
%case 7: restr. dla srednich po obserwacjach, wszystkie
%case 8: restr. dla srednich po obserwacjach, wszystkie
%        dla 'm' pierwszych obiektow
        %disp(war);
        %pause;
        %war.case = 4; %sposob narzuczenia restrykcji
        %war.m = 27; %liczba objektow podlegajaca restrykcjom (dla case=3)
        %war.zm = 2; %zmienna na ktorej maja byc restrykcje (case 2 i 3)
        [restrict, indeks] = restrykcje(n, t, vectors, vectors_obs, b_C, war, dim);
        
        v_restrict = restrict;
        %v_restrict = zeros(2*a(3),a(2));
        %v_restrict(1:a(3),1:a(2))= one(1:a(3),1:a(2)); %przepisuje pierwszy zestaw restrykcji
        v_indeks = indeks;
        
        %to jest recznie wpisany warunek gwarantujacy nieujemnosc
        %elastycznosci pracy, warunek sprawdza caly wektor elastycznosci
        %pracy, do zmiany i do wprowadzenia do wydzielonego miejsca na
        %restrykcje
        %if opcja == 1 || opcja == 2
        %    v_restrict = vectors(:,:,2);
        %end
return