function ret = decompose2(x,n,t,folder,pliki,ile_danych,model,factors,vecs,xRaw)
% Funkcja powinna byc wywolana w nastepujacy sposob:
% ret = decompose2(user.x,user.n,user.t,user.folder,user.pliki,user.ile_danych,user.model,user.factors,user.vecs,user.xRaw);
% gdzie struktura user pochodzi z pliku wyniki.mat. Wystarczy wiec wczytac
% w/w zmienna strukturalna i wywolac funkcje decompose wedlu podanego
% schematu. 

%-- Przygotowuje zmienne
ret = struct;
l_plikow = size(ile_danych,1);
folder = strcat(folder,'\decomp');
if ~isdir(folder)
    mkdir(folder);
    disp('Folder for decompostion results has been created');
end
for a2 = 1:l_plikow
    str = strrep(sprintf('%2.3f', a2/1000),'.','0');
    
    pliki(a2).logEC = strcat(folder,'\','logEC',str);
    pliki(a2).logAEC = strcat(folder,'\','logAEC',str);
    pliki(a2).logIC = strcat(folder,'\','logIC',str);
    pliki(a2).logAIC = strcat(folder,'\','logAIC',str);
    pliki(a2).logOC = strcat(folder,'\','logOC',str);
    pliki(a2).logAOC = strcat(folder,'\','logAOC',str);
    pliki(a2).logTC = strcat(folder,'\','logTC',str);
    pliki(a2).logATC = strcat(folder,'\','logATC',str);
    pliki(a2).logPC = strcat(folder,'\','logPC',str);
    pliki(a2).logAPC = strcat(folder,'\','logAPC',str);
    pliki(a2).logSSC = strcat(folder,'\','logSSC',str);
    pliki(a2).log_aSSC = strcat(folder,'\','log_aSSC',str);
    pliki(a2).logSEC = strcat(folder,'\','logSEC',str);
    pliki(a2).log_aSEC = strcat(folder,'\','log_aSEC',str);
    
end
ret.EC = 0; ret.AEC = 0;
ret.OC = 0; ret.AOC = 0;
ret.IC = 0; ret.AIC = 0;
ret.TC = 0; ret.ATC = 0;
ret.PC = 0; ret.APC = 0;
% Nowe
ret.SSC = 0; ret.aSSC = 0;
ret.SEC = 0; ret.aSEC = 0;
% Nowe
ret.PIC = 0; ret.aPIC = 0;
ret.RTC = 0; ret.aRTC = 0;
ret.TSC = 0; ret.aTSC = 0;
ret.TPC = 0; ret.aTPC = 0;

%----DEKOMPOZYCJA DLA EFFICIENCY CHANGE I OUTPUT CHANGE
[ret.EC, ret.AEC, ret.OC, ret.AOC] = licz_EiO(x,n,t,pliki,l_plikow,ile_danych);

%----DEKOMPOZYCJA DLA INPUT CHANGE I TECHNICAL CHANGE
switch model
    case {1,2,3,4}
        [ret.IC, ret.AIC, ret.TC, ret.ATC] = input_CD(x,n,t,pliki,ile_danych,factors, l_plikow, model);
        %[ret.PIC, ret.APIC, ret.PTC, ret.APTC, ret.SEC, ret.ASEC, ret.SSC, ret.ASSC] = input2_CD(x,n,t,pliki,ile_danych,factors, l_plikow, model);
    %case 4
    %    [ret.IC, ret.AIC] = input_CDTL(x,n,t,pliki,ile_danych,factors, l_plikow);
    case {5,6,7,8,9,10,11,12}
        [ret.IC, ret.AIC, ret.TC, ret.ATC, ret.SSC, ret.aSSC, ret.SEC, ret.aSEC, ret.RTC, ret.aRTC, ret.TSC, ret.aTSC, ret.PIC, ret.aPIC] = input_TR(xRaw, x,n,t,pliki,ile_danych,factors, l_plikow, vecs, model);
        %msgbox('Translog');
    otherwise
        msgbox('Unknown model');
end

[ret.PC, ret.APC, ret.TPC, ret.aTPC] = licz_PC(n,t,pliki,l_plikow,ile_danych,model);


end


function [PC, APC, TPC, aTPC] = licz_PC(n,t,pliki,l_plikow,ile_danych,model)
PC = zeros(n*(t-1),2);
APC = zeros(n,2);
TPC = zeros(n*(t-1),2);
aTPC = zeros(n,2);

h = waitbar(0,'Sted 3/3: PRODUCTIVITY',...
            'Name','DECOMPOSITION',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
assignin('base', 'Waitbar_handle2', h);

switch model
    case {2,3,6,8}
       for b = 1:l_plikow
            %disp('HALO petla');
            if ile_danych(b,1) ~= 0
                var = load(pliki(b).logEC,'LogEC');
                logEC = var.LogEC;
                var = load(pliki(b).bety,'b_B');
                rozmiar = size(var.b_B);
                zacznij = rozmiar(2) - ile_danych(b,1)+1;
                time_const = var.b_B(rozmiar(1)-1,zacznij:rozmiar(2));
                var = load(pliki(b).logSEC,'LogSEC');
                LogSEC = var.LogSEC;
                clear var;
                %disp(size(time_const));
                %disp(size(logEC));
                LogPC = logEC+ones(n*(t-1),1)*time_const;
                LogTPC = LogPC+LogSEC;
                clear logEC LogSEC;
                save(pliki(b).logPC, 'LogPC');
                %tu moznna uzyc PC do wyznaczenia APC
                LogPC = 100*(exp(LogPC)-1);
                PC(:,1) = PC(:,1)+ mean(LogPC,2).*ile_danych(b,2);
                PC(:,2) = PC(:,2)+ std(LogPC,1,2).*ile_danych(b,2);
                LogTPC = 100*(exp(LogTPC)-1);
                TPC(:,1) = TPC(:,1)+ mean(LogTPC,2).*ile_danych(b,2);
                TPC(:,2) = TPC(:,2)+ std(LogTPC,1,2).*ile_danych(b,2);
                clear LogPC LogTPC;
                var = load(pliki(b).logAEC,'LogAEC');
                LogAEC = var.LogAEC;
                var = load(pliki(b).log_aSEC,'Log_aSEC');
                Log_aSEC = var.Log_aSEC;
                clear var;
                LogAPC = LogAEC+ones(n,1)*time_const;
                clear logAEC logATC;
                save(pliki(b).logAPC, 'LogAPC');
                Log_aTPC = LogAPC+Log_aSEC;
                clear Log_aSEC;
                %tu moznna uzyc PC do wyznaczenia APC
                LogAPC = 100*(exp(LogAPC)-1);
                APC(:,1) = APC(:,1)+ mean(LogAPC,2).*ile_danych(b,2);
                APC(:,2) = APC(:,2)+ std(LogAPC1,2).*ile_danych(b,2);
                Log_aTPC = 100*(exp(Log_aTPC)-1);
                aTPC(:,1) = aTPC(:,1)+ mean(Log_aTPC,2).*ile_danych(b,2);
                aTPC(:,2) = aTPC(:,2)+ std(Log_aTPC,1,2).*ile_danych(b,2);
                
                clear LogAPC;
            end
            counter = b/l_plikow;
            msg = sprintf('Step 3/3: %d%%',round(counter*100));
            waitbar(counter,h,msg);
        end
    case {4,9,10,11,12}
        for b = 1:l_plikow
            %disp('HALO petla');
            if ile_danych(b,1) ~= 0
                var = load(pliki(b).logEC,'LogEC');
                disp(var);
                logEC = var.LogEC;
                var = load(pliki(b).logTC,'LogTC');
                logTC = var.LogTC;
                var = load(pliki(b).logSEC,'LogSEC');
                LogSEC = var.LogSEC;
                clear var;
                LogPC = logEC+logTC;
                LogTPC = LogPC+LogSEC;
                clear logEC logTC LogSEC;
                save(pliki(b).logPC, 'LogPC');
                %tu moznna uzyc PC do wyznaczenia APC
                LogPC = 100*(exp(LogPC)-1);
                PC(:,1) = PC(:,1)+ mean(LogPC,2).*ile_danych(b,2);
                PC(:,2) = PC(:,2)+ std(LogPC,1,2).*ile_danych(b,2);
                LogTPC = 100*(exp(LogTPC)-1);
                TPC(:,1) = TPC(:,1)+ mean(LogTPC,2).*ile_danych(b,2);
                TPC(:,2) = TPC(:,2)+ std(LogTPC,1,2).*ile_danych(b,2);
                clear LogPC LogTPC;

                var = load(pliki(b).logAEC,'LogAEC');
                LogAEC = var.LogAEC;
                var = load(pliki(b).logATC,'LogATC');
                logATC = var.LogATC;
                var = load(pliki(b).log_aSEC,'Log_aSEC');
                log_aSEC = var.Log_aSEC;
                clear var;
                LogAPC = LogAEC+logATC;
                Log_aTPC = LogAPC + log_aSEC;
                clear logAEC logATC log_aSEC;
                save(pliki(b).logAPC, 'LogAPC');
                %tu moznna uzyc PC do wyznaczenia APC
                LogAPC = 100*(exp(LogAPC)-1);
                APC(:,1) = APC(:,1)+ mean(LogAPC,2).*ile_danych(b,2);
                APC(:,2) = APC(:,2)+ std(LogAPC,1,2).*ile_danych(b,2);
                
                Log_aTPC = 100*(exp(Log_aTPC)-1);
                aTPC(:,1) = aTPC(:,1)+ mean(Log_aTPC,2).*ile_danych(b,2);
                aTPC(:,2) = aTPC(:,2)+ std(Log_aTPC,1,2).*ile_danych(b,2);
                clear LogAPC Log_aTPC;
            end
            counter = b/l_plikow;
            msg = sprintf('Step 3/3: %d%%',round(counter*100));
            waitbar(counter,h,msg);
        end
    otherwise
        PC = 0;
        APC = 0;
end
delete(h);
end

function [EC, AEC, OC, AOC] = licz_EiO(x,n,t,pliki,l_plikow,ile_danych)
EC = zeros(n*(t-1),2);
AEC = zeros(n,2);
OC = zeros(n*(t-1),2);
AOC = zeros(n,2);

h = waitbar(0,'Step 1/3: output and efficiency',...
            'Name','DECOMPOSTION',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
assignin('base', 'Waitbar_handle2', h);


for b = 1:l_plikow
    %disp('HALO petla');
    if ile_danych(b,1) ~= 0
        var = load(pliki(b).u,'u_B');
        %disp(size(var.b_B));
        %disp(size(vectors));
        %wyznaczam koniec i poczatek sczytywania danych
        rozmiar = size(var.u_B,2);
        zacznij = rozmiar - ile_danych(b,1)+1;
        ineff = var.u_B(:,zacznij:rozmiar);
        
        %wyliczam OC
        var = load(pliki(b).bety,'b_B');
        bety = var.b_B(:,zacznij:rozmiar);
        clear var;
        gdp_logest = x*bety-ineff;
        clear bety;
        LogOC = gdp_logest((n+1):(n*t),:) - gdp_logest(1:((t-1)*n),:);
        save(pliki(b).logOC, 'LogOC');
        LogOC = 100*(exp(LogOC)-1);
        OC(:,1) = OC(:,1)+ mean(LogOC,2).*ile_danych(b,2);
        OC(:,2) = OC(:,2)+ std(LogOC')'.*ile_danych(b,2);
        clear LogOC;
        LogAOC = (gdp_logest((n*(t-1)+1):n*t,:) - gdp_logest(1:n,:))./(t-1);
        save(pliki(b).logAOC, 'LogAOC');
        LogAOC = 100*(exp(LogAOC)-1);
        AOC(:,1) = AOC(:,1)+ mean(LogAOC,2).*ile_danych(b,2);
        AOC(:,2) = AOC(:,2)+ std(LogAOC')'.*ile_danych(b,2);
        clear LogAOC;
        
        %wyznaczam OC
        %disp(size(eff));
        %logEC  = u_ti - u_(t+1)i; bo "-u_ti";
        LogEC =  ineff(1:((t-1)*n),:) - ineff((n+1):(n*t),:); %tu logEC przestaje byc juz logarytmem
        save(pliki(b).logEC, 'LogEC');
        LogEC = 100*(exp(LogEC)-1);
        EC(:,1) = EC(:,1)+ mean(LogEC,2).*ile_danych(b,2);
        EC(:,2) = EC(:,2)+ std(LogEC')'.*ile_danych(b,2);
        clear LogEC;
        LogAEC = (ineff(1:n,:) - ineff((n*(t-1)+1):n*t,:))/(t-1);
        save(pliki(b).logAEC, 'LogAEC');
        LogAEC = 100*(exp(LogAEC)-1);
        AEC(:,1) = AEC(:,1)+ mean(LogAEC,2).*ile_danych(b,2);
        AEC(:,2) = AEC(:,2)+ std(LogAEC')'.*ile_danych(b,2);
        clear LogAEC;
    end
    counter = b/l_plikow;
    msg = sprintf('Step 1/3: %d%%',round(counter*100));
    waitbar(counter,h,msg);
end
delete(h);
end


function [IC, AIC, TC, ATC, SSC, aSSC, SEC, aSEC, RTC, aRTC, TSC, aTSC, PIC, aPIC] = input_TR(x,xTR,n,t,pliki,ile_danych,factors, l_plikow, vecs, model)
TC = 0; ATC = 0;
x_delta = x((n+1):(n*t),:) - x(1:((t-1)*n),:);
x_delta_av = (x((n*(t-1)+1):n*t,:,:) - x(1:n,:,:))./(t-1);
dim(1) = size(x_delta,1);
dim(2) = size(vecs.all,2);
dim(3) = factors;
vecs_delta = zeros(dim(1),dim(2),dim(3)); %tu trzeba dodac dwa wymiary ostani na poziomie
vecs_delta_av = zeros(n,dim(2),dim(3));   %tu trzeba dodac dwa wymiary 
for i = 1:dim(3)
    vecs_delta(:,:,i) = (vecs.all((n+1):(n*t),:,i) + vecs.all(1:((t-1)*n),:,i)).*0.5;
    vecs_delta_av(:,:,i) = (vecs.all((n*(t-1)+1):n*t,:,i) + vecs.all(1:n,:,i)).*0.5;
    for a=1:dim(1)
        vecs_delta(a,:,i) = vecs_delta(a,:,i).*x_delta(a,i);
    end
    for a=1:n
        vecs_delta_av(a,:,i) = vecs_delta_av(a,:,i).*x_delta_av(a,i);
    end
end
%dla nowej dekompozycji
vecs_delta(:,:,dim(3)+1) = (vecs.all((n+1):(n*t),:,dim(3)+1) + vecs.all(1:((t-1)*n),:,dim(3)+1)).*0.5; %wektore dla sr_x_RTS
vecs_delta(:,:,dim(3)+2) = (vecs.all((n+1):(n*t),:,dim(3)+1) - vecs.all(1:((t-1)*n),:,dim(3)+1)); %wektor dla dif_x_RTS
vecs_delta_av(:,:,dim(3)+1) = (vecs.all((n*(t-1)+1):n*t,:,dim(3)+1) + vecs.all(1:n,:,dim(3)+1)).*0.5; %wektor dla sr_x_RTS dla srednich
vecs_delta_av(:,:,dim(3)+2) = (vecs.all((n*(t-1)+1):n*t,:,dim(3)+1) - vecs.all(1:n,:,dim(3)+1))./(t-1); %wektor dla diff_x_RTS dla srednich
disp(size(vecs_delta(:,:,dim(3)+2)));
switch model
    case {5,6,7,8} %zwykle translogi
        [IC, AIC] = IC_static(vecs_delta,vecs_delta_av, n, t, pliki, ile_danych, l_plikow, dim);
    case {9,10} %TR_LT
        %docinam macierz wedlug konstrukcji modelu wyrzucajac kolumny t*x
        vecs_delta(:,1:(0.5*dim(2)),:) = [];
        vecs_delta_av(:,1:(0.5*dim(2)),:) = [];
        %disp(vecs_delta_av);
        x_mid = (xTR((n+1):(n*t),:) + xTR(1:((t-1)*n),:)).*0.5;
        x_mid(:,1:(0.5*dim(2))) = [];
        x_mid_av = (xTR((n*(t-1)+1):n*t,:,:) + xTR(1:n,:,:)).*0.5;
        x_mid_av(:,1:(0.5*dim(2))) = [];
        %disp(size(x_mid));
        %disp(x_delta_av_comp);
        [IC, AIC, TC, ATC, SSC, aSSC, SEC, aSEC, RTC, aRTC, TSC, aTSC, PIC, aPIC] = IC_LT(vecs_delta,vecs_delta_av, x_mid, x_mid_av, n, t, pliki, ile_danych, l_plikow, dim);
    case {11,12} %TR_QT
        to_cut = 2*dim(2)/3;
        vecs_delta(:,1:to_cut,:) = [];
        vecs_delta_av(:,1:to_cut,:) = [];
        %disp(vecs_delta_av);
        x_mid = (xTR((n+1):(n*t),:) + xTR(1:((t-1)*n),:)).*0.5;
        x_mid(:,1:to_cut) = [];
        x_mid_av = (xTR((n*(t-1)+1):n*t,:,:) + xTR(1:n,:,:)).*0.5;
        x_mid_av(:,1:to_cut) = [];
        %disp(size(x_mid));
        %disp(x_delta_av_comp);
        [IC, AIC, TC, ATC] = IC_QT(vecs_delta,vecs_delta_av, x_mid, x_mid_av, n, t, pliki, ile_danych, l_plikow, dim);
    otherwise
        msgbox('Unknown model');
        %disp(vecs_delta);
        
end
%assignin('base','vecs_delta',vecs_delta);
%assignin('base','vecs_dl_av',vecs_delta_av);
%[IC, AIC] = IC_static(vecs_delta, vecs_delta_av, n, t, pliki, ile_danych, l_plikow, dim);

end


function [IC, AIC, TC, ATC] = input_CD(x,n,t,pliki,ile_danych,factors, l_plikow, model)
TC = 0; ATC = 0;
%---DEKOMPOZYCJA DLA INPUT CHANGE
%--Przygotowuje macierze
x_delta = x((n+1):(n*t),:) - x(1:((t-1)*n),:);
x_delta_av = (x((n*(t-1)+1):n*t,:) - x(1:n,:))./(t-1);
dim = size(x_delta);
dim(3) = factors;
x_delta_comp = zeros(dim(1),dim(2),dim(3));
x_delta_av_comp = zeros(n,dim(2),dim(3));

if model == 4 %CD-LT
    for i = 1:dim(3)
        x_delta_comp(:,i,i) = x_delta(:,i+0.5*dim(2));%bo w LT x znajduja sie od polowy macierzy
        x_delta_av_comp(:,i,i) = x_delta_av(:,i+0.5*dim(2));
    end
    %docinam macierz wedlug konstrukcji modelu wyrzucajac kolumny t*x
    x_delta_comp(:,(0.5*dim(2)+1):dim(2),:) = [];
    x_delta_av_comp(:,(0.5*dim(2)+1):dim(2),:) = [];
    
    x_mid = (x((n+1):(n*t),:) + x(1:((t-1)*n),:)).*0.5;
    x_mid(:,1:(0.5*dim(2))) = [];
    x_mid_av = (x((n*(t-1)+1):n*t,:,:) + x(1:n,:,:)).*0.5; 
    x_mid_av(:,1:(0.5*dim(2))) = [];
    %disp(x_delta_av_comp);
    
    [IC, AIC, TC, ATC] = IC_LT(x_delta_comp,x_delta_av_comp,x_mid,x_mid_av, n, t, pliki, ile_danych, l_plikow, dim);
else
    for i = 1:dim(3)
        x_delta_comp(:,i,i) = x_delta(:,i);
        x_delta_av_comp(:,i,i) = x_delta_av(:,i);
    end
    [IC, AIC] = IC_static(x_delta_comp,x_delta_av_comp, n, t, pliki, ile_danych, l_plikow, dim);
end

end


function [IC, AIC, TC, ATC, SSC, aSSC, SEC, aSEC, RTC, aRTC, TSC, aTSC, PIC, aPIC] = IC_LT(x_delta_comp,x_delta_av_comp, x_mid, x_mid_av, n, t, pliki, ile_danych, l_plikow, dim)
%to jest kod wylaczne dla LT! Powodem jest to ze bety po pobraniu dzielone
%sa na dwie czesci (bety stale i te z trendem) wedlug konstrukcji macierzy
%x w funkcji PrepairData
IC = zeros(n*(t-1),2,dim(3)+1);
AIC = zeros(n,2,dim(3)+1);
TC = zeros(n*(t-1),2);
ATC = zeros(n,2);
SSC = zeros(n*(t-1),2);
aSSC = zeros(n,2);
SEC = zeros(n*(t-1),2);
aSEC = zeros(n,2);
RTC = zeros(n*(t-1),2);
aRTC = zeros(n,2);
TSC = zeros(n*(t-1),2);
aTSC = zeros(n,2);
PIC = zeros(n*(t-1),2);
aPIC = zeros(n,2);

h = waitbar(0,'Step 2/3: INPUT, TECH',...
            'Name','DECOMPOSITION',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
assignin('base', 'Waitbar_handle2', h);

for b = 1:l_plikow
    %disp('HALO petla');
    if ile_danych(b,1) ~= 0
        %rozm(2) = ile_danych(b,1);
        var = load(pliki(b).bety,'b_B');
        %disp(size(var.b_B));
        %disp(size(vectors));
        %wyznaczam koniec i poczatek sczytywania danych
        rozmiar = size(var.b_B,2);
        zacznij = rozmiar - ile_danych(b,1)+1;
        bety = var.b_B(:,zacznij:rozmiar);
        %rozm = size(bety);
        beta_const = bety((0.5*dim(2)+1):dim(2),:); %bo bety dla x sa na koncu a dla xt na poczatku 
        beta_var = bety(1:(0.5*dim(2)),:);
        %disp(size(beta_const));
        clear bety;
        
        LogTC = x_mid*beta_var;
        save(pliki(b).logTC, 'LogTC');
        %LogTC = 100*(exp(LogTC)-1);
        TC(:,1) = TC(:,1)+ mean(100*(exp(LogTC)-1),2).*ile_danych(b,2);
        TC(:,2) = TC(:,2)+ std(100*(exp(LogTC)-1),1,2).*ile_danych(b,2);
        
        %DLA NOWEJ DEKOMPOZYCJI
        LogSSC = x_delta_comp(:,:,dim(3)+1)*beta_var;
        save(pliki(b).logSSC, 'LogSSC');
        %LogSSC = 100*(exp(LogSSC)-1);
        SSC(:,1) = SSC(:,1)+mean(100*(exp(LogSSC)-1),2).*ile_danych(b,2);
        SSC(:,2) = SSC(:,2)+std(100*(exp(LogSSC)-1),1,2).*ile_danych(b,2);
        
        LogRTC = LogTC-LogSSC;
        RTC(:,1) = RTC(:,1)+mean(100*(exp(LogRTC)-1),2).*ile_danych(b,2);
        RTC(:,2) = RTC(:,2)+std(100*(exp(LogRTC)-1),1,2).*ile_danych(b,2);
        
        clear LogTC LogRTC;
        
        c=1;
        LogSEC = zeros(dim(1),ile_danych(b,1));
        for a = 1:(t-1)
            LogSEC(c:(c+n-1),:) = x_delta_comp(c:(c+n-1),:,dim(3)+2)*(beta_const+(a+0.5).*beta_var);
            c = c+n;
        end
        save(pliki(b).logSEC, 'LogSEC');
        %LogSEC = 100*(exp(LogSEC)-1);
        SEC(:,1) = SEC(:,1)+mean(100*(exp(LogSEC)-1),2).*ile_danych(b,2);
        SEC(:,2) = SEC(:,2)+std(100*(exp(LogSEC)-1),1,2).*ile_danych(b,2);
        
        LogTSC = LogSEC+LogSSC;
        TSC(:,1) = TSC(:,1)+mean(100*(exp(LogTSC)-1),2).*ile_danych(b,2);
        TSC(:,2) = TSC(:,2)+std(100*(exp(LogTSC)-1),1,2).*ile_danych(b,2);
        clear LogSSC LogTSC;
        
        LogIC = zeros(dim(1),ile_danych(b,1),dim(3)+1); %"+1" poniewaz ostatni wymiar zawiera sume, czyli laczny wplyw dla IC
        %disp(size(eff));
        %logEC  = u_ti - u_(t+1)i; bo "-u_ti";
        for i = 1:dim(3)
            c = 1;
            for a = 1:(t-1)
                LogIC(c:(c+n-1),:,i) = x_delta_comp(c:(c+n-1),:,i)*(beta_const+(a+0.5).*beta_var);
                LogIC(c:(c+n-1),:,dim(3)+1) = LogIC(c:(c+n-1),:,dim(3)+1)+LogIC(c:(c+n-1),:,i);
                c = c+n;
            end
        end
        save(pliki(b).logIC, 'LogIC', '-v7.3');
        LogICc = 100*(exp(LogIC)-1);
        %disp(size(LogIC));
        %disp(size(ret.IC));
        for i =1:(dim(3)+1)
            IC(:,1,i) = IC(:,1,i)+ mean(LogICc(:,:,i),2).*ile_danych(b,2);
            IC(:,2,i) = IC(:,2,i)+ std(LogICc(:,:,i),1,2).*ile_danych(b,2);
        end
        clear LogICc;
        LogPIC = LogIC(:,:,dim(3)+1) - LogSEC;
        PIC(:,1) = PIC(:,1)+mean(100*(exp(LogPIC)-1),2).*ile_danych(b,2);
        PIC(:,2) = PIC(:,2)+std(100*(exp(LogPIC)-1),1,2).*ile_danych(b,2);
        
        clear LogIC LogSEC LogPIC;
        
        %AVERAGES
        
        LogATC = x_mid_av*beta_var;
        save(pliki(b).logATC, 'LogATC');
        %LogATC = 100*(exp(LogATC)-1);
        ATC(:,1) = ATC(:,1)+ mean(100*(exp(LogATC)-1),2).*ile_danych(b,2);
        ATC(:,2) = ATC(:,2)+ std(100*(exp(LogATC)-1),1,2).*ile_danych(b,2);
        %clear LogATC;
        
        %Dla nowej dekompozycji
        Log_aSSC = x_delta_av_comp(:,:,dim(3)+1)*beta_var;
        save(pliki(b).log_aSSC, 'Log_aSSC');
        %Log_aSSC = 100*(exp(Log_aSSC)-1);
        aSSC(:,1) = aSSC(:,1)+mean(100*(exp(Log_aSSC)-1),2).*ile_danych(b,2);
        aSSC(:,2) = aSSC(:,2)+std(100*(exp(Log_aSSC)-1),1,2).*ile_danych(b,2);
        
        Log_aRTC = LogATC-Log_aSSC;
        aRTC(:,1) = aRTC(:,1)+mean(100*(exp(Log_aRTC)-1),2).*ile_danych(b,2);
        aRTC(:,2) = aRTC(:,2)+std(100*(exp(Log_aRTC)-1),1,2).*ile_danych(b,2);
        clear LogATC Log_aRTC;
        
        Log_aSEC = x_delta_av_comp(:,:,dim(3)+2)*(beta_const+(0.5*t+0.5).*beta_var);
        save(pliki(b).log_aSEC, 'Log_aSEC');
        %Log_aSEC = 100*(exp(Log_aSEC)-1);
        aSEC(:,1) = aSEC(:,1)+mean(100*(exp(Log_aSEC)-1),2).*ile_danych(b,2);
        aSEC(:,2) = aSEC(:,2)+std(100*(exp(Log_aSEC)-1),1,2).*ile_danych(b,2);
        
        Log_aTSC = Log_aSEC+Log_aSSC;
        aTSC(:,1) = aTSC(:,1)+mean(100*(exp(Log_aTSC)-1),2).*ile_danych(b,2);
        aTSC(:,2) = aTSC(:,2)+std(100*(exp(Log_aTSC)-1),1,2).*ile_danych(b,2);
        clear Log_aSSC Log_aTSC;
        
        LogAIC = zeros(n,ile_danych(b,1),dim(3)+1);
        for i = 1:dim(3)
                LogAIC(:,:,i) = x_delta_av_comp(:,:,i)*(beta_const+(0.5*t+0.5).*beta_var);
                LogAIC(:,:,dim(3)+1) = LogAIC(:,:,dim(3)+1)+LogAIC(:,:,i);
        end
        save(pliki(b).logAIC, 'LogAIC');
        LogAICc = 100*(exp(LogAIC)-1);
        for i =1:(dim(3)+1)
            AIC(:,1,i) = AIC(:,1,i)+ mean(LogAICc(:,:,i),2).*ile_danych(b,2);
            AIC(:,2,i) = AIC(:,2,i)+ std(LogAICc(:,:,i)')'.*ile_danych(b,2);
        end
        clear LogAICc;
        
        Log_aPIC = LogAIC(:,:,dim(3)+1) - Log_aSEC;
        aPIC(:,1) = aPIC(:,1)+mean(100*(exp(Log_aPIC)-1),2).*ile_danych(b,2);
        aPIC(:,2) = aPIC(:,2)+std(100*(exp(Log_aPIC)-1),1,2).*ile_danych(b,2);
        %Dla nowej dekompozycji
        clear Log_aPIC LogAIC Log_aSEC; 
        
        
    end
    counter = b/l_plikow;
    msg = sprintf('Step 2/3: %d%%',round(counter*100));
    waitbar(counter,h,msg);
end
delete(h);
end

function [IC, AIC] = IC_static(x_delta_comp,x_delta_av_comp, n, t, pliki, ile_danych, l_plikow, dim)

h = waitbar(0,'Step 2/3: INPUT, (TECH=const)',...
            'Name','DECOMPOSITION',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
assignin('base', 'Waitbar_handle2', h);

IC = zeros(n*(t-1),2,dim(3)+1);
AIC = zeros(n,2,dim(3)+1);
for b = 1:l_plikow
    %disp('HALO petla');
    if ile_danych(b,1) ~= 0
        var = load(pliki(b).bety,'b_B');
        %disp(size(var.b_B));
        %disp(size(vectors));
        %wyznaczam koniec i poczatek sczytywania danych
        rozmiar = size(var.b_B,2);
        zacznij = rozmiar - ile_danych(b,1)+1;
        bety = var.b_B(:,zacznij:rozmiar);
        rozm = size(bety,2);
        LogIC = zeros(dim(1),ile_danych(b,1),dim(3)+1);%"+1" poniewaz ostatni wymiar zawiera sume, czyli laczny wplyw dla IC
        %disp(size(eff));
        %logEC  = u_ti - u_(t+1)i; bo "-u_ti";
        for i = 1:dim(3)
            LogIC(:,:,i) = x_delta_comp(:,:,i)*bety;
            LogIC(:,:,dim(3)+1) = LogIC(:,:,dim(3)+1)+LogIC(:,:,i);
        end
        save(pliki(b).logIC, 'LogIC');
        LogIC = 100*(exp(LogIC)-1);
        %disp(size(LogIC));
        %disp(size(ret.IC));
        for i =1:(dim(3)+1)
            IC(:,1,i) = IC(:,1,i)+ mean(LogIC(:,:,i),2).*ile_danych(b,2);
            IC(:,2,i) = IC(:,2,i)+ std(LogIC(:,:,i)')'.*ile_danych(b,2);
        end
        clear LogIC;
        LogAIC = zeros(n,ile_danych(b,1),dim(3)+1);
        for i = 1:dim(3)
            LogAIC(:,:,i) = x_delta_av_comp(:,:,i)*bety;
            LogAIC(:,:,dim(3)+1) = LogAIC(:,:,dim(3)+1)+LogAIC(:,:,i);
        end
        save(pliki(b).logAIC, 'LogAIC');
        LogAIC = 100*(exp(LogAIC)-1);
        for i =1:(dim(3)+1)
            AIC(:,1,i) = AIC(:,1,i)+ mean(LogAIC(:,:,i),2).*ile_danych(b,2);
            AIC(:,2,i) = AIC(:,2,i)+ std(LogAIC(:,:,i)')'.*ile_danych(b,2);
        end
    end
    counter = b/l_plikow;
    msg = sprintf('Step 2/3: %d%%',round(counter*100));
    waitbar(counter,h,msg);
end
delete(h);

end


function [IC, AIC, TC, ATC] = IC_QT(x_delta_comp,x_delta_av_comp, x_mid, x_mid_av, n, t, pliki, ile_danych, l_plikow, dim)
%to jest kod wylaczne dla LT! Powodem jest to ze bety po pobraniu dzielone
%sa na dwie czesci (bety stale i te z trendem) wedlug konstrukcji macierzy
%x w funkcji PrepairData
IC = zeros(n*(t-1),2,dim(3)+1);
AIC = zeros(n,2,dim(3)+1);
TC = zeros(n*(t-1),2);
ATC = zeros(n,2);

%----------
cut1 = dim(2)/3;
cut2 = 2*dim(2)/3;
%disp(dim(2)); disp(cut1); disp(cut2);
%----------

h = waitbar(0,'Step 2/3: INPUT, TECH',...
            'Name','DECOMPOSITION',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
assignin('base', 'Waitbar_handle2', h);

for b = 1:l_plikow
    %disp('HALO petla');
    if ile_danych(b,1) ~= 0
        var = load(pliki(b).bety,'b_B');
        %disp(size(var.b_B));
        %disp(size(vectors));
        %wyznaczam koniec i poczatek sczytywania danych
        rozmiar = size(var.b_B,2);
        zacznij = rozmiar - ile_danych(b,1)+1;
        bety = var.b_B(:,zacznij:rozmiar);
        %srozm = size(bety); %moze byc zbedne
        
        %----------
        beta_t2 = bety(1:cut1,:);
        beta_t = bety((cut1+1):cut2,:);
        beta_const = bety((cut2+1):dim(2),:); %bo bety dla x sa na koncu a dla xt na poczatku 
        %---------
        
        %disp(size(beta_const));
        clear bety;
        
        %to sie sporo zmieni...
        %----------
        c=1;
        for a = 1:(t-1)
            LogTC(c:(c+n-1),:) = x_mid(c:(c+n-1),:)*(beta_t+(2*a+1)*beta_t2);
            c = c+n;
        end
        %---------
        
        
        save(pliki(b).logTC, 'LogTC');
        LogTC = 100*(exp(LogTC)-1);
        TC(:,1) = TC(:,1)+ mean(LogTC,2).*ile_danych(b,2);
        TC(:,2) = TC(:,2)+ std(LogTC')'.*ile_danych(b,2);
        clear LogTC;
        %------
        %LogATC = x_mid_av*beta_var;

        LogATC = x_mid_av*(beta_t+((t^2-1)/(t-1))*beta_t2); %Jesli tu jest zle
        %------
        
        save(pliki(b).logATC, 'LogATC');
        LogATC = 100*(exp(LogATC)-1);
        ATC(:,1) = ATC(:,1)+ mean(LogATC,2).*ile_danych(b,2);
        ATC(:,2) = ATC(:,2)+ std(LogATC')'.*ile_danych(b,2);
        clear LogATC;
        LogIC = zeros(dim(1),ile_danych(b,1),dim(3)+1); %"+1" poniewaz ostatni wymiar zawiera sume, czyli laczny wplyw dla IC
        %disp(size(eff));
        %logEC  = u_ti - u_(t+1)i; bo "-u_ti";
        
        %----------
        for i = 1:dim(3)
            c = 1;
            for a = 1:(t-1)
                LogIC(c:(c+n-1),:,i) = x_delta_comp(c:(c+n-1),:,i)*(beta_const+(a+0.5).*beta_t+(a^2+a+0.5).*beta_t2);
                LogIC(c:(c+n-1),:,dim(3)+1) = LogIC(c:(c+n-1),:,dim(3)+1)+LogIC(c:(c+n-1),:,i);
                c = c+n;
            end
        end
        %-----------
        
        save(pliki(b).logIC, 'LogIC');
        LogIC = 100*(exp(LogIC)-1);
        %disp(size(LogIC));
        %disp(size(ret.IC));
        for i =1:(dim(3)+1)
            IC(:,1,i) = IC(:,1,i)+ mean(LogIC(:,:,i),2).*ile_danych(b,2);
            IC(:,2,i) = IC(:,2,i)+ std(LogIC(:,:,i)')'.*ile_danych(b,2);
        end
        clear LogIC;
        LogAIC = zeros(n,ile_danych(b,1),dim(3)+1);
        %------------
        for i = 1:dim(3)
                LogAIC(:,:,i) = x_delta_av_comp(:,:,i)*(beta_const+(0.5*t+0.5).*beta_t+(0.5*t^2+0.5).*beta_t2);
                LogAIC(:,:,dim(3)+1) = LogAIC(:,:,dim(3)+1)+LogAIC(:,:,i);
        end
        save(pliki(b).logAIC, 'LogAIC');
        LogAIC = 100*(exp(LogAIC)-1);
        for i =1:(dim(3)+1)
            AIC(:,1,i) = AIC(:,1,i)+ mean(LogAIC(:,:,i),2).*ile_danych(b,2);
            AIC(:,2,i) = AIC(:,2,i)+ std(LogAIC(:,:,i)')'.*ile_danych(b,2);
        end
    end
    counter = b/l_plikow;
    msg = sprintf('Step 2/3: %d%%',round(counter*100));
    waitbar(counter,h,msg);
end
delete(h);
end


