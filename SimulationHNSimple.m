
function  [wyniki, pliki] = SimulationHNSimple(y, x, n, t, u, beta, iteracje, hand, objetosc, spalonych, restrykcje, mediana, priors)
t1 = clock;

%tworze folder z losowaniami - czas startu badania
przyjetych = iteracje - spalonych;
if (przyjetych <= 0)
    disp('Check the number of cycles!');
    return;
end
wyniki = datestr(now);
wyniki = strrep(wyniki,':','-');
if ~isdir(wyniki)
    mkdir(wyniki);
    disp('Folder with results has been created. Its name is the current date and time of the simulation. ');
end
%wyznaczam skok dla waitbara
skok_bar = round(iteracje/100);
[r, c] = size(x);
h = waitbar(0,'Running, this may take some time...',...
            'Name','MCMC',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
assignin('base', 'Waitbar_handle', h);
      
pozycja = [4, 30, 270, 70];
set(h,'Position', pozycja);

if iteracje <= objetosc
    %to odpala kiedy tylko jeden plik do zapisu
    b_B = zeros(c,iteracje);
    u_B = zeros(r,iteracje);
    fi_B = zeros(1,iteracje);
    s_B = zeros(1,iteracje);    %uzmienniajac moge uogolnic cala procedure losowania
    mdd_d = zeros(1,iteracje);
    liczba_plikow = 1;
    dlugosc_mat = iteracje;
else
    b_B = zeros(c,objetosc);
    u_B = zeros(r,objetosc);
    fi_B = zeros(1,objetosc);
    s_B = zeros(1,objetosc);    %uzmienniajac moge uogolnic cala procedure losowania
    mdd_d = zeros(1,objetosc);
    %wyznaczam ilosc przejsc petli zewnetrznej
    liczba_plikow = floor(iteracje/objetosc);
    if rem(iteracje,objetosc) ~= 0
        liczba_plikow = liczba_plikow + 1; %jesli cos zostalo to dodaje cykl
    end
    dlugosc_mat = objetosc;
end
    
disp(liczba_plikow);
disp(dlugosc_mat);
%wyznaczanie nazw plikow do zapisu
pliki = struct([]);
for a2 = 1:liczba_plikow
    str = strrep(sprintf('%2.3f', a2/1000),'.','0');
    
    pliki(a2).bety = strcat(wyniki,'\','bety',str);
    pliki(a2).u = strcat(wyniki,'\','u',str);
    pliki(a2).fi = strcat(wyniki,'\','fi',str);
    pliki(a2).s = strcat(wyniki,'\','s',str);
    pliki(a2).mdd = strcat(wyniki,'\','mdd',str);
end
%koniec wyznaczania nazw


%a1 to laczna liczba wszystkich iteracji
a1 = 0;
%rest_from wyznacza od kad ida restrykcje na parametry
restr_from = round(spalonych/2);
restr = restrykcje.restrict;
%sta?e parametry symulacji wyprowadzone przed petle
%par.cholaski = (chol(inv(x'*x))');
par.xtrx = (x'*x);
par.x2dim = size(x,2);


par.pr_ex = priors.ex;
par.pr_prec = priors.diag/priors.cov;
par.pr_diag = priors.diag;
par.pr_precXex = par.pr_prec*priors.ex;

control = 1;
przepelnien = 0;
for a2 = 1:liczba_plikow
    for a3 = 1:dlugosc_mat
        a1 = a1 + 1;
        [fi_B(1,a3), s_B(1,a3) , u, beta, mdd_d(1,a3), przepelnien] = GibbsHNSimple(y, x, n, t, u, beta, przepelnien, a1, restr_from, restr, par, mediana);
        b_B(:,a3) = beta;
        u_B(:,a3) = u;
        if (rem(a1,skok_bar)==0)
            counter = a1/iteracje;
            msg = sprintf('%d%%    cycles: %d',round(counter*100), a1);
            waitbar(counter,h,msg);
            %robie wypasny wykres
            plot(hand, b_B(c,1:a3));
            title(hand,'Sequential plot: the intercept');
        end
        if  getappdata(h,'canceling')
            control = 0;
            break
        end
    end
    if control == 0
        break
    end
    if isdir(wyniki)
        disp('I am adding new files to the results folder');
    else % ten warunek nie zostanie wykonany nigdy
        mkdir(wyniki);
        disp('Results folder has been created, which is odd...');
    end
    %cd(wyniki);
    %str = sprintf('%2.3f', a2/1000);
    %str = strrep(str,'.','0');
    save(pliki(a2).bety, 'b_B');
    save(pliki(a2).u, 'u_B');
    save(pliki(a2).fi, 'fi_B');
    save(pliki(a2).s,'s_B');
    save(pliki(a2).mdd,'mdd_d');
    
    %save(strcat(wyniki,'\',str), 'b_B', 'u_B', 'fi_B', 's_B', 'mdd_d');
    if a2 == (liczba_plikow - 1)
        dlugosc_mat = iteracje - a2*dlugosc_mat;
    end
    clear b_B u_B fi_B s_B mdd_d;
    b_B = zeros(c,dlugosc_mat);
    u_B = zeros(r,dlugosc_mat);
    fi_B = zeros(1,dlugosc_mat);
    s_B = zeros(1,dlugosc_mat);
    mdd_d = zeros(1,dlugosc_mat);
end 
delete(h);
disp(etime(clock,t1)/60);
if(control == 0)
    wyniki = 'lipa';
end

return;