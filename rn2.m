%Kamil Makie�a
%Katedra Ekonometrii i Bada� Operacyjnych
%Uniwersytet Ekonomiczny w Krakowie

%funkcja przyjmuje jeden argument: wektor logarytm�w 10'ych warto�ci 
%funkcji wiarygodno�ci dla (zaakceptowanych) przebieg�w MCMC
function wynik = rn2(mdd)
%sortuje wektor ze zlogarytmowanymi wartosci funkcji wiarygodnosci
mdd = sort(mdd);
wynik = licz_rn2(mdd);
%miejsce na korekte Lenk'i, "in progress"...
return

%funkcja lokalna z kt�rej korzysta rn2
function w = licz_rn2(mdd)
rozstep = range(mdd);
if rozstep < 300
    %przesuni�cie o sta�� 
    mdd_shift = mean(mdd);
    mdd = mdd-mdd_shift;
    %licz� warto��
    w = log10(1/mean(1./power(10,mdd)))+mdd_shift;
else
     dim = size(mdd,1);
     a = round(dim/2);
     b = dim-a;
     %przesuniecie o sta��
     w_l = licz_rn2(mdd(1:a));
     w_h = licz_rn2(mdd(a+1:dim));
     sr = (w_l+w_h)/2;
     w_l = w_l-sr;
     w_h = w_h-sr;
     %wyznaczam warto�� w oparciu wz�r na �redni� harmoniczn� wa�on� oraz
     %dodaje sta��
     w = log10(1/(((a/dim)/power(10,w_l))+((b/dim)/power(10,w_h))))+sr;
end