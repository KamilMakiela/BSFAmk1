%Kamil Makie³a
%Katedra Ekonometrii i Badañ Operacyjnych
%Uniwersytet Ekonomiczny w Krakowie

%funkcja przyjmuje jeden argument: wektor logarytmów 10'ych wartoœci 
%funkcji wiarygodnoœci dla (zaakceptowanych) przebiegów MCMC
function wynik = rn2(mdd)
%sortuje wektor ze zlogarytmowanymi wartosci funkcji wiarygodnosci
mdd = sort(mdd);
wynik = licz_rn2(mdd);
%miejsce na korekte Lenk'i, "in progress"...
return

%funkcja lokalna z której korzysta rn2
function w = licz_rn2(mdd)
rozstep = range(mdd);
if rozstep < 300
    %przesuniêcie o sta³¹ 
    mdd_shift = mean(mdd);
    mdd = mdd-mdd_shift;
    %liczê wartoœæ
    w = log10(1/mean(1./power(10,mdd)))+mdd_shift;
else
     dim = size(mdd,1);
     a = round(dim/2);
     b = dim-a;
     %przesuniecie o sta³¹
     w_l = licz_rn2(mdd(1:a));
     w_h = licz_rn2(mdd(a+1:dim));
     sr = (w_l+w_h)/2;
     w_l = w_l-sr;
     w_h = w_h-sr;
     %wyznaczam wartoœæ w oparciu wzór na œredni¹ harmoniczn¹ wa¿on¹ oraz
     %dodaje sta³¹
     w = log10(1/(((a/dim)/power(10,w_l))+((b/dim)/power(10,w_h))))+sr;
end