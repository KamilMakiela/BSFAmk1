function [mdd, sr, sr_kor, sr10, zer] = gestosc(str, l_plikow, poczatek, koniec)
pliki = struct([]);
for a = 1:l_plikow
    if a < 10
        pliki(a).mdd = strcat(str,'\','mdd0000',num2str(a),'.mat');
    else 
        pliki(a).mdd = strcat(str,'\','mdd000',num2str(a),'.mat');
    %disp(pliki(a).mdd);
    end
end

var(1:l_plikow) = struct('mdd_d', 0);
for b = 1:l_plikow
    var(b) = load(pliki(b).mdd,'mdd_d');
end
mdd = [var(:).mdd_d]';
mdd = mdd(poczatek:koniec);
sr = 1/mean(mdd);
zer = sum(mdd==0);
mdd(mdd==0) = NaN;
sr_kor = 1/nanmean(mdd);
sr10 = 1/trimmean(mdd,5);