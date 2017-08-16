function [beta, elast, skala, eff, lambda, s, mdd, SSE, ile_danych] = stats(folder, spalonych, model, H, x, y, vec, pliki, n,t)

[beta, elast, skala, eff, lambda, s, mdd, ile_danych] = stats_CDnew2(folder, spalonych, H, vec, pliki, n,t);


try
    SSE = ((x*beta(:,1) - y - eff.all(:,3))'*(x*beta(:,1) - y - eff.all(:,3)));
catch exception
    disp('SSE does not include inefficiency component');
    SSE = (x*beta(:,1) - y)'*(x*beta(:,1) - y);
    %SSE = 0;
end