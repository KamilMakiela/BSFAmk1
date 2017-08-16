%Prosty Gibbs bez restrykcji
function [fi, s, u_n, beta_n, mdd, przep] = GibbsVED1Simple(y, x, n, t, u, beta, par_rozk, st_fi, przep, a1, restr_from, restr, par, mediana)

%{
    u0 = 0.3;             %srednia nieefektywnosc (dla startu algorytmu)
    r0=0.75;              %mediana a priori
    n0=10^(-6);           %dla precyzji
    s0=10^(-6);
%}
u_n = u;
fi = zeros(2,1);
niesk = 1;

while niesk == 1;
    niesk = 0;
    %losuje fi
    %fi= randg(n*t+1)/(sum(u)-log(0.75));
    fi(1,1) = randg(1+n*t)/(sum(1*u.*(st_fi(2,1).^(par_rozk)))-mediana);   
    fi(2,1) = randg(1+sum(par_rozk))/(par_rozk'*(u.*fi(1,1)) + 1);
    %fi(2,1) = randg(1+sum(par_rozk))/(par_rozk'*(u.*fi(1,1)*fi(3,1)) + 1);
    
    
    %losuje s^2
    %s0 = 10^(-6) a priori
    %n0 = 10^(-6)
    s2= 1/(randg(0.5*(10^(-4)+ n*t))/(0.5*(10^(-4)+(y+u-x*beta)'*(y+u-x*beta))));
    s=sqrt(s2);
    
    %losuj beta
    %mi=(x'*x)\(x'*(y+u));
    %sigma=s*(chol(inv(x'*x))');
    %flaga = 0;
    
    
    %DLA ROZK£ADU REFERENCYJNEGO
    %mi = par.xtrx\(x'*(y+u));
    %sigma = s*par.cholaski;
    
    %DLA ROZK. W£AŒCIWEGO
    sigma2 = par.pr_diag/(par.pr_prec + par.xtrx/s2);
    mi = sigma2*(par.pr_precXex + x'*(y+u)./s2);
    sigma = chol(sigma2)';
    
    if isnan(restr)
        beta_n = mi + sigma*randn(par.x2dim,1);
    else
        if a1 <= restr_from
            beta_n = mi + sigma*randn(par.x2dim,1);
        else
            beta_n = mi+sigma*randn(par.x2dim,1);
            %disp(restr*beta_n);
            war = min(restr*beta_n);
            %if war < 0
            %    flaga = 1;
            %else
            %    flaga = 0;
            %end
            while war <= 0
                %disp(war);
                %disp(a1);
                beta_n = mi+sigma*randn(par.x2dim,1);
                war = min(restr*beta_n);
            end
        end
    end

    %losuje u
    try 
%        r = rand(n*t,1);
        lbd = fi(1,1).^(-1) .* fi(2,1).^(-par_rozk);
        mi = x*beta_n - y - s2./lbd;
        u_n = kmdraw(mi,s,n*t);
%         u_n = mi+s*norminv((r+(1-r).*normcdf(-mi./s,0,1)),0,1);
%         if sum(isinf(u_n)) >= 1
%             dupa;
%         end
    catch e;
        przep = przep +1;
        fprintf('Repeat: %6.1f at: %6.1f \n', przep, a1);
        niesk = 1;
    end

end
%liczenie skladnika gestosci a posteriori
%y_s = (y+u_n-x*(beta_n))./s;
%mdd = -sum(log((normpdf(y_s))));
%mdd = 1/prod(normpdf(y_s));
mdd = 1/prod(normpdf(y, x*(beta_n)-u_n, ones(size(y))*s));
return