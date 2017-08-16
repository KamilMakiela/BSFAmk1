%Prosty Gibbs bez restrykcji
function [fi_n, s_n, u_n, beta_n, mdd, przep] = GibbsHNAlt(y, x, n, t, u, beta, fi, s, przep, a1, restr_from, restr, par, mediana)

%{
    u0 = 0.3;             %srednia nieefektywnosc (dla startu algorytmu)
    r0=0.75;              %mediana a priori
    n0=10^(-6);           %dla precyzji
    s0=10^(-6);
%}
%u_n = u;
nt=n*t;
s2 = s^2;
omega_kw = 1/fi;

niesk = 1;
while niesk == 1;
    niesk = 0;
    
    %losuje u
    mi = omega_kw*(x*beta-y)/(omega_kw+s2);
    s_u = sqrt((s2*omega_kw)/(s2+omega_kw));
    u_n = kmdraw(mi,s_u,nt);
    
    %losuje fi
    fi_n= randg(0.5*nt+5)/(0.5*sum(u_n.^2)+10*(log(mediana)).^2); % precyzja rozkladu polnormalnego
    
    %losuje s^2
    %s0 = 10^(-4) a priori
    %n0 = 10^(-4)
    s2_n= 1/(randg(0.5*(10^(-4)+ nt))/(0.5*(10^(-4)+(y+u_n-x*beta)'*(y+u_n-x*beta))));
    s_n=sqrt(s2_n);
    
    %BETA
    
    %DLA ROZK£ADU REFERENCYJNEGO
    %mi = par.xtrx\(x'*(y+u_n));
    %sigma = s_n*par.cholaski;
    
    %DLA ROZK. W£AŒCIWEGO
    sigma2 = par.pr_diag/(par.pr_prec + par.xtrx/s2_n);
    mi = sigma2*(par.pr_precXex + x'*(y+u_n)./s2_n);
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
            if war <= 0
                niesk =1;
                if war > min(restr*beta)
                    beta = beta_n;
                    fi = fi_n;
                    s = s_n;
                    s2 = s2_n;
                    omega_kw = 1/fi_n;
                    %u = u_n;
                    %disp(war);
                end
                %disp(war);
                %disp(a1);
            else
                beta = beta_n;
                fi = fi_n;
                s = s_n;
                s2 = s2_n;
                omega_kw = 1/fi_n;
            end
        end
    end



end
%liczenie skladnika gestosci a posteriori
%y_s = (y+u_n-x*(beta_n))./s_n;
%mdd = -sum(log((normpdf(y_s))));
%mdd = 1/prod(normpdf(y_s));
%mdd = 1/prod(normpdf(y, x*(beta_n)-u_n, ones(size(y))*s_n));
mdd = sum(log10(normpdf(y, x*(beta_n)-u_n, ones(size(y))*s_n)));
return

