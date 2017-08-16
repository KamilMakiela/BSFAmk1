%Prosty Gibbs bez restrykcji
function [fi_n, s_n, u_n, beta_n, mdd, przep] = GibbsVED2Alt(y, x, n, t, u, beta, par_rozk, fi, s, przep, a1, restr_from, restr, par, mediana)
% disp(beta);
% disp(fi);
% disp(s);
% disp('pass');
%{
    u0 = 0.3;             %srednia nieefektywnosc (dla startu algorytmu)
    r0=0.75;              %mediana a priori
    n0=10^(-6);           %dla precyzji
    s0=10^(-6);
%}
u_n = u;
fi_n = zeros(3,1);
niesk = 1;
while niesk == 1;
    niesk = 0;
    
    %losuje u
%     try 
        %r = rand(n*t,1);
        lbd = fi(1,1)^(-1).* fi(2,1).^(-par_rozk(:,1)).* fi(3,1).^(-par_rozk(:,2));
        mi = x*beta - y - (s^2)./lbd;
        u_n = kmdraw(mi,s,n*t);
%         u_n = mi+s*norminv((r+(1-r).*normcdf(-mi./s,0,1)),0,1);
%         if sum(isinf(u_n)) >= 1
%             dupa;
% %         end
%     catch e;
%         przep = przep +1;
%         u_n = u;
%         fprintf('Powtorka: %6.1f przy: %6.1f \n', przep, a1);
%         niesk = 1;
%     end
    
    %losuje fi
    %fi= randg(n*t+1)/(sum(u)-log(0.75));
    fi_n(1,1) = randg(1+n*t)/(sum(1*u_n.*(fi(2,1).^(par_rozk(:,1))).*(fi(3,1).^(par_rozk(:,2))))-mediana);    
    fi_n(2,1) = randg(1+sum(par_rozk(:,1)))/(sum(par_rozk(:,1).*(u_n.*fi_n(1,1).*(fi(3,1).^(par_rozk(:,2))))) + 1);
    fi_n(3,1) = randg(1+sum(par_rozk(:,2)))/(sum(par_rozk(:,2).*(u_n.*fi_n(1,1).*(fi_n(2,1).^(par_rozk(:,1))))) + 1);
    
    %losuje s^2
    %s0 = 10^(-4) a priori
    %n0 = 10^(-4)
    s2_n= 1/(randg(0.5*(10^(-4)+ n*t))/(0.5*(10^(-4)+(y+u_n-x*beta)'*(y+u_n-x*beta))));
    s_n=sqrt(s2_n);
    
    %losuj beta
    %mi=(x'*x)\(x'*(y+u));
    %sigma=s*(chol(inv(x'*x))');
    %flaga = 0;
    
    sigma2 = par.pr_diag/(par.pr_prec + par.xtrx/s2_n);
    mi = sigma2*(par.pr_precXex + x'*(y+u)./s2_n);
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
                if war > min(min(restr*beta))
                    beta = beta_n;
                    fi = fi_n;
                    s = s_n;
                    
                    %s2 = s2_n;
                    %u = u_n;
                    %disp(war);
                end
                %disp(war);
                %disp(a1);
            else
                beta = beta_n;
                fi = fi_n;
                s = s_n;
                %disp(war);
            end
        end
    end

end
%liczenie skladnika gestosci a posteriori
%y_s = (y+u_n-x*(beta_n))./s_n;
%mdd = -sum(log((normpdf(y_s))));
%mdd = 1/prod(normpdf(y_s));
mdd = 1/prod(normpdf(y, x*(beta_n)-u_n, ones(size(y))*s_n));
return