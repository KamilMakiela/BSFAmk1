function u = kmdraw(mi,s,n)
flaga = 0;
while flaga == 0
    r=rand(n,1);
    u = mi+s.*norminv((r+(1-r).*normcdf(-mi./s,0,1)),0,1);
    if sum(isinf(u)) >= 1
        a=isinf(u);
        u(a==1) = ltruncn2(0,mi(a==1),s,0);
    end
    if isempty(find(u<0, 1))
        flaga = 1;
    end
end

function x=ltruncn2(a,mu,sigma, count)
count = count+1;
%fprintf('Poziom %4.1f \n', count);
%Generates a draw of univariate truncated normal random variable x,
%with mean mu and std dev sigma, where x>a
% lambda = c (a* u Ani)
n = size(mu);
c=(a-mu)./sigma;
z=-log(1-rand(n))./c;
check = rand(n) - exp(-0.5.*(z.^2));
arg = find(check>0);
if isempty(arg)
    z=z+c;
    x=mu+sigma*z;
else
    %fprintf('Poziom %4.1f \n', count);
    arg2 = find(check<=0);
    z(arg2) = z(arg2)+c(arg2);
    x(arg2) = mu(arg2)+sigma*z(arg2);
    x(arg) = ltruncn2(0,mu(arg),sigma, count);
end