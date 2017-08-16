delta = 0.001;
x_kl = exp(user.xRaw);
x_kdelta_l = exp(user.xRaw);
x_kdelta_l(:,1) = x_kdelta_l(:,1) + delta;

x_k2delta_l = exp(user.xRaw);
x_k2delta_l(:,1) = x_k2delta_l(:,1) + 2*delta;

x_k_ldelta = exp(user.xRaw);
x_k_ldelta(:,2) = x_k_ldelta(:,2)+delta;

x_k_l2delta = exp(user.xRaw);
x_k_l2delta(:,2) = x_k_l2delta(:,2)+2*delta;

x_kdelta_ldelta = exp(user.xRaw);
x_kdelta_ldelta = x_kdelta_ldelta + delta;

[~, ~, ~, ~, arg_x_kl, ~, ~, ~, ~, ~, ~, ~] = PrepareData(user.y, log(x_kl), 30, 15, user.xlabel, 11, user.war);
[~, ~, ~, ~, arg_x_kdelta_l, ~, ~, ~, ~, ~, ~, ~] = PrepareData(user.y, log(x_kdelta_l), 30, 15, user.xlabel, 11, user.war);
[~, ~, ~, ~, arg_x_k2delta_l, ~, ~, ~, ~, ~, ~, ~] = PrepareData(user.y, log(x_k2delta_l), 30, 15, user.xlabel, 11, user.war);
[~, ~, ~, ~, arg_x_k_ldelta, ~, ~, ~, ~, ~, ~, ~] = PrepareData(user.y, log(x_k_ldelta), 30, 15, user.xlabel, 11, user.war);
[~, ~, ~, ~, arg_x_k_l2delta, ~, ~, ~, ~, ~, ~, ~] = PrepareData(user.y, log(x_k_l2delta), 30, 15, user.xlabel, 11, user.war);
[~, ~, ~, ~, arg_x_kdelta_ldelta, ~, ~, ~, ~, ~, ~, ~] = PrepareData(user.y, log(x_kdelta_ldelta), 30, 15, user.xlabel, 11, user.war);

f1 = exp(arg_x_kl*b_B).*(arg_x_kdelta_l*b_B-arg_x_kl*b_B)./delta;
f2 = exp(arg_x_kl*b_B).*(arg_x_k_ldelta*b_B-arg_x_kl*b_B)./delta;
f3 = exp(arg_x_kdelta_l*b_B).*(arg_x_k2delta_l*b_B-arg_x_kdelta_l*b_B)./delta;
f4 = exp(arg_x_k_ldelta*b_B).*(arg_x_k_l2delta*b_B-arg_x_k_ldelta*b_B)./delta;
%f5 = exp(arg_x_kdelta_l*b_B).*(arg_x_kdelta_ldelta*b_B-arg_x_kdelta_l*b_B)./delta;
%f6 = exp(arg_x_k_ldelta*b_B).*(arg_x_kdelta_ldelta*b_B-arg_x_k_ldelta*b_B)./delta;

f_kk = (f3-f1)./delta;
f_ll = (f4-f2)./delta;
f_kl = exp(arg_x_kl*b_B).*(arg_x_kdelta_ldelta*b_B-arg_x_kl*b_B)./(delta^2);
%f_lk = (f6-f2)./delta;

for i=1:450
    elast(i,:) = (x_kl(i,1).*f1(i,:)+ x_kl(i,2).*f2(i,:)).*f1(i,:).*f2(i,:)./...
        (x_kl(i,1)*x_kl(i,2)*(2*f_kl(i,:).*f1(i,:).*f2(i,:) - f_kk(i,:).*f2(i,:).^2 - f_ll(i,:).*f1(i,:).^2));
end

a(:,1) = mean(f1,2);
a(:,2) = mean(f2,2);
a(:,3) = mean(f3,2);
a(:,4) = mean(f4,2);
a(:,5) = mean(f_kk,2);
a(:,6) = mean(f_ll,2);
%a(:,7) = mean(f_kl,2);
a(:,7) = mean(f_lk,2);

mean(elast,2);