function [y,z] = vect2mat(x1,x2,a,b)
y = zeros(a,b);
help = zeros(a,b);
z = zeros(2*a,b);
c = size(x1);

k=1;
for i = 1:b
    y(:,i) = x1(k:(k+a-1));
    help(:,i) = x2(k:(k+a-1));
    k = k+a;
end

k=1;
for i=1:a
    z(k,:) = y(i,:);
    z(k+1,:) = help(i,:);
    k = k+2;
end
return