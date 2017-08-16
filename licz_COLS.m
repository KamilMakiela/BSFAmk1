function [elast_C, b_C, u_C, sse_C, e]= licz_COLS(y, x, vecs)


dim = size(x);
b_C = (x'*x)\(x'*y);
e = y - x*b_C;
e_max = max(e);
u_C = e_max - e;
b_C(dim(2),1) = b_C(dim(2),1) + e_max;
sse_C = e'*e;

elast_C = 0;
wymiar = size(vecs);
%disp(wymiar);
%disp(size(b_C));
if wymiar(1) > 1
	elast_C = zeros(wymiar(1),wymiar(3));
    for a = 1:wymiar(3)
        elast_C(:,a) = vecs(:,:,a)*b_C;
    end
end
return