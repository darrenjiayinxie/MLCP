function Phi = Phi_operator(F,z,l)
n = size(z,1);
Phi = zeros(n,1);
Fz = F(z);
for i = 1:n
    if l(i) == 0
        Phi(i) = FB_function(z(i)-l(i),Fz(i));
    else
        Phi(i) = -Fz(i);
    end
end



end