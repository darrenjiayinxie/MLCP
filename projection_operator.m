function z = projection_operator(x,l)
n = size(x,1);
z = zeros(n,1);
for i = 1:n
    if x(i) >= l(i)
        z(i) = x(i);
    else
        z(i) = l(i);
    end
    
    
end

end