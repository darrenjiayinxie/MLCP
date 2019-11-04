function F_B = normal_map(x,F,l)

z = projection_operator(x,l);

F_B = F(z) + x - z;


end