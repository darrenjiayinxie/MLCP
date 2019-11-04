function Psi = merit_function(F,z,l)
Phi = Phi_operator(F,z,l);
Psi = 0.5*(Phi'*Phi);
end