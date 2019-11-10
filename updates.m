function A = updates(A)
%% identifying
% identify the leaving type
A.type_leave = A.table(A.i_table,2);
% indetify the i_leave
A.i_leave = A.table(A.i_table,1); % index of the leaving variable

%% updating 
% update basic variables and Basis

if A.direction == -1
    A.Vec_B = A.Vec_B + A.theta*A.coeff;
    if A.type_enter <=1
        A.Vec_B(A.i_table) = A.z(A.i_enter)-A.theta; 
    elseif A.type_enter ==2
        A.Vec_B(A.i_table) = A.w(A.i_enter)-A.theta; 
    end        
else
    A.Vec_B = A.Vec_B - A.theta*A.coeff;
    if A.type_enter == 3
        A.Vec_B(A.i_table) =  A.theta; 
    elseif A.type_enter <=1
        A.Vec_B(A.i_table) = A.z(A.i_enter)+ A.theta; 
    elseif A.type_enter == 2
        A.Vec_B(A.i_table) = A.w(A.i_enter)+ A.theta; 
    end
    
end
% u = A.Be-A.Basis(:,A.i_table);
% v = zeros(A.n,1);
% v(A.i_table) = 1;
% 
% A.i_Basis = A.i_Basis - (A.i_Basis*u*v'*A.i_Basis)/(1+v'*A.i_Basis*u);

A.Basis(:,A.i_table) = A.Be;

% update the table
A.table(A.i_table,2) = A.type_enter;
A.table(A.i_table,1) = A.i_enter;

% update the i_enter for next pivot
A.i_enter = A.i_leave;
%update the var
if A.Vec_B(A.i_t) > A.t+1e-6
    A.var = zeros(A.n,3);
    A.t= A.Vec_B(A.i_t);
end


% update A.z
i_z = A.table((A.table(:,2) <=1),1);
i_w = A.table((A.table(:,2) ==2),1);
i_zB = A.table(:,2) <=1;
i_wB = A.table(:,2) ==2;
A.z(i_z,1) = A.Vec_B(i_zB);
A.z(i_w,1) = 0;
A.w(i_w) = A.Vec_B(i_wB);
A.w(i_z,1) = 0;
end