function A = pivot_rules(A)

%% determining entering type
if A.type_leave == 1 % z type
    A.type_enter = 2;
elseif A.type_leave == 2 % w type 
    A.type_enter = 1;
elseif A.type_leave == 3 % t type %% solution found
    return
elseif A.type_leave == 4 % a type
    if A.l(A.i_enter)== 0
        A.type_enter = 1;
    else
        A.type_enter = 0;
    end
end


%% determining the covering vector
if A.type_enter <= 1
    A.Be = A.M(:,A.i_enter);
    A.var(A.i_enter,1) = A.var(A.i_enter,1)+1;

elseif A.type_enter == 2
    I_w = -eye(A.n);
    A.Be = I_w(:,A.i_enter);
    A.var(A.i_enter,2) = A.var(A.i_enter,2)+1;
elseif A.type_enter == 4
    I_a = eye(A.n);
    A.Be = I_a(:,A.i_enter);
    A.var(A.i_enter,3) = A.var(A.i_enter,3)+1;
end


A.coeff = A.Basis\A.Be;

end