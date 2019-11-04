function [zk_n,wk_n,A] = MLCP(M,q,r,zk,wk,l)


A.zer_tol = 1e-16;
A.piv_tol = 1e-10;
A.minor_iter=1;


A.n = size(M,1);
A.zk = zk;
A.wk = wk;
A.r = r;
A.M = M;
A.q = q;
A.l = l;

A.Basis = zeros(A.n,A.n);
I_w = -eye(A.n);
I_a = eye(A.n);
A.table = zeros(A.n,2); % container for the index of basic variables and variable type
A.var = zeros(A.n,3);

A.z = A.zk;
A.w = A.wk;
A.s = 1/norm(A.r); % scale variable
N_c = 30;
%% Initialize the basis
for i = 1:A.n
    
    if (A.l(i)<A.zk(i)) % case 1
        if A.l(i) == 0
             A.table(i,2) = 1; % the type index for  bounded z_i is 1
        else
             A.table(i,2) = 0; % the type index for  z_i is 1
        end   
        A.table(i,1) = i;
        A.Basis(:,i) = A.M(:,i);
        A.var(i,1) = 1;
    elseif (A.l(i) == A.zk(i)) % case 2
        A.table(i,1) = i;
        A.table(i,2) = 2; % the index for w_i is 2
        A.Basis(:,i) = I_w(:,i);
        A.var(i,2) =1;
    elseif (A.l(i)>A.zk(i)) % error
        error('zk not in the box');
    end   
end


%% basic variables
i_M = find((A.table(:,2) <=1)); % index of the z variables in basis
i_w = find(A.table(:,2) ==2);
A.Vec_B = zeros(A.n,1);
A.Vec_B(i_M) = A.zk(i_M);
A.Vec_B(i_w) = A.wk(i_w);

%% Rank deficiency
Mzz = A.M(i_M,i_M);

while det(Mzz) == 0
    [~,idx]=ldcols(Mzz); %idx index of the columns that are l.d
    k_z = i_M(idx(1)); % index of one of the z variables in basis that are l.d
    i_table = find(A.table(:,1) == k_z); % table index for k_z
    if A.l(k_z) == A.zk(k_z)
        A.table(i_table,2) = 2; % z type switch to w type
        A.Basis(:,i_table) = I_w(:,k_z);
        A.var(i_table,1) = A.var(i_table,1) -1;
        A.var(i_table,2) = A.var(i_table,2) +1;
        A.Vec_B(i_table) = 0;
    else
        A.table(i_table,2) = 4; % z type switch to a type
        A.Basis(:,i_table) = I_a(:,k_z);
        A.var(i_table,1) = A.var(i_table,1) -1;
        A.var(i_table,3) = A.var(i_table,3) +1;
        A.z(k_z) = A.Vec_B(i_table);
        A.Vec_B(i_table) = 0;
    end
    Mzz(:,idx(1)) = [];
    Mzz(idx(1),:) = [];
    i_M(idx(1)) =[];
end




%% Start of the Pivot, let t be entering variable, and residual column r entering the basis
A.type_enter = 3; % t type
A.type_leave = 10; % unknown know
A.i_enter = A.n+1; % index of t type. i_var: index of entering variables
A.Be = -A.r*A.s;


%A.i_Basis = inv(A.Basis);


A.coeff = A.Basis\A.Be;
A = ratio_test(A); 
A.i_t = A.i_table;
A.t = A.theta;
if A.theta >= 1/A.s
    A.t = 1/A.s;
    A.theta = 1/A.s;
    A.Vec_B = A.Vec_B - A.theta*A.coeff;
    % update A.z
    i_z = A.table((A.table(:,2) <= 1),1);
    A.z(i_z,1) = A.Vec_B(i_z);
    zk_n = A.z;
    wk_n = A.w;
    return
end

A = updates(A);

A = pivot_rules(A);



%% Pivot loop
while (A.type_leave ~= 3)
    
    A = ratio_test(A); 
    A = updates(A);
    A = pivot_rules(A);
    
    % check cycling
    if (find(A.var(:,:)>N_c))
        error('cycling');
    end
    A.minor_iter = A.minor_iter+1;    
end
zk_n = A.z;
wk_n = A.w;
end