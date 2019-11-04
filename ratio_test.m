%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M_B(z_B-t*X_z)+I_wB(w_B-t*X_w)+I_aB(a_B-t*X_w) 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = ratio_test(A)

x = A.Vec_B;
%% determine the moving direction when entering variables as z type not on boundary
A.direction = 1; % direction by default
flag = 0;
% idenitify 
if (A.type_leave == 4)&&(A.type_enter <= 1)
    flag = 1;
end
% determining the directions
if flag == 1
    jd = find(A.table(:,2)==3&(-A.coeff<-A.piv_tol)); % first test in opposite direction
    if isempty(jd)
        A.direction = 1; % t does not increase
        theta_z = Inf;
    else
        A.direction = -1;
        if A.type_enter == 1
            theta_z = A.z(A.i_enter);
        else
            theta_z = Inf; 
        end
    end
else
    theta_z = Inf; 
end

d = A.direction*A.coeff;
%% identify leaving variables as z and w type

j1=find((d>=A.piv_tol)&(A.table(:,2)==1));
j2=find((d>=A.piv_tol)&(A.table(:,2)==2));
j12 = [j1;j2];
if isempty(j12)
    theta12 = Inf;
else
    theta12=min((x(j12)+A.zer_tol)./d(j12)); % minimum ratio
end
j12=j12((x(j12)./d(j12))<=theta12); 
%% identify leaving variables as a type
j4 = find((abs(d)>=A.piv_tol)&(A.table(:,2)==4));
if isempty(j4)
    theta4 = Inf;
else
    theta4 = 0;
end
%j4=j4((x(j4)./d(j4))<=(theta4)); 
%% identify leaving variables as t type
j3 = find(A.table(:,2)==3&(d<=-A.piv_tol));
if isempty(j3)
    theta3 = Inf;
else
    theta3 = (x(j3)+A.zer_tol/A.s-1/A.s)/d(j3);
end

%% identify the switch position in table according to priority table 
Theta = [theta3,theta4,theta12];

[theta,index] = min(Theta);
A.theta = theta;
if (theta == Inf)||((flag == 1)&&(theta>=theta_z))
    error('ray termination');
else 
    if index(1) ==1 % t type
        A.i_table = j3;
    elseif index(1) == 2 % a type
        A.i_table = j4(1); % choose the first one
    elseif index(1) == 3 % w,z type
        A.i_table = j12(1); % choose the first one
    end
 
end

end