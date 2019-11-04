function A = parameters(A)
N=24;
A.n = N;
A.wk = zeros(N,1);

global e_t e_o e_r mu;
e_t = 1;
e_o = 1;
e_r = 1;
mu =0.5;

% mass of the cuboid

global m g;
m = 5;
g = 9.8;

% dimensions of the cuboid
global len wid heg I_xx I_yy I_zz;
len = 1;  % in fixed body frame's x direction
wid = 1;  % in fixed body frame's y direction
heg = 1;  % in fixed body frame's z direction

I_xx = (m/12)*(wid^2+heg^2);
I_yy = (m/12)*(heg^2+len^2);
I_zz = (m/12)*(wid^2+len^2);
% 
% I_xx = 5;
% I_yy = 5;
% I_zz = 1;

% time-step length
global h;
h = 0.01;

% applied wrenches
global p_x p_y p_z p_xt p_yt p_zt;

p_x = 1;
p_y = 1;
p_z = 1;
p_xt = 1;
p_yt = 0;
p_zt = 0;

% q_old - position and orientation vector at l, q_old=[q_xo;q_yo;q_zo;q0_o;q1_o;q2_o;q3_o]
global q_old;

q_old = [0;0;heg/2;1;0;0;0];

% nu_old - generalized velocity vector at l, nu_old=[v_xo;v_yo;v_zo;w_xo;w_yo;w_zo]
global nu_old;

nu_old = [0.1;0;0;0;0;0];

%% defining the unknown variables 

% Nu - generalized valocity at l+1, Nu = [v_x;v_y;v_z;w_x;w_y;w_z]
nu =[0.1;0;0;0;0;0]; 

% ECP - pair of position vector of ECP at l+1, ECP = [a1_x;a1_y;a1_z;a2_x;a2_y;a2_z]
ECP = [0;0;0;0;0;0]; 

% Con_wrench - contact wrench at l+1, Con_wrench = [p_t;p_o;p_r]
Con_wrench = [1;0;0];

% La - lagrange multiplier at l+1, La =[l1;l2;l3;l4;l5;l6;l7]
La =[0;0;0;0;0;0;0]; % at least provide one non-zero value

% p_n,sig - normal impulse and contact velocity at l+1
p_n =  g*m*h;
sig = 0;

% Z - initial guess total unknown variables
A.zk = [nu;ECP;Con_wrench;sig;La;p_n]; 



%% defining the infinity constant, lower bound and upper bound


% n_1 - size of unbounded variables
n_1 = size([nu;ECP;Con_wrench]);

% n_2 - size of bounded variables
n_2 = size([sig;La;p_n]); 

% l - lower bound 
A.l(1:n_1,1) = -Inf; 
A.l(n_1+1:n_1+n_2,1) = 0;

A.u(1:n_1+n_2,1) = Inf;

f= @fcn;
Df = @jacobian;

A.M =Df(A.zk);
A.q = f(A.zk) - A.M*A.zk;

A.r = -(A.M*A.zk-eye(N)*A.wk+A.q);