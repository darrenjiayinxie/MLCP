%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x - initial guess
% l - lower bound
% F - functions 
% dF - differential of functions
% the path solver contains d-steps, m-steps, watchdog steps and projected
% gradient steps
function PATH_solver(x,l,F,dF)


%% Initialization 
z0 = projection_operator(x,l);
n_bar = 5; % For every n_bar steps, an m-steps are performed
delta = 1; % caritia for d-steps
beta = 0.5;
merit_tol = 1e-10;

k = 0;
check_point = 0;
best_point = 0;
j = 0;
b = 0;
delta0 = delta;
R0 = merit_function(z0,F,l);

%% major iteration
z = z0;

while merit_function(z,F,l) <= merit_tol
    
    
    
end

























end