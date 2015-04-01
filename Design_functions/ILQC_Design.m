function [Controller,cost] = ILQC_Design(Model,Task,Controller,Simulator)
% ILQC_Design Implements the iterative Linear Quadratic Controller (ILQC)
%
%    see script section 1.6 for a formal description of the algorithm.


%% Define functions that return the quadratic approximations of the cost 
% function at specific states and inputs.
% Example usage (Eq.(1.78)):
%     xn = [ x y z ... ]'; 
%     un = [ Fx Mx My Mz]';
%     t  = t;
%     Qm(xn,un) = Qm_fun(t,xn,un); 

% intermediate cost (l) quadratizations
l_ = Task.cost.l*Task.dt;
q_fun   = matlabFunction(  l_,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% dl/dx
l_x = jacobian(Task.cost.l,Task.cost.x)'*Task.dt; % cont -> discr. time
Qv_fun  = matlabFunction( l_x,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% ddl/dxdx
l_xx = jacobian(l_x,Task.cost.x);
Qm_fun  = matlabFunction(l_xx,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% dl/du
l_u = jacobian(Task.cost.l,Task.cost.u)'*Task.dt; % cont -> discr. time
Rv_fun  = matlabFunction( l_u,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% ddl/dudu
l_uu = jacobian(l_u,Task.cost.u);
Rm_fun  = matlabFunction(l_uu,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});
% ddl/dudx
l_xu = jacobian(l_x,Task.cost.u)';
Pm_fun  = matlabFunction(l_xu,'vars',{Task.cost.t,Task.cost.x,Task.cost.u});

% final cost (h) quadratizations
h_ = Task.cost.h;
qf_fun  = matlabFunction(  h_,'vars',{Task.cost.x});
% dh/dx
h_x = jacobian(Task.cost.h,Task.cost.x)';
Qvf_fun = matlabFunction( h_x,'vars',{Task.cost.x});
% ddh/dxdx
h_xx = jacobian(h_x,Task.cost.x);
Qmf_fun = matlabFunction(h_xx,'vars',{Task.cost.x});


n_x = length(Task.cost.x); % dimension of state space
n_u = length(Task.cost.u); % dimension of control input
n_t  = (Task.goal_time-Task.start_time)/Task.dt + 1; % number of time steps

% desired value function V* is of the form (Eq.(1.79))
% V*(dx,n) = s + dx'*Sv + 1/2*dx'*Sm*dx
s    = zeros(1  ,n_t);
Sv   = zeros(n_x,n_t);
Sm   = zeros(n_x,n_x,n_t);

duff = zeros(n_u,1  ,n_t-1);
K    = zeros(n_u,n_x,n_t-1);

% Shortcuts for function pointers to linearize systems dynamics:
% e.g. Model_Alin(x,u,Model_Param)
Model_Param = Model.param.syspar_vec;
Model_Alin  = Model.Alin{1}; 
Model_Blin  = Model.Blin{1}; 


% Each ILQC iteration approximates the cost function as quadratic around the
% current states and inputs and solves the problem using DP.
i  = 1;
while ( i <= Task.max_iteration && ( norm(squeeze(duff)) > 0.01 || i == 1 ))
     
    %% Problem 2.1.1: forward pass / "rollout" of the current policy
    % sim_out = ...
    
    cost(i) = Calculate_Cost(sim_out, q_fun, qf_fun);
    fprintf('Cost of Iteration %2d (metric: ILQC cost function!): %6.4f \n', i-1, cost(i));
    
    if ( i > 1 && cost(i) > 2*cost(i-1) )
        fprintf('It looks like the solution may be unstable. \n')
        fprintf('Press ctrl+c to interrupt iLQG, or any other key to continue. \n')
        pause
    end
    
    % define nominal state and control input trajectories
    % X0 = ...
    % U0 = ...
    % T0 = ...
    

      
    %% Problem 2.1.2: Solve Riccati-like equations backwards in time
    % Initialize the value function elements starting at final time step 
    % (Eq.(1.87)
    xf = X0(:,end); % final state when using current controller   
    % Sm(:,:,n_t) = ...
    % Sv(:,n_t)   = ...
    % s(n_t)      = ...
    
    % "Backward pass": Calculate the coefficients (s,Sv,Sm) for the value 
    % functions at earlier times by proceeding backwards in time (DP-approach)
    for n = (length(sim_out.t)-1):-1:1
        
        % state of system at time step n
        x0 = X0(:,n);
        u0 = U0(:,n);
              
        % Discretize and linearize continuous system dynamics Alin around
        % specific pair (xO,uO). See exercise sheet Eq.(4) for details
        % A = ...;
        % B = ...;
        

        % 2nd order approximation of cost function at time step n (Eq.(1.78))
        % q   = ...
        % Qv  = ...
        % Qm  = ...
        % Rv  = ...
        % Rm  = ...
        % Pm  = ...
        
        % control dependent terms of cost function (Eq.(1.81)) 
        % g = ...                    % linear control dependent
        % G = ...                    % control and state dependent
        % H = ...                    % quadratic control dependent
        
        % ensuring H remains symmetric
        H = (H+H')/2; % important, do not delete!
             
        % the optimal change of the input trajectory du = duff + K*dx (Eq.(1.82)) 
        % duff(:,:,n) = ...
        % K(:,:,n)    = ...
               
        % Solve Riccati-like equations for current time step n (Eq.(1.84)
        % Sm(:,:,n) = ...
        % Sv(:,n) = ...
        % s(n) = ...
              
    end % of backward pass for solving Riccati equation
    
    % define theta_ff in this function
    Controller.theta = Update_Controller(X0,U0,duff,K);
    
    i = i+1;
end

% simulating for the last update just to calculate the final cost
sim_out    = Simulator(Model,Task,Controller);
cost(i) = Calculate_Cost(sim_out, q_fun, qf_fun);
fprintf('Cost of Iteration %2d: %6.4f \n', i-1, cost(i));
end



function theta = Update_Controller(X0,U0,dUff,K)
% UPDATE_CONTROLLER Updates the controller after every ILQC iteration
% (Eq.(1.88))
%
%  X0  - state trajectory generated with controller from previous 
%        ILQC iteration.
%  UO  - control input generated from previous ILQC iteration.
%  dUff- optimal update of feedforward input found in current iteration
%  K   - optimal state feedback gain found in current iteration
%
%  The updated control policy has the following form:
%  U1 = U0 + dUff + K(X - X0)
%     = U0 + dUff - K*X0 + K*X
%     =      Uff         + K*x
%  
%  This must be brought into the form 
%  U1 = theta' * [1,x']   


% feedforward control input theta_ff = U0 + dUff - K*X0
% theta_ff = ...

% feedback gain of control input
theta_fb = permute(K,[2 1 3]);      

% puts below (adds matrices along first(=row) dimension)
theta = [theta_ff;        % size: (n_x+1) * n_u * n_t-1
         theta_fb];  

end


function cost = Calculate_Cost(sim_out, q_fun, qf_fun)
% calcules the cost of current state and input trajectory for the current
% ILQC cost function. Not neccessarily the same as the LQR cost function.

X0 = sim_out.x(:,1:end-1);
xf = sim_out.x(:,end);
U0 = sim_out.u;
T0 = sim_out.t(1:end-1);

cost = sum(q_fun(T0,X0,U0)) + qf_fun(xf);
end
