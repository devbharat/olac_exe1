function Task = Task_Design( )
%TASK_DESIGN Defines a task for a controller to execute

Task = struct;

Task.dt             = 0.02;        % sampling time period

Task.start_time     = 0;

Task.goal_time      = 10;          % Time to reach goal

Task.start_x        = [ 0; 0; 0;   % position x,y,z
                        0; 0; 0;   % roll, pitch, yaw
                        0; 0; 0;   % velocity x,y,z
                        0; 0; 0 ]; % angular rates roll, pitch, yaw
                    
Task.goal_x         = [10; 0; 0;   % position x,y,z
                        0; 0; 0;   % roll, pitch, yaw
                        0; 0; 0;   % velocity x,y,z
                        0; 0; 0 ]; % angular rates roll, pitch, yaw
                    
Task.vp1            = [ 5;  0; 0;  % via-point state to pass through
                        0; 0; 0;   % roll, pitch, yaw
                        0; 0; 0;   % velocity x,y,z
                        0; 0; 0 ]; % angular rates roll, pitch, yaw 
                    
Task.vp2 = [5;-5;-5;0;0;0;         % difficult waypoint only possible for
            0;0;0;0;0;0;];         % ILQC   
                    
Task.vp_time = Task.goal_time/3;   % time to pass through via-point
                 
Task.max_iteration  = 5;           % Maximum ILQC iterations  

Task.input_noise_mag = 0.186;      % Adds noise on input to simulation

Task.cost = [];                    % cost encodes performance criteria of 
                                   % how to execute the task. Filled later.
end

