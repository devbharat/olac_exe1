%% OLCAR - Exercise 1 - ILQC
close all; clearvars; clc;
       
addpath(genpath(pwd)); % adds folders and subfolders to path
% To generate plots of LQR/ILQG rollout
plotvec = {'quad_pos_noLoad','quad_angles_noLoad','control_input_noLoad','rotor_thrust_noLoad'};
create_pdf = [0 0 0 0 ];    % for which plots should a pdf be created
plot_ind   = [1 2 3 4 ];    % which data to plot on screen


%% Task definition 
Task = Task_Design();

%Load the dynamic model of the quadcopter
load('Quadrotor_Model.mat','Model'); % save as structure "Model" 

% Define cost function
Task.cost = Cost_Design( Model.param.mQ, Task );



%% Problem 1: Initial controller design 
[Initial_Controller, Cost_LQR] = LQR_Design(Model, Task);

% Visualization of LQR controller
sim_out_lqr = Quad_Simulator(Model,Task,Initial_Controller);
disp('LQR controller performance:');
fprintf('Cost with LQR controller (metric: LQR cost function!): J* = %.3f \n', Cost_LQR);
fprintf('Start Quadcopter position: x = %.3f, y = %.3f, z = %.3f \n', sim_out_lqr.x(1:3,1));
fprintf('Final Quadcopter position: x = %.3f, y = %.3f, z = %.3f \n\n', sim_out_lqr.x(1:3,end));
Visualize(sim_out_lqr,Model.param,'plot_mode',1);
% Plot_Result(sim_out_lqr,Model.param,'plots',plotvec(plot_ind),'file',create_pdf(plot_ind),'path',pwd)


%% comment out to proceed to Problem 2
return;


%% Problem 2: ILQC controller design 
t_cpu = cputime; 
[ILQC_Controller, Cost] = ILQC_Design(Model,Task,Initial_Controller,@Quad_Simulator);
t_cpu = cputime - t_cpu;

% Visualization of ILQC controller
sim_out_ilqc = Quad_Simulator(Model,Task,ILQC_Controller);
fprintf('The ILQC algorithm found a solution in %fs \n\n',t_cpu);
fprintf('Final Quadcopter Position: xQ = %.3f, yQ = %.3f, zQ = %.3f \n',sim_out_ilqc.x(1:3,end));
fprintf('Final Quadcopter Velocity: xQ = %.3f, yQ = %.3f, zQ = %.3f \n',sim_out_ilqc.x(7:9,end));
Visualize(sim_out_ilqc,Model.param,'plot_mode',1);
% Plot_Result(sim_out_ilqc,Model.param,'plots',plotvec(plot_ind),'file',create_pdf(plot_ind),'path',pwd)




