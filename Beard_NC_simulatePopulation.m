function [ X1, X2, X3, X4] = Beard_NC_simulatePopulation(expt, Par_vary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code adapted from crosst.m of Huber 2011 (PMID: 21364572) 
%%% and regulation.m of Huber 2012 (PMID: 22218564). 
%%% Model originally developed in Beard 2005 (PMID: 16163394)
%%% This code simulates a population for one expt type
%%% - Runs multiple times, varying parameters by +/- X% each time
%%% Can simulate all drug additions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define model parameters
xpar = define_model_parameters;
otherpar = define_other_parameters;

%%% Define simulation time-frames and other time settings
printRequest = 0;   %Set to 1 if you want to print info on time settings
[t_prior, t_start,t_final,t_no_time,stepsize,tt,time]...
    = defineSimulationTimeFrames(printRequest);

protons_per_ATPase = 3;     % three protons per ATP synthase (for output)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set what figures to plot
plot_ss = 0;    % Set to 1 to plot steady-state simulations
plotStateVar = 1;  % Set to 1 to plot all raw state variables
plotStateVar_Choice = [4 8 13 19 23 24];   % Choose which outputs you want to plot (listed in getLegends.m)
plotStateVarFC = 0;  %...state variable fold changes
plotStateVarFC_Choice = [4 8 13 19 23 24]; 
plotOutput = 1;       % ...raw outputs  
plotOutput_Choice = [1 2 3 4 5 6 13 15 17 27];   
plotOutputFC = 0;  %...output fold changes
plotOutputFC_Choice = [1 2 3 4 5 6 13 15 17 27]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ctot0       = xpar(2);           % Cyt-c from Beard. Total IMS Cyt-c, Cred+Cox, molar
Qtot0       = xpar(3);           % total IMS Ubiquinol, Q+QH2, molar
ADTP_tot    = xpar(4);           % total Adenosine phosphates in cytosol, calculated from Evans_77

%%% Basal Respiratory State
state_fact  = 10/11; 
ATP_e       = state_fact*ADTP_tot;  % ATP level
ADP_e       = ADTP_tot-ATP_e;       % ADP level
fprintf('Respiratory State (ATP:ADP): %0.1i\n', state_fact)
fprintf('Glycolytic Capacity (K_ADTP_dyn): %0.1i\n', otherpar(4))

xo_single_cell  = initial(ADP_e, ATP_e, Ctot0, Qtot0); 

%%% ODE options
options = odeset('RelTol',1e-5, 'AbsTol',1e-8, 'MaxStep',10e-1, ...
    'InitialStep',1e-1, 'MaxOrder',5, 'BDF','on');


%%% Set seed of random number generator (for reproducible simulations)
setRNG = 1;
fprintf('\n********************')
fprintf('\n********************')
if setRNG == 1
    rng(1)
    fprintf('\nInitialised random number generator seed [rng(1)]')
end

%%% Set number of simulations (= number of cells in population)
numSims = 50;
printSim = 1;   % Set to 1 to output simulation numbers to screen
fprintf('\n%i simulations', numSims)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set conditions to simulate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define experiment to simulate (based on input, expt)
% 1: Rotenone, Oligo
% 2: AA, Oligo
% 3: Seahorse (Oligo, FCCP, Rot+AA)
% 4: Increased energy demand (IncEnDem, Oligo, FCCP, Rot+AA) 
% 5: FCCP, Rotenone 
printRequest = 1; %Set to 1 if you want to print expt info
fprintf('\nRunning expt %i', expt)
[rotenone,AA,oligo,CIV,FCCP,energy] = defineExptsToSimulate(expt,t_no_time,printRequest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set impairments to simulate
% Set parameters to 1 for 100% condition, and to <1 to
% simulate impairment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_xpar7 = Par_vary(1);%0.5;%0.42;%    % x_DH impairment 
v_xpar9 = Par_vary(2);%2e-1;%9e-2;%4.5e-2;%      % x_C1  
v_xpar10 = Par_vary(3);%1.2e-1;%7.2e-2;%5.2e-2;%     % x_C3
v_xpar13 = Par_vary(4);%6e-3;%3.2e-3%2.2e-3;%     % x_C4 impairment
v_xpar15 = Par_vary(5);%1e-4;%6e-5%4.2e-5;%    % x_F1 impairment
v_xpar21 = Par_vary(6);%1.5;%2.0;%     % x_Hle  
v_otherpar4 = Par_vary(7);%0.3;%0.82;%  % K_ADTP_dyn impairment
v_otherpar6 = Par_vary(8); % K_ADTP_cons (redundant)
v_xpar41 = Par_vary(9); % ROS Scavenging  


if v_xpar7 ~= 1, fprintf('Impairing x_DH [xpar(7)=%0.2i]\n',v_xpar7), end
if v_xpar9 ~= 1, fprintf('Impairing x_C1 [xpar(9)=%0.2i]\n',v_xpar9), end
if v_xpar10 ~= 1, fprintf('Impairing x_C3 [xpar(10)=%0.2i]\n',v_xpar10), end
if v_xpar13 ~= 1, fprintf('Impairing x_C4 [xpar(13)=%0.2i]\n',v_xpar13), end
if v_xpar15 ~= 1, fprintf('Impairing x_F1 [xpar(15)=%0.2i]\n',v_xpar15), end
if v_xpar21 ~= 1, fprintf('Impairing x_Hle [xpar(21)=%0.2i]\n',v_xpar21), end
if v_otherpar4 ~= 1, fprintf('Impairing K_ADTP_dyn [otherpar(4)=%0.2i]\n',v_otherpar4), end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run multiple simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xparAll, stateVarAll,stateVarFCAll,stateVarAll1,...
    stateVarAll8,stateVarAll30,stateVarAll45,...
    stateVarAll70, OutputAll,OutputFCAll,OutputAll1,...
    OutputAll8,OutputAll30,OutputAll45,OutputAll70,...
    median_stateVarAll,median_stateVarFCAll,median_stateVarAll1,...
    median_stateVarAll8,median_stateVarAll30,median_stateVarAll45,...
    median_stateVarAll70,median_OutputAll,median_OutputFCAll,...
    median_OutputAll1,median_OutputAll8,median_OutputAll30,...
    median_OutputAll45,median_OutputAll70,...
    mean_stateVarAll,mean_stateVarFCAll,mean_stateVarAll1,...
    mean_stateVarAll8,mean_stateVarAll30,mean_stateVarAll45,...
    mean_stateVarAll70, mean_OutputAll,mean_OutputFCAll,mean_OutputAll1,...
    mean_OutputAll8,mean_OutputAll30,mean_OutputAll45,mean_OutputAll70] ...
    = populationRun(numSims, oligo,rotenone,AA,CIV,FCCP,energy,time,...
    v_xpar7,v_xpar9,v_xpar10,v_xpar13,v_xpar15,v_xpar21,...
    v_otherpar4,v_otherpar6,v_xpar41,...
     t_prior,t_start,t_final,tt,stepsize,...
     options, xo_single_cell, plot_ss, printSim);
 
 
%Return required Outputs and state variables
X1 = OutputAll1(:, 10);
X2 = stateVarAll1(:,25);
X3 = stateVarAll1(:,19);
X4 = OutputAll1(:,25);

 
 

% Print again just to make sure you don't miss it!
fprintf('\n********************')
if v_xpar7 ~= 1, fprintf('\nImpairing x_DH [xpar(7)=%0.2i]\n',v_xpar7), end
if v_xpar9 ~= 1, fprintf('\nImpairing x_C1 [xpar(9)=%0.2i]\n',v_xpar9), end
if v_xpar10 ~= 1, fprintf('\nImpairing x_C3 [xpar(10)=%0.2i]\n',v_xpar10), end
if v_xpar13 ~= 1, fprintf('\nImpairing x_C4 [xpar(13)=%0.2i]\n',v_xpar13), end
if v_xpar15 ~= 1, fprintf('\nImpairing x_F1 [xpar(15)=%0.2i]\n',v_xpar15), end
if v_xpar21 ~= 1, fprintf('\nImpairing x_Hle [xpar(21)=%0.2i]\n',v_xpar21), end
if v_otherpar4 ~= 1, fprintf('\nImpairing K_ADTP_dyn [otherpar(4)=%0.2i]\n',v_otherpar4), end
if v_otherpar6 ~= 1, fprintf('\nImpairing K_ADTP_cons [otherpar(6)=%0.2i]\n',v_otherpar6), end
fprintf('********************\n')

if setRNG == 1
    fprintf('\nRandom number generator seed was reset [rng(1)] prior to simulations.\n')
end
fprintf('\nRunning expt %i\n', expt)
end