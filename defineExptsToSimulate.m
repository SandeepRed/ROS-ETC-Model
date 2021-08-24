function [rotenone,AA,oligo,CIV,FCCP,energy] = ...
    defineExptsToSimulate(expt,t_no_time, printRequest)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define experiments that should be simulated, and set conditions for these
% experiments (time and extent of drug additions/complex inhibition)
% 1: Rotenone, Oligo
% 2: AA, Oligo
% 3: Seahorse (Oligo, FCCP, Rot+AA)
% 4: Increased energy demand (IncEnDem, Oligo, FCCP, Rot+AA) 
% 5: FCCP, Rotenone 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define default values to return if values not defined explicitly below
% Onset of drug addition/complex inhibition (minutes)
[rotenone, AA, oligo, CIV, FCCP, energy] = defineDefaultDrugCond(t_no_time);

% Based on input s, define experimental conditions to simulate
switch expt
    % Rotenone (Rot, Oligo)
    case 1  
        rotenone.t = 10; 
        oligo.t = 400;
        FCCP.t = 700;
    
    % Antimycin A (AA, Oligo)
    case 2  
        AA.t = 10; 
        oligo.t = 400;
    
    % Seahorse (Oligo, FCCP, Rot+AA)
    case 3  
        rotenone.t = 54; 
        AA.t = 54; 
        oligo.t = 18; 
        FCCP.t = 36;
        % Set to ~0 to align with expts where measurements after Rot/AA 
        % drugs are considered non-mitochondrial respiration 
        % Setting = 0 returns NaN and Inf values
        rotenone.percent = 1e-20;   
        AA.percent = 1e-20;
    
    % Seahorse with increased energy demand (IncEnDem, Oligo, FCCP, Rot+AA)
    % Not used in Theurey et al 2018
    case 4 
        rotenone.t = 54; 
        AA.t = 54; 
        oligo.t = 15;   
        FCCP.t = 40; 
        energy.t = 3;   
        rotenone.percent = 1e-20;   
        AA.percent = 1e-20; 
    
    % FCCP, Rotenone
    case 5 
        rotenone.t = 40; 
        FCCP.t = 10;   
        AA.t = t_no_time; 
        oligo.t = t_no_time; 
end

if printRequest == 1
    fprintf('\n*********')
    fprintf('\nRot at %i; AA at %i; Oligo at %i; FCCP at %i; energy at %i.\n',...
            rotenone.t,AA.t,oligo.t,FCCP.t,energy.t)
    if energy.t < t_no_time
    fprintf('(Simulating increased energy demand)\n')
    end
end
end