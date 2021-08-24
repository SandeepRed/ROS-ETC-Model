function [rotenone, AA, oligo, CIV, FCCP, energy] = defineDefaultDrugCond(t_no_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define default drug conditions (inhibition of complexes) if values not 
% defined explicitly elsewhere
%
% Define time and extent of drug addition
%%% Oligomyin (oligo): ATP synthase inhibition
%%% Rotenone: Complex I inhibition
%%% Antimycin A (AA): Complex III inhibition
%%% CIV: Complex IV inhibition
%%% FCCP: Proton leak increase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Onset of drug addition/complex inhibition (minutes)
% Set to t_no_time for no drug addition
    oligo.t         = t_no_time; 
    FCCP.t          = t_no_time;
    rotenone.t      = t_no_time;
    AA.t            = t_no_time;    
    CIV.t           = t_no_time;
    energy.t        = t_no_time;    % time to initiate increased energy demand

% Extent of inhibition (standard values) 
% Determined from drug response curves from literature, 
% Recorded in "ICs and Ks from lit.xls" 

    oligo.percent       = 0.0004/100; 
    rotenone.percent    = 0.14/100;       
    AA.percent          = 0.8/100;       
    CIV.percent         = 100/100;
    % Increase proton leak to max. Any higher values gives ~ same result
    % for JC_4 (i.e. simulating max respiration)
    % (see "ICs and Ks from lit.xls")
    FCCP.factor         = 13;    
    energy.factor       = 1; % If using this, need to re-calibrate
    
end