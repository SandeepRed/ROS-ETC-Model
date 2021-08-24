function [xo_single_cell, xpar, otherpar, oligo, rotenone, AA, CIV, FCCP]...
    = populationGenerate(xo_single_cell_orig, xpar_orig, otherpar_orig,...
    oligo_orig, rotenone_orig, AA_orig, CIV_orig, FCCP_orig)

% Generate population of 'cells' by varying parameters, initial
% concentrations and drug applications

%%% Define range of variations (0.2 => +/- 20%)
varyP = 0.2;    % Parameters 
varyI = 0.2;    % Initial concs
varyD = 3;      % Drug effect (3 => can vary 1/3-3 times defined value)
varyFCCP = 0.2;    % FCCP factor treated separately - vary between +/-20%
varyDt = 0;     % Time of drug effect
    
numStateVars = size(xo_single_cell_orig,2);  % Number of state variables
numParams = size(xpar_orig,1);   % Number of parameters
numOtherParams = size(otherpar_orig,1);
    
%%% Vary model parameters and initial concentrations 
%%% Note: Only varying parameters defined in xpar 

% Apply normally distributed random value between +/- varyP to each
% parameter (normal distribution has mean = 0 and sd = 0.5*varyP)
randParams = (0.5*varyP).*randn(numParams,1) +0;
randOtherParams = (0.5*varyP).*randn(numOtherParams,1) +0;

xpar = xpar_orig + xpar_orig.*randParams;
otherpar = otherpar_orig + otherpar_orig.*randOtherParams;
    
randInitC = (0.5*varyI).*randn(numStateVars,1) +0;
randInitC = randInitC';

xo_single_cell = xo_single_cell_orig + xo_single_cell_orig.*randInitC;
  
%Vary extent of drug effect (5 drugs)
% Normally distribute variations between 1/3 - 3 times original drug value
% First, obtain normal distrib with mean 0 and sd = 1/4 of varyD (=3/4)
% =) 2sd ~1.5

randDrugs = (0.25*varyD).*randn(4,1) +0; %(-1.58->+1.58)
randFCCP = (0.5*varyFCCP).*randn(1,1) +0;

% Then 2^ to centre distribution around 1, ranging from 1/3 to 3 (all values positive)
% (2^-1.5 ~1/3; 2^1.5 ~ 3)
randDrugs = 2.^randDrugs; %(1/3->3) 

oligo.percent       = oligo_orig.percent*randDrugs(1);
rotenone.percent    = rotenone_orig.percent*randDrugs(2);
AA.percent          = AA_orig.percent*randDrugs(3);      
CIV.percent         = CIV_orig.percent*randDrugs(4);
% FCCP factor is bigger number
FCCP.factor         = FCCP_orig.factor + FCCP_orig.factor*randFCCP;
    
% Vary drug time (to allow for diffusion effects etc.)
% VaryDt currently set to 0, so no variation here
randDTime = ((varyDt + varyDt).*rand(5,1)) - varyDt;
oligo.t       = abs(oligo_orig.t + oligo_orig.t*randDTime(1));
rotenone.t    = abs(rotenone_orig.t + rotenone_orig.t*randDTime(2));
AA.t          = abs(AA_orig.t + AA_orig.t*randDTime(3));      
CIV.t         = abs(CIV_orig.t + CIV_orig.t*randDTime(4));
FCCP.t         = abs(FCCP_orig.t + FCCP_orig.t*randDTime(5));

end