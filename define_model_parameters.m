function [xpar] = define_model_parameters
%*************************************************************************
% Define model parameters
% Updated from Huber 2011 (PMID: 21364572) and Huber 2012 (PMID: 22218564)
% Refs for parameter values in Supp Table 1-4
%*************************************************************************
%
% Metabolites
%
xpar    = zeros(31,1);
xpar(1) = 726e-6;     % NADtot    
xpar(2) = 2.70e-3;    % Cyt-c (total IMS Cyt-C, Cred+Cox, M)
xpar(3) = 1.35e-3;    % total IMS Ubiquinol, Q+QH2, M
xpar(4) = 2.6e-3;     % ADTPtot; total Adenosine phosphates in cytosol
                                            
%*************************************************************************
% FLUX 1: Parameter for input function
%*************************************************************************
xpar(5) = 0.13413e-3;    % k_Pi1      % Dehydrogenase flux input
xpar(6) = 0.677e-3;      % k_Pi2      % Dehydrogenase flux input
xpar(7) = 0.8*0.0737;    % x_DH       % Dehydrogenase activity
xpar(8) =  4.3;          % r_DH       % Input-flux: Initial disturbance of equilibrium
%*************************************************************************
% FLUX 2: Parameter for complex I
%**************************************************x_c4***********************
xpar(9) = 1.0200e+003;          % x_C1       % Complex I activity
%*************************************************************************
% FLUX 3: Parameter for complex III
%*************************************************************************
xpar(10) = 0.2241;              % x_C3      % Complex III activity
xpar(11) = 0.192e-3;            % k_Pi3     % Complex III / Pi parameter
xpar(12) = 25.31e-3;            % k_Pi4     % Complex III / Pi parameter
%*************************************************************************
% FLUX 4: Parameter for complex IV
%*************************************************************************
xpar(13) = 3.2131e-004;         % x_C4      % Complex IV activity
xpar(14) = 1.2e-4;              % k_O2      % kinetic constant for complex IV (M)
%*************************************************************************
% FLUX 5: Parameter for ATP synthase and Mg-binding to ATP/ADP
%*************************************************************************
xpar(15) = 6.8294e+003;         % x_F1      % F1Fo ATPase activity
xpar(16) = 8*24e-6;             % K_DT      % Mg/ATP binding constant (M)
xpar(17) = 347e-6;              % K_DD      % Mg/ADP binding constant (M)
xpar(18) = 1e6;                 % x_MgA     % Mg2+ binding activity
%*************************************************************************
% FLUX 6: Parameter for Adenosine Transferase
%*************************************************************************
xpar(19) = 0.0020;              % x_ANT     % ANT activity
xpar(20) = 3.50e-6;             % k_mADP    % ANT Michaelis-Menten constant
%*************************************************************************
% FLUX 7: Proton leak activity
%*************************************************************************
xpar(21) = 150.0;               % x_Hle     % Proton leak activity  
%*************************************************************************
% FLUX 8: Parameter for OM transporters
%*************************************************************************
xpar(22) = 200e1;               % x_Ht      % MOM permeability to protons (micron s^{-1})
xpar(23) = 5.99;                % gamma     % MOM area per unit volume micron^{-1}
xpar(24) = 85.0;                % x_A       % MOM permeability to nucleotides (micron s^{-1}) 
xpar(25) = 327;                 % x_Pi2     % MOM permeability to phosphate (micron s^{-1})
%*************************************************************************
% FLUX 9: Phosphate-Hydrogen-Cotransport
%*************************************************************************
xpar(26) = 10^(-6.75);          % k_dHPi    % H/Pi co-transpor Molar (form factor 1) (binding)
xpar(27) = 0.45082e-3;          % k_PiH     % H+/Pi co-transport activity (form factor 2) (Michaelis constant)
xpar(28) = 3.85e5;              % x_Pi1     % H+/Pi co-transport activity
%*************************************************************************
% FLUX 10: Potassium-Hydrogen-Antiport
%*************************************************************************
xpar(29) = 2.9802e7;            % x_KH      % K+ / H+ antiporter activity

xpar(30) = 100;                 % x_buff    % Inner Matrix hydrogen buffer capacity
xpar(31) = 6.7568e-6;           % CIM       % Inner Membrane capacitance

%*************************************************************************
% Reactive oxygen species
%*************************************************************************

xpar(32) = 0.046176;          %a_Qo  Complex III regression constant 1   
xpar(33) = 6.03612;           %m_Qo  Complex III regression constant 2 
xpar(34) = 0.0796868;         %c_Qo  Complex III regression intercept
xpar(35) = 1e-8;              %x_Qo  Complex III H2O2 scaling factor 
xpar(36)= 0.142;              %a_If  Complex I regression constant 1  
xpar(37)= 2.7;                %m_If  Complex I regression constant 2 
xpar(38) = 1.0003;            %c_If  Complex I regression intercept 
xpar(39) = 1.22;              %k_If  Boltzmann slope factor 
xpar(40) = 1.2e-7;            %x_If  Complex I H2O2 scaling factor
xpar(41) = 120.0;             %k_H2O2 First order H2O2 degradation kinetic 
xpar(42) = 0.63;              %e_If   delPsi-dependent H2O2 production   

end