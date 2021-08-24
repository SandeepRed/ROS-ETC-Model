function [Init]=initial(ADP_e, ATP_e, Ctot0, Qtot0)
%   Set initial concentrations for  state variables
%   Values may change following calculation of ss
%   Values set to non-zero values here (unlike Huber, Beard) to speed up 
%   ss calculations

Init(1)  = 10^(-8);   % H_x, Matrix hydrogen concentration pH=8 
Init(2)  = 0.050;     % K_x, Matrix potassium concentration [M]  
Init(3)  = 0.001;     % Mg_x, Matrix magnesium concentration [M] 
Init(4)  = 397e-6;    % NADH_x, Reduced NADH in matrix [M]        
Init(5)  = 0.0008;    % QH2, Reduced ubiquinol in matrix [M]    
Init(6)  = 1.0e-3;    % Cred, Reduced cyt-c in IMS [M]           
Init(7)  = 2.6e-5;    % O2, unused as state variable, parameter only
Init(8)  = ATP_e;     % ATP_x, Matrix free ATP (set to external as initial)  
Init(9)  = ADP_e;     % ADP_x, Matrix free ADP (set to external as initial)   
Init(10) = 778e-6;    % ATP_mx, Matrix ATP bound to magnesium          
Init(11) = 200e-6;    % ADP_mx, Matrix ADP bound to magnesium     
Init(12) = 10e-3;     % Pi_x, Matrix inorganic phosphate [M]    
Init(13) = 4.6e-3;    % ATP_i, IMS free ATP                               
Init(14) = 230e-6;    % ADP_i, IMS free ADP                           
Init(15) = 0;         % AMP_i, IMS free AMP                              
Init(16) = 4.6e-3;    % IMS ATP bound to magnesium
Init(17) = 226e-6;    % IMS ADP bound to magnesium
Init(18) = 10e-3;     % Pi_i, IMS inorganic phosphate [M]            
Init(19) = 150;       % dPsi, Mitochondrial membrane potential [mV]    
Init(20) = Ctot0;     % Ctot, total cyt-c in IMS                        
Init(21) = Qtot0;     % Qtot, total ubiquinol                            
Init(22) = 10^(-8.0); % H_i, IMS hydrogen concentration pH=8          
Init(23) = ATP_e;     % ATP_c, cytosolic ATP 
Init(24) = ADP_e;     % ADP_c, cytosolic ADP 
Init(25)= 1e-8;  

end