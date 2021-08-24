function [otherpar] = define_other_parameters
% Define parameters not defined in xpar or sub_energetic (e.g.
% parameters you might want to vary!)
% Refs for parameters in Supp Table 1-4

otherpar    = zeros(3,1);
otherpar(1) = 120e-3;       % K_ei; cytosolic and IM potassium-concentration
otherpar(2) = 20e-3;        % Mg_tot; cytosolic and IM Mg-concentration
otherpar(3) = 20e-3;        % Pi_e; cytosolic and IM phosphate-concentration
otherpar(4) = 0.9*3.8;      % K_ADTP_dyn; Cytosolic ATP production
otherpar(5) = 0.8*0.63;     % x_ATPK; 
otherpar(6) = 1;            % K_ADTP_cons; Cytosolic ATP consumption (redundant parameter)

end