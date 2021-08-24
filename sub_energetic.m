function [f] = sub_energetic(t,x,xpar,otherpar, time,oligo,rotenone,AA,CIV,FCCP,energy)
%==================================================================
%==================================================================
% The bioenergetic model is derived from the great work of
% Daniel A. Beard, Department of Physiology, Medical College of Wisonsin,
% Milwaukee, WI.
%
% Beard DA. (2005) A Biophysical Model of the Mitochenondrial Respiratory
% System and Oxidative Phosphorylation. PLoS Comp Bio. 1 (4) e36
%
% Including the ATP/ADP transfer modelling from:
% Korzeniewski B. (1998), Biophys Chem 75(1), 73-80
%
% Huber 2011 incorporated cytosolic ATP production and consumption
%==========================================================================
%==========================================================================

global C pp_casp3
               
%*************************************************************************
% Concentration state variables: (rename them for clarity)
%*************************************************************************
H_x    = x(1);                      % Matrix hydrogen concentration
K_x    = x(2);                      % Matrix potassium concentration
Mg_x   = x(3);                      % Matrix magnesium concentration
NADH_x = x(4);                      % Reduced NADH in matrix
QH2    = x(5);                      % Reduced ubiquinol in matrix
Cred   = x(6);                      % Reduced cyt-c in IMS
O2     = x(7);                      % O2 unused as state variable, parameter only
ATP_x  = x(8);                      % Matrix free ATP (set to external as initial)
ADP_x  = x(9);                      % Matrix free ADP (set to external as initial)
ATP_mx = x(10);                     % Matrix ATP bound to magnesium
ADP_mx = x(11);                     % Matrix ADP bound to magnesium
Pi_x   = x(12);                     % Matrix inorganic phosphate
ATP_i  = x(13);                     % IMS free ATP
ADP_i  = x(14);                     % IMS free ADP
AMP_i  = x(15);                     % IMS free AMP
ATP_mi = x(16);                     % IMS ATP bound to magnesium
ADP_mi = x(17);                     % IMS ADP bound to magnesium
Pi_i   = x(18);                     % IMS inorganic phosphate
dPsi   = x(19);                     % Mitochondrial membrane potential
Ctot   = x(20);                     % total cyt-c in IMS
Qtot   = x(21);                     % total ubiquinol
H_i    = x(22);                     % IMS hydrogen concentration
ATP_e  = x(23);                     % cytosolic ATP
ADP_e  = x(24);                     % cytosolic ADP
AMP_e  = 0;                         % cytosolic AMP concentration (Molar), Nicholls and others
H2O2 = x(25);                       % cytosolic H2O2 concentration (Molar)


%*************************************************************************
%   Model parameters (fixed parameters)
%   References in Supplementary Table 1-4
%*************************************************************************
%   Thermodynamic parameters
F       = 0.096484;                  % kJ mol^{-1} mV^{-1}
R       = 8314e-6;                   % Universal gas constant (kJ * mol^{-1} * K^{-1})
T       = (273.15 + 37);             % Temperature (K), 37 degree
RT      = R*T;                       % kJ  mol^{-1}

%   Gibbs free energy parameters
dG_C1o  = -69.37;                    % kJ mol^{-1}
dG_C3o  = -32.53;                    % kJ mol^{-1}
dG_C4o  = -122.94;                   % kJ mol^{-1}
dG_F1o  = 36.03;                     % kJ mol^{-1}
n_A     = 3.0;                       % numbers of proteins used by ATP synthase

%   Cytosolic Ion concentration and pH
pH_e    = 7.4;             % External pH (cytosol)
H_e     = 10^(-pH_e);      % cytosolic hydrogen concentration (Molar)
K_ei    = otherpar(1);     % cytosolic and IM potassium-concentration 
Mg_tot  = otherpar(2);     % cytosolic and IM magnesium-concentration 
Pi_e    = otherpar(3);     % cytosolic and IM phosphate-concentration
                           
c_inc   = 1;               % Cyt-c after release set to 0.1%
V_mito  = 0.06;            % 6% mitochondrial fraction Ward 2007
W_c     = 1/V_mito;        % Volume fraction cytosol/mitochondria
        
%   Matabolites
NADtot      = xpar(1);     % See define_model_parameters for xpar
Ctot0       = xpar(2);           
Qtot0       = xpar(3);          

% Parameters for input function
k_Pi1  =  xpar(5);                 
k_Pi2  =  xpar(6);                
x_DH   =  xpar(7);                 
r_DH   =  xpar(8);                

% Parameters for complex I
x_C1   =  xpar(9);                 

% Parameters for complex III
x_C3   =  xpar(10);                
k_Pi3  =  xpar(11);              
k_Pi4  =  xpar(12);               

% Parameters for complex IV
x_C4   =  xpar(13);              
k_O2    = xpar(14);              

% Parameters for ATP synthase and Mg-binding to ATP?ADP
x_F1   =  xpar(15);              
K_DT    = xpar(16);             
K_DD    = xpar(17);             
x_MgA  =  xpar(18);              

% Parameters for Adenosine Transferase
x_ANT  =  xpar(19);               
k_mADP =  xpar(20);              

% Proton leak activity
x_Hle  =  xpar(21);              

% Parameters for OM transporters
x_Ht    = xpar(22);               
gamma   = xpar(23);               
x_A     = xpar(24);              
x_Pi2   = xpar(25);               

% Phosphate-Hydrogen-Cotransport
k_dHPi  = xpar(26);              
k_PiH   = xpar(27);                
x_Pi1   = xpar(28);             

% Potassium-Hydrogen-Antiport
x_KH   =  xpar(29);                

% Matrix buffer and membrane capacitance
x_buff =  xpar(30);          
CIM     = xpar(31);          


% ROS metabolism parameters                
a_Qo   =   xpar(32);
m_Qo   =   xpar(33);            
c_Qo   =   xpar(34);          
x_Qo   =   xpar(35);           
a_If   =   xpar(36);  
m_If   =   xpar(37);
c_If   =   xpar(38);
k_If   =   xpar(39);   
x_If   =   xpar(40);              
k_H2O2 =   xpar(41); 
e_If   =   xpar(42);       

% Release and compartment parameters
W_m     = 0.143/0.20;               % mitochondrial water space (ml water / ml mito)
W_x     = 0.9*W_m;                  % Matrix water space (ml water / ml mito)
W_i     = 0.1*W_m;                  % IM water space (ml water / ml mito)
t12cyto = 90;                       % cyt-c release half time Waterhouse 2000 (seconds)
Cfin    = Ctot0*W_i/W_c*c_inc;      % final cytochrome c after release in IMS given by re-equilibration of the
                                    % IMS water space with the rest of the cell (0.1%)
                                    % c_inc is used to modify OXPHOS-contributing cyt-c after release

                                    
                                    
                             
                                    
                                    
                                    
                                    
%*************************************************************************
% potassium uniporter and adenylate kinase neglected
%*************************************************************************
x_K    = 0;                         % Passive potassium transporter activity
x_AK   = 0e6;                       % AK activity
K_AK    = 0;                        % Adenelyte Kinase switched off

%*************************************************************************
% Balancing moeities
%*************************************************************************
NAD_x  = NADtot - NADH_x;
Q      = Qtot   - QH2;
Cox    = Ctot   - Cred;
null = 0;

%*************************************************************************
% Other concentrations computed from the state variables
%*************************************************************************
ATP_fx  = ATP_x - ATP_mx;
ADP_fx  = ADP_x - ADP_mx;
ATP_fi  = ATP_i - ATP_mi;
ADP_fi  = ADP_i - ADP_mi;

%*************************************************************************
% Preventing numerical instabilities
% Sets to 0 if value < 0
%*************************************************************************
        Q       = Q*(Q>0);
        QH2     = QH2*(QH2>0);
        Qtot    = Qtot*(Qtot>0);
        Cox     = Cox*(Cox>0);
        Cred    = Cred*(Cred>0);
        Ctot    = Ctot*(Ctot>0);
        ATP_fx  = ATP_fx*(ATP_fx>0);
        ATP_fi  = ATP_fi*(ATP_fi>0);
        ADP_fx 	= ADP_fx*(ADP_fx>0);
        ADP_fi  = ADP_fi*(ADP_fi>0);

%*************************************************************************
% ADP/Mg/K binding in E space
%*************************************************************************
        ADP_me  = ( (K_DD+ADP_e+Mg_tot) - sqrt((K_DD+ADP_e+Mg_tot)^2-4*(Mg_tot*ADP_e)) )/2;
        ADP_fe  = ADP_e - ADP_me;
        Mg_e    = Mg_tot - ADP_me;
        Mg_i    = Mg_e;                         % Mg,K in IM space same as external/cytosolic
        K_i     = K_ei; 

%*************************************************************************
% Simulating drug additions by reducing kinetic constants (i.e. 
% activity) of relevant complexes
% Modelling Complex I block (e.g. Rotenone)
% Modelling Complex III block (e.g. Antimycin A)
% Modelling Complex IV block (e.g. CN)
% Modelling F1Fo ATP synthase block (e.g. oligomycin)
% Modelling increased proton leak (e.g. FCCP)
%*************************************************************************
        if t > rotenone.t,  x_C1 = rotenone.percent * x_C1;  end
        if t > AA.t,  x_C3 = AA.percent * x_C3;  end
        if t > CIV.t,  x_C4 = CIV.percent * x_C4;  end
        if t > oligo.t, x_F1 = oligo.percent * x_F1; end   
        if t > FCCP.t,  x_Hle = FCCP.factor * x_Hle;  end
                
%*************************************************************************
% Calculating Membrane proton motive force and respiration fluxes
%*************************************************************************
        dG_H    = F*dPsi + 1*RT*log(H_i/H_x);   % Protomotive force
        dG_C1op = dG_C1o - 1*RT*log(H_x/1e-7);
        dG_C3op = dG_C3o + 2*RT*log(H_x/1e-7);
        dG_C4op = dG_C4o - 2*RT*log(H_x/1e-7);
        dG_F1op = dG_F1o - 1*RT*log(H_x/1e-7);
        
        J_DH    = x_DH*(r_DH*NAD_x-NADH_x)*((1+Pi_x/k_Pi1)/(1+Pi_x/k_Pi2));
        J_C1    = x_C1*(exp(-(dG_C1op+4*dG_H)/RT)*NADH_x*Q - NAD_x*QH2);
        J_C3    = x_C3*((1+Pi_x/k_Pi3)/(1+Pi_x/k_Pi4))*...
            (exp(-(dG_C3op+4*dG_H-2*F*dPsi)/(2*RT))*Cox*QH2^0.5 - Cred*Q^0.5);
        J_C4    = x_C4*(O2/(O2+k_O2))*(Cred/Ctot)*...
            (exp(-(dG_C4op+2*dG_H)/(2*RT))*Cred*(O2^0.25) - Cox*exp(F*dPsi/RT));
        J_F1    = x_F1*(exp(-(dG_F1op-n_A*dG_H)/RT)*(K_DD/K_DT)*ADP_mx*Pi_x - ATP_mx);
        
%*************************************************************************
% Modelling ATP transferase Korzeniewski 1998
%*************************************************************************
        Psi_x = -0.65*dPsi;
        Psi_i = +0.35*dPsi;
        if (ADP_fi > null) || (ATP_fi > null)
            J_ANT  = x_ANT*( ADP_fi/(ADP_fi+ATP_fi*exp(-F*Psi_i/RT))...
                - ADP_fx/(ADP_fx+ATP_fx*exp(-F*Psi_x/RT)) )*(ADP_fi/(ADP_fi+k_mADP));
        else
            J_ANT  = 0;
        end
                
%*************************************************************************
% Calculating ionic fluxes
%*************************************************************************

        H2PIi       = Pi_i*H_i/(H_i+k_dHPi); 
        H2PIx       = Pi_x*H_x/(H_x+k_dHPi);
        J_Pi1       = x_Pi1*(H_x*H2PIi - H_i*H2PIx)/(H2PIi+k_PiH);
        J_Hle       = x_Hle*dPsi*(H_i*exp(F*dPsi/RT)-H_x)/(exp(F*dPsi/RT)-1);
        J_KH        = x_KH*( K_i*H_x - K_x*H_i );
        J_K         = x_K*dPsi*(K_i*exp(F*dPsi/RT)-K_x)/(exp(F*dPsi/RT)-1);
        J_AKi       = x_AK*( K_AK*ADP_i*ADP_i - AMP_i*ATP_i );
        J_AMP       = gamma*x_A*(AMP_e-AMP_i);
        J_ADP       = gamma*x_A*(ADP_e-ADP_i);
        J_ATP       = gamma*x_A*(ATP_e-ATP_i);
        J_Pi2       = gamma*x_Pi2*(Pi_e-Pi_i);
        J_Ht        = gamma*x_Ht*(H_e-H_i);

        J_MgATPx    = x_MgA*(ATP_fx*Mg_x-K_DT*ATP_mx);
        J_MgADPx    = x_MgA*(ADP_fx*Mg_x-K_DD*ADP_mx);
        J_MgATPi    = x_MgA*(ATP_fi*Mg_i-K_DT*ATP_mi);
        J_MgADPi    = x_MgA*(ADP_fi*Mg_i-K_DD*ADP_mi);

%*************************************************************************
%   Cytosolic Energy balance for single cell model (see Huber 2011)
%*************************************************************************
        K_ADTP_dyn = otherpar(4);
        % Increase cytosolic ATP consumption to simulate increased ATP
        % demand when t > energy.t
        if t > energy.t
            K_ADTP_cons = energy.factor*otherpar(6);
        else, K_ADTP_cons = otherpar(6);
        end
        
        if  t > time.ATPcyto    % Cytosolic ATP production and consumption
            % x_ATPK varied to establish desired ss equilibrium between
            % ATP:ADP (See x_ATPK and K_ADTP_dyn calcs)
            x_ATPK   = otherpar(5);    
        else
            x_ATPK   = 0;
        end
                    
        J_ATPK   = x_ATPK * (K_ADTP_cons*ATP_e- K_ADTP_dyn*ADP_e);

        
        
        

%*************************************************************************
%   Reactive oxygen species
%*************************************************************************          
        
         V_H2O2_Qo = ( a_Qo * exp( m_Qo * (QH2 /Qtot)) - c_Qo)*(x_Qo) ;         
                     
         V_H2O2_If = (a_If* exp(m_If*(NADH_x/NADtot)) + c_If)*((1-e_If) + e_If*( 1 - (  1 /  (1 + exp((dPsi - 169)/k_If)))) )* (x_If);   
         
         
%*************************************************************************
%   Calculating derivatives for next time-step
%*************************************************************************
        f(1)  = x_buff*H_x*( +1*J_DH - (4+1)*J_C1 - (4-2)*J_C3 - (2+2)*J_C4...
            + (n_A-1)*J_F1 + 2*J_Pi1 + J_Hle - J_KH )/W_x; % H_x
        f(2)  = (J_KH + J_K)/W_x; % K_x
        f(3)  = (-J_MgATPx - J_MgADPx)/W_x; % Mg_x
        f(4)  = (+J_DH - J_C1)/W_x; % NADH
        f(5)  = (+J_C1 - J_C3)/W_x; % QH2
        f(6)  = (+2*J_C3 - 2*J_C4)/W_i; % Cred
        f(7)  = 0; % O2 stays constant
        f(8)  = (+J_F1 - J_ANT)/W_x; % ATP_x
        f(9)  = (-J_F1 + J_ANT)/W_x; % ADP_x
        f(10) = (J_MgATPx)/W_x; % ATP_mx
        f(11) = (J_MgADPx)/W_x; % ADP_mx
        f(12) = (-J_F1 + J_Pi1 )/W_x;  % Pi_x 
        f(13) = (+J_ATP + J_ANT + J_AKi )/W_i; % ATP_i
        f(14) = (+J_ADP - J_ANT - 2*J_AKi )/W_i; % ADP_i
        f(15) = (+J_AMP + J_AKi)/W_i; % AMP_i
        f(16) = (J_MgATPi)/W_i; % ATP_mi
        f(17) = (J_MgADPi)/W_i; % ADP_mi
        f(18) = (-J_Pi1 + J_Pi2 )/W_i; % Pi_i
        f(19) = ( 4*J_C1 + 2*J_C3 + 4*J_C4 - n_A*J_F1 - J_ANT - J_Hle  - J_K )/CIM; %MOD
        f(20) = 0; % Cytochrome-c loss with half time t12cyto (used in Huber 2011)
        f(21) = 0; % Cleavage of complex I/II total (used in Huber 2011)
        f(22) = (- J_DH +(4+1)*J_C1 + (4-2)*J_C3 + (2+2)*J_C4 - (n_A-1)*J_F1 ...
            - 2*J_Pi1 - J_Hle + J_KH + J_Ht)/W_i;  % H_i

        if  t>time.ATPvary
            f(23) = (-J_ATP - J_ATPK)/W_c;      % ATP_c; ATP levels established by mitochondria & glycolysis
            f(24) = (-J_ADP + J_ATPK)/W_c;      % ADP_c
        else
            f(23) = 0;                          % Fixed external/cytosolic ATP for e.g. isolated mitochondria
            f(24) = 0;
        end
        
            f(25) =   (V_H2O2_If + V_H2O2_Qo- k_H2O2*H2O2)/W_c;


%*************************************************************************
% Scale Time to minutes (parameters in seconds, but output in minutes)
%*************************************************************************
        f = 60*f';                        

%*************************************************************************
% Output variables
%*************************************************************************
C(1)    = J_DH;           % Dehydrogenase flux (input function)
C(2)    = J_C1;           % Respiration/electron flux through Complex I/II
C(3)    = J_C3;           % Respiration/electron flux through Complex III 
C(4)    = J_C4;           % Respiration/electron flux through Complex IV
C(5)    = J_F1;           % ATP (FoF1) synthase flux
C(6)    = J_ANT;          % ATP transfer matrix IMS (ATP transferase inner membrane)
C(7)    = J_ATP;          % ATP transfer IMS cytosol
C(8)    = J_ADP;          % ADP transfer IMS cytosol
C(9)    = H_x;            % Matrix Hydrogen
C(10)   = V_H2O2_If;           % Total available cyt-c (IMS)
C(11)   = QH2/Qtot;           % Total available Complex1 (Total ubiquinol)
C(12) = 0;                % free active caspase3 in nM (used in Huber 2011)
C(13) = dPsi;                   % Mitochondrial Membrane potential
C(14) = dG_H/F;                 % Membrame Protomotive Force
C(15) = J_Hle;                  % Proton leak flux (Hydrogen leaks)
C(16) = ATP_e/(ADP_e+ATP_e);    % Cytosolic (external Medium for isolated mitochondria) ATP ratio
C(17) = -log10(H_x) - pH_e;     % pH difference between matrix and cytosolic pH
        
toth    = (4+1)*J_C1 +(4-2)*J_C3 +(2+2)*J_C4 + J_KH ;           % total proton extrusion

C(18)   = 100*((4+1)*J_C1 + (4-2)*J_C3 + (2+2)*J_C4)/toth;      % proton extrusion by respiration
C(19)   = 100*(J_KH)/toth;              % proton generation by potassium hydrogen
C(20)   = 100*(J_DH)/toth;              % proton consumption by input flux
C(21)   = 100*(n_A-1)*J_F1/toth;        % proton consumption by ATP synthase
C(22)   = 100*2*J_Pi1/toth;             % proton consumption by phosphate-hydrogen
C(23)   = 100*J_Hle/toth;               % proton consumption by proton leaks
C(24)   = ATP_e;                        % cytosolic ATP
C(25)   = V_H2O2_Qo;     %NADH_x/NADtot;                    % Matrix QH2 Percent
C(26)   = -log10(H_x);                  % Matrix pH
if t > time.ATPcyto  % so long as cyto ATP processes are switched on
    C(27)   = x_ATPK * K_ADTP_dyn * ADP_e;  % Cytosolic ATP production
else C(27) = 0;
end
if t > time.ATPvary     % If cytosolic ATP is 'allowed' to vary
    C(28) = x_ATPK * K_ADTP_cons * ATP_e;   % Cytosolic ATP consumption
else C(28) = 0;
end

end