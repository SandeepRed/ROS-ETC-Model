%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Takes experimental conditions as input, simulates population using 
%%% these experimental conditions, returns outputs and
%%% the mean and median of all the outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xparAll, stateVarAll,stateVarFCAll,stateVarAll1,...
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
     options, xo_single_cell, plot_ss, printSim)

global C
 
for s = 1:numSims
  %%% Initialise model parameters for each simulation
  xpar = define_model_parameters;
  otherpar = define_other_parameters;
 
  %%% Apply impairments defined outside for loop (if any)
  xpar(7) = xpar(7) * v_xpar7;    % x_DH
  xpar(9) = xpar(9) * v_xpar9;    % x_C1
  xpar(10) = xpar(10) * v_xpar10; % x_C3 
  xpar(13) = xpar(13) * v_xpar13; % x_C4
  xpar(15) = xpar(15) * v_xpar15; % x_F1
  xpar(21) = xpar(21) * v_xpar21; % x_Hle  
  otherpar(4) = otherpar(4) * v_otherpar4;  % K_ADTP_dyn
  otherpar(6) = otherpar(6) * v_otherpar6;  % K_ADTP_cons
  xpar(41) = xpar(41) * v_xpar41; % ROS Scavenging  

         
  % Run first simulation with 'original' parameters
  if s == 1 
    % For the first simulation, save original values so can vary around
    % these original values in subsequent simulations
    xpar_orig = xpar;
    otherpar_orig = otherpar;
    oligo_orig = oligo;
    rotenone_orig = rotenone;
    AA_orig = AA;
    CIV_orig = CIV;
    FCCP_orig = FCCP;
    xo_single_cell_orig = xo_single_cell; 
    % Save matrices to map variation of parameters
    xparAll = xpar;
    otherparAll = otherpar;
    xo_single_cellAll = xo_single_cell';
    
  elseif s > 1  
    % Call function to apply parameter variation
    [xo_single_cell, xpar, otherpar, oligo, rotenone, AA, CIV, FCCP]...
    = populationGenerate(xo_single_cell_orig, xpar_orig, otherpar_orig, ...
    oligo_orig, rotenone_orig, AA_orig, CIV_orig, FCCP_orig);
    
    % Save all params (to plot variation if required)
    xparAll = [xparAll xpar]; 
    otherparAll = [otherparAll otherpar];
    xo_single_cellAll = [xo_single_cellAll xo_single_cell'];
  end
      
  % Want to catch integration errors (so these values don't get saved)
    warnId = 'MATLAB:ode15s:IntegrationTolNotMet';
    warning('error', warnId);
    
  try   % Implementing error catching
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Steady-state Calculations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  oligo.check = 0;    % To ensure oligo doesn't get reset in ss calcs
  
  % Run sub-routine for set period to get steady-state
  [t0, y_ss]= ode15s(@sub_energetic,[t_prior t_start],xo_single_cell,...
      options,xpar,otherpar, time,oligo,rotenone,AA,CIV,FCCP,energy);
  if printSim == 1
      fprintf('Steady-state calculated. (s = %i) \n', s)
  end

  if plot_ss == 1
    figure('Name','SS calculations (state vars)')
    for i=1:size(y_ss,2)
        subplot(4,6,i),plot(t0,y_ss(:,i));
        title(statevar_legend(i))
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% ODE Calculations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  oligo.check     = 1;    % Oligo trigger (stores J_F1 value prior to Oligo)
  [t, y] = ode15s(@sub_energetic,[t_start t_final],y_ss(end,:),options,...
      xpar,otherpar, time,oligo,rotenone,AA,CIV,FCCP,energy);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Calculate outputs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Output = zeros(length(t),size(C,2));
  for i = 1:length(t)
    t_out = t(i); 
    y_out = y(i,:)'; 
    [f] = sub_energetic(t_out,y_out,xpar,otherpar,time,oligo,...
        rotenone,AA,CIV,FCCP,energy);
    Output(i,:) = C;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Calculate fold changes of state variables and outputs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Convert all solutions to set number of datapoints
  % Careful here because this samples the solution - may miss very rapid
  % changes
  for i = 1:size(y,2)    
    stateVarSpline(:,i) = spline(t, y(:,i), tt);    
  end
  for i = 1:size((Output),2)    
    OutputSpline(:,i) = spline(t, Output(:,i), tt);    
  end

  % Calculate fold change (divide by baseline value)
  % Baseline is second value of solution (to ignore initial fluctuations)
  for i = 1:size(y,2)
    stateVarFC(:,i) = (stateVarSpline(:,i)./stateVarSpline(2,i)); 
  end
  for i = 1:size(Output,2)
    OutputFC(:,i) = (OutputSpline(:,i)./OutputSpline(2,i)); 
  end
  
  % Save solutions into 3d matrices 
  stateVarAll(:,:,s) = stateVarSpline;
  stateVarFCAll(:,:,s) = stateVarFC;
  OutputAll(:,:,s) = OutputSpline;
  OutputFCAll(:,:,s) = OutputFC;
       
  % Take values at t=1, 30, 45, 70 
  tSample1 = 1/stepsize;
  tSampleEnergy = 8/stepsize;
  tSample2 = 29/stepsize; %sometimes fluctuations happen at t =30
  tSample3 = 44/stepsize;
  tSample4 = 70/stepsize;
  stateVarAll1(s,:) = stateVarSpline(tSample1,:);
  OutputAll1(s,:) = OutputSpline(tSample1,:);
  stateVarAll8(s,:) = stateVarSpline(tSampleEnergy,:);
  OutputAll8(s,:) = OutputSpline(tSampleEnergy,:);
  stateVarAll30(s,:) = stateVarSpline(tSample2,:);
  OutputAll30(s,:) = OutputSpline(tSample2,:);
  stateVarAll45(s,:) = stateVarSpline(tSample3,:);
  OutputAll45(s,:) = OutputSpline(tSample3,:);
  stateVarAll70(s,:) = stateVarSpline(tSample4,:);
  OutputAll70(s,:) = OutputSpline(tSample4,:);
  
  catch ME    % If the integration tolerance warning is generated, don't save values!
      fprintf('Simulation %i: Integration error (result saved as NaN).\n',s)
      % Set to 0 so this simulation will still be counted
      stateVarAll(:,:,s) = zeros(t_final/stepsize+1,size(xo_single_cell_orig,2),1);
      stateVarFCAll(:,:,s) = zeros(t_final/stepsize+1,size(xo_single_cell_orig,2),1);
      stateVarAll1(s,:) = zeros(1,size(xo_single_cell_orig,2));
      stateVarAll8(s,:) = zeros(1,size(xo_single_cell_orig,2));
      stateVarAll30(s,:) = zeros(1,size(xo_single_cell_orig,2));
      stateVarAll45(s,:) = zeros(1,size(xo_single_cell_orig,2));
      stateVarAll70(s,:) = zeros(1,size(xo_single_cell_orig,2));

      OutputAll(:,:,s) = zeros(t_final/stepsize+1,28,1);
      OutputFCAll(:,:,s) = zeros(t_final/stepsize+1,28,1);
      OutputAll1(s,:) = zeros(1,28);
      OutputAll8(s,:) = zeros(1,28);
      OutputAll30(s,:) = zeros(1,28);
      OutputAll45(s,:) = zeros(1,28);
      OutputAll70(s,:) = zeros(1,28);
  end
    
end     % End of population simulations 
% Dimensions of stateVarAll (and others)
% 1: time (~800 rows)
% 2: each column is a state variable (24 columns)
% 3: individual simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Any simulations that did not integrate - replace with NaN
for i = 1:size(stateVarAll,3)
  if stateVarAll(1,1,i) == 0
    stateVarAll(:,:,i) = NaN;
    stateVarFCAll(:,:,i) = NaN;
    stateVarAll1(i,:) = NaN;
    stateVarAll8(i,:) = NaN;
    stateVarAll30(i,:) = NaN;
    stateVarAll45(i,:) = NaN;
    stateVarAll70(i,:) = NaN;
    OutputAll(:,:,i) = NaN;
    OutputFCAll(:,:,i) = NaN;
    OutputAll1(i,:) = NaN;
    OutputAll8(i,:) = NaN;
    OutputAll30(i,:) = NaN;
    OutputAll45(i,:) = NaN;
    OutputAll70(i,:) = NaN;
  end
end

% Calculate mean and median of all simulations
% Calculate along 3rd dimension of 3d matrices
mean_stateVarAll = nanmean(stateVarAll,3);
mean_stateVarFCAll = nanmean(stateVarFCAll,3);
mean_stateVarAll1 = nanmean(stateVarAll1);
mean_stateVarAll8 = nanmean(stateVarAll8);
mean_stateVarAll30 = nanmean(stateVarAll30);
mean_stateVarAll45 = nanmean(stateVarAll45);
mean_stateVarAll70 = nanmean(stateVarAll70);
mean_OutputAll = nanmean(OutputAll,3);
mean_OutputFCAll = nanmean(OutputFCAll,3);
mean_OutputAll1 = nanmean(OutputAll1);
mean_OutputAll8 = nanmean(OutputAll8);
mean_OutputAll30 = nanmean(OutputAll30);
mean_OutputAll45 = nanmean(OutputAll45);
mean_OutputAll70 = nanmean(OutputAll70);

median_stateVarAll = nanmedian(stateVarAll,3);
median_stateVarFCAll = nanmedian(stateVarFCAll,3);
median_stateVarAll1 = nanmedian(stateVarAll1);
median_stateVarAll8 = nanmedian(stateVarAll8);
median_stateVarAll30 = nanmedian(stateVarAll30);
median_stateVarAll45 = nanmedian(stateVarAll45);
median_stateVarAll70 = nanmedian(stateVarAll70);
median_OutputAll = nanmedian(OutputAll,3);
median_OutputFCAll = nanmedian(OutputFCAll,3);
median_OutputAll1 = nanmedian(OutputAll1);
median_OutputAll8 = nanmedian(OutputAll8);
median_OutputAll30 = nanmedian(OutputAll30);
median_OutputAll45 = nanmedian(OutputAll45);
median_OutputAll70 = nanmedian(OutputAll70);

end