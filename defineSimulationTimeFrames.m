function [t_prior, t_start,t_final,t_no_time,stepsize,tt,time]...
    = defineSimulationTimeFrames(printRequest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set simulation time frames and other time settings
%%% Set printRequest to 1 if you want to print outputs as below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_prior = -80;          % Time point prior to steady-state calcs
t_start= 0;             % Simulation start time (after ss calcs)
t_final = 80;           % Simulation time in min
t_no_time = t_final+10; % Dummy time for neglected effects
stepsize = 0.1;         % Simulation Step Size [min] to sample the solution
tt = 0:stepsize:t_final;    % Allocate same # of data points to each result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time.ATPvary:  Time at which dydt for ATP_c and ADP_c are defined 
%%% (otherwise f(23), f(24) = 0).  
%%% ATP_c, ADP_c vary when t > time.ATPvary
%%% Set to 0 to allow cytosolic ATP to vary (single-cell conditions); 
%%% Set to t_no_time to fix cytosolic ATP and ADP (e.g. isolated 
%%% mitochondria with external media buffered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.ATPvary      = -80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time.ATPcyto: Time when cytosolic ATP production and consumption starts  
%%% Set to t_no_time if e.g. you want to model isolated mitochondria, or 
%%% measure/calculate respiratory control ratio (RCR; State3/State4).
%%% Set to 0 to model single (whole) cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.ATPcyto      = -80; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cytochrome-c release (after MOMP)
%%% time.release = time of cyt-c release (depolarisation) in minutes
%%% Set to t_no_time for no depolarisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.release      = t_no_time;        

if printRequest == 1
    fprintf('ATP set to vary (time.ATPvary) at t = %0.1i min\n', time.ATPvary)
    fprintf('Cytosolic ATP prod.& cons. (time.ATPcyto) initiated at t = %0.1i min\n',...
    time.ATPcyto)
end
end