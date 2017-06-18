%% See pattern
clear all, close all, clc
destin  = 'p:\1230829-marshall_risk_assessment\08_papers\Samoa_2017_impact_typhoons\runs\RMI_enhanced\version004\_clean\'
cd(destin); load 'tcstructure.mat'; load spwstucture.mat;
vmaxs = linspace(10,70, 100);
lats  = linspace(0,0, 100);
times = linspace(0,600,100)/24;
x     = linspace(1,20, 100);
y     = linspace(1,20, 100);
averagepressure     = 1.0102e+05; 

% Create tc
addpath('p:\1230829-marshall_risk_assessment\03_data_analysis\typhoons\functions');
for it = 1:length(times)

    [vmax,dp] = marshall_typhoon_umax_dp(vmaxs(it),[]);   
    if dp.marshall > 0;
        tc.track.pc(it,1) = averagepressure/100 - real((dp.marshall));
    else
        tc.track.pc(it,1) = averagepressure/100;
    end
    
    % Real values
    tc.track.time(it,1)     = times(it);
    tc.track.x(it,1)        = x(it);
    tc.track.y(it,1)        = y(it);
    tc.track.vmax(it,1)     = vmaxs(it);


    % Determine radius
    [rmax_mean,dr35_mean, rmax_RMSE, dr35_RMSE]           = rmax_r35(vmaxs(it),y(it), 3)
    tc.track.rmax(it,1)     = rmax_mean;
    
    % Starte r35 if possible;
    r35                     = rmax_mean+dr35_mean;
    if isnan(r35); r35 = -999; end
    tc.track.r35ne(it,1) = r35 ;    tc.track.r35se(it,1) = r35;     tc.track.r35sw(it,1) = r35;    tc.track.r35nw(it,1) = r35;
    tc.track.r50ne(it,1) = -999 ;    tc.track.r50se(it,1) = -999;     tc.track.r50sw(it,1) = -999;    tc.track.r50nw(it,1) = -999;
    tc.track.r65ne(it,1) = -999 ;    tc.track.r65se(it,1) = -999;     tc.track.r65sw(it,1) = -999;    tc.track.r65nw(it,1) = -999;
    tc.track.r100ne(it,1) = -999 ;    tc.track.r100se(it,1) = -999;     tc.track.r100sw(it,1) = -999;    tc.track.r100nw(it,1) = -999;
end

% Create tc and spw
spw.reference_time  = datenum(2017,01,01);
tc.wind_speed_unit  = 'ms';
tc.radius_unit      = 'km';
tc.radius_velocity  = [34,50,64,100] * 0.514;
spw.pn = averagepressure/100;
spw.cut_off_speed = 5;
tc.factor           = 0.88; % 1 minute to 10 minutes
tc      = wes3(tc,'tcstructure',spw,['tst.spw']); 

figure; hold on; 
for ii = 1:100
    for jj = 1:36
        plot([tc.track(ii).wind_speed(jj,:)])
    end
end
