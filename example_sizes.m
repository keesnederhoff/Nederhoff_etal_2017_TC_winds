%% Hurricane examples
clear all, close all, clc

%% Load data
% Big data set
load 'd:\data\hurricanes\hurricanes.mat'

% convert to m/s and km
kts2ms = 0.514;
nm2km  = 1.852;

% Find something you are looking for
for ii = 1:length(tc);
    names{ii,:} = tc(ii).name;
    year(ii) = tc(ii).time(1);
end
IndexC          = strfind(names,'MATTHEW') 
Index           = find(not(cellfun('isempty', IndexC)));
idhurricane     = 2244;
datestr(year(idhurricane))

% Directory
destout = 'd:\papers\probabilistic_forecasting_tropical_cyclones\figures_hurricanes\';

%% Different hurricanes
for ii = 1:4
    
% Katrina (2005)
if ii == 1;
idhurricane = 2244
fhurricane  = 'Katrina (2005)'
fdir        = 'katrina2005';
tcchoosen   = tc(idhurricane);
landfall    = datenum(2005, 8, 29, 11, 0, 0);
xstart      = datenum(2005, 8, 26, 0, 0, 0);
xend        = datenum(2005, 8, 30, 0, 0, 0);
end

% Ike (2008)
if ii == 2;
idhurricane = 2355
fhurricane  = 'Ike (2008)'
fdir        = 'ike2008';
tcchoosen   = tc(idhurricane);
landfall    = datenum(2008, 9, 13, 7, 0, 0);
xstart      = datenum(2008, 9, 3, 0, 0, 0);
xend        = datenum(2008, 9, 13, 0, 0, 0);
end

% Sandy (2012)
if ii == 3;
idhurricane = 2499
fhurricane = 'Sandy (2012)'
fdir        = 'sandy2012';
tcchoosen   = tc(idhurricane);
landfall    = datenum(2012, 10, 29, 21, 0, 0);
xstart      = datenum(2012, 10, 25, 0, 0, 0);
xend        = datenum(2012, 10, 30, 12, 0, 0);
end

% Matthew (2016)
if ii == 4;
fname       = 'd:\data\hurricanes\archive_server\2016\bal142016.dat'
fhurricane = 'Matthew (2016)'
fdir        = 'matthew2016';
tcchoosen   = tc_read_jtwc_best_track(fname);
landfall    = datenum(2016, 10, 8, 12, 0, 0);
xstart      = datenum(2016, 10, 1, 0, 0, 0);
xend        = datenum(2016, 10, 9, 0, 0, 0);
end

%% Convert units
nt=length(tcchoosen.x);
for it=1:nt     
    
    %% Main
    % Wind speed
    if tcchoosen.vmax(it) ~= -999
        tcchoosen.vmax(it) = tcchoosen.vmax(it) * kts2ms;   % knots to m/s
    end
    % RMAX
    if tcchoosen.rmax(it) ~= -999
        tcchoosen.rmax(it) = tcchoosen.rmax(it) * nm2km;   % nautical miles to km
    end

    %% R35 a
    if tcchoosen.r35nw(it) ~= -999
        tcchoosen.r35nw(it) = tcchoosen.r35nw(it) * nm2km;   % nautical miles to km
    end
    % R35 b 
    if tcchoosen.r35ne(it) ~= -999
        tcchoosen.r35ne(it) = tcchoosen.r35ne(it) * nm2km;   % nautical miles to km
    end
    % R35 c
    if tcchoosen.r35nw(it) ~= -999
        tcchoosen.r35sw(it) = tcchoosen.r35sw(it) * nm2km;   % nautical miles to km
    end
    % R35 d
    if tcchoosen.r35ne(it) ~= -999
        tcchoosen.r35se(it) = tcchoosen.r35se(it) * nm2km;   % nautical miles to km
    end

    %% R50 a
    if tcchoosen.r50nw(it) ~= -999
        tcchoosen.r50nw(it) = tcchoosen.r50nw(it) * nm2km;   % nautical miles to km
    end
    % R50 b 
    if tcchoosen.r50ne(it) ~= -999
        tcchoosen.r50ne(it) = tcchoosen.r50ne(it) * nm2km;   % nautical miles to km
    end
    % R50 c
    if tcchoosen.r50sw(it) ~= -999
        tcchoosen.r50sw(it) = tcchoosen.r50sw(it) * nm2km;   % nautical miles to km
    end
    % R50 d
    if tcchoosen.r50se(it) ~= -999
        tcchoosen.r50se(it) = tcchoosen.r50se(it) * nm2km;   % nautical miles to km
    end
    %% R65 a
    if tcchoosen.r65nw(it) ~= -999
        tcchoosen.r65nw(it) = tcchoosen.r65nw(it) * nm2km;   % nautical miles to km
    end
    % R65 b 
    if tcchoosen.r65ne(it) ~= -999
        tcchoosen.r65ne(it) = tcchoosen.r65ne(it) * nm2km;   % nautical miles to km
    end
    % R65 c
    if tcchoosen.r65nw(it) ~= -999
        tcchoosen.r65sw(it) = tcchoosen.r65sw(it) * nm2km;   % nautical miles to km
    end
    % R65 d
    if tcchoosen.r65ne(it) ~= -999
        tcchoosen.r65se(it) = tcchoosen.r65se(it) * nm2km;   % nautical miles to km
    end

    %% R100 a
    if tcchoosen.r100nw(it) ~= -999
        tcchoosen.r100nw(it) = tcchoosen.r100nw(it) * nm2km;   % nautical miles to km
    end
    % R100 b 
    if tcchoosen.r100ne(it) ~= -999
        tcchoosen.r100ne(it) = tcchoosen.r100ne(it) * nm2km;   % nautical miles to km
    end
    % R100 c
    if tcchoosen.r100sw(it) ~= -999
        tcchoosen.r100sw(it) = tcchoosen.r100sw(it) * nm2km;   % nautical miles to km
    end
    % R100 d
    if tcchoosen.r100se(it) ~= -999
        tcchoosen.r100se(it) = tcchoosen.r100se(it) * nm2km;   % nautical miles to km
    end

end

%% Calculate values
rmax_HURDAT = []; r35_mean_HURDAT = []; dr35_mean_HURDAT = []; pn = []; dpc = []; vmax = [];
for ii = 1:length(tcchoosen.r35se)
   
    % Wind speed
    if tcchoosen.vmax(ii) > 0;
        vmax(ii) = tcchoosen.vmax(ii);
    else
        vmax(ii) = NaN;
    end
    
    % Pressure drop
    if ~isnan(tcchoosen.pressure_last_closed_isobar(ii))
        pn(ii) = tcchoosen.pressure_last_closed_isobar(ii);
    else
        pn(ii) = 1015;
    end
    dpc(ii) = pn(ii) - tcchoosen.pc(ii);
    if dpc(ii) < 0; 
        dpc(ii) = NaN
    end
    
    % Rmax
    if tcchoosen.rmax(ii) > 0
        rmax_HURDAT(ii) = tcchoosen.rmax(ii);
    else
        rmax_HURDAT(ii) = NaN;
    end
    
    % R35
    if tcchoosen.r35se(ii) > 0 & tcchoosen.r35sw(ii) > 0 & tcchoosen.r35ne(ii) > 0 & tcchoosen.r35nw(ii) > 0
        r35_mean_HURDAT(ii) = mean([tcchoosen.r35nw(ii), tcchoosen.r35se(ii), tcchoosen.r35sw(ii), tcchoosen.r35ne(ii)]);
    else
        r35_mean_HURDAT(ii) = NaN;
    end
    
    % Delta R35
    dr35_mean_HURDAT(ii) = r35_mean_HURDAT(ii)-rmax_HURDAT(ii);
end

%% 2. Plot
% General
Cliner = linspecer(6)
[rmax_mean,dr35_mean, rmax_RMSE, r35_RMSE] = rmax_r35(vmax,tcchoosen.y);
r35_mean = dr35_mean + rmax_mean;
destoutTMP = [destout, '\', fdir]; mkdir(destoutTMP); cd(destoutTMP)

% Pressure and wind speed
xtime   = [tcchoosen.time fliplr(tcchoosen.time)];
yvalue  = [rmax_mean-2*rmax_RMSE fliplr(rmax_mean+2*rmax_RMSE)]
id = ~isnan(yvalue);
Y = 29.0/2;   X = 20.5;  Margin = 0;             %# A4 paper size
xSize = X - 2*Margin;   ySize = Y - 2*Margin;        %# figure size on paper (width & height)
hFig = figure; hold on;
set(hFig, 'PaperUnits','centimeters')
set(hFig, 'PaperSize',[X Y])
set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize])
set(hFig, 'PaperOrientation','portrait')
hold on
yyaxis left
h1 = plot(tcchoosen.time,vmax)
ylim([0 max(vmax)*1.1])
set(h1, 'color', Cliner(2,:), 'linewidth', 2);
ylabel('maximum sustained winds [m/s]')

yyaxis right
h2 = plot(tcchoosen.time,dpc)
ylim([0 max(dpc)*1.1])
xlim([xstart xend])
set(h2, 'color', Cliner(1,:), 'linewidth', 2);
ylabel('pressure drop in the eye of the storm [hPa]')

h3 = plot([landfall landfall], [0 500], '--k')
title({'Pressure drop and maximum sustained wind in HURDAT2', ['for hurricane ' fhurricane]})
grid on; box on;
datetick('x','mm/dd','keepticks')
print('-dpng','-r300', ['1_pressure_wind.png']); close all



% Rmax
xtime   = [tcchoosen.time fliplr(tcchoosen.time)];
yvalue  = [rmax_mean-2*rmax_RMSE fliplr(rmax_mean+2*rmax_RMSE)]
id = ~isnan(yvalue);
Y = 29.0/2;   X = 20.5;  Margin = 0;             %# A4 paper size
xSize = X - 2*Margin;   ySize = Y - 2*Margin;        %# figure size on paper (width & height)
hFig = figure; hold on;
set(hFig, 'PaperUnits','centimeters')
set(hFig, 'PaperSize',[X Y])
set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize])
set(hFig, 'PaperOrientation','portrait')
hfill = fill(xtime(id), yvalue(id), 'r'); set(hfill, 'facealpha',.5)
h1 = plot(tcchoosen.time, rmax_HURDAT); set(h1, 'color', Cliner(2,:), 'linewidth', 2);
h2 = plot(tcchoosen.time, rmax_mean); set(h2, 'color', Cliner(1,:), 'linewidth', 2);
h3 = plot([landfall landfall], [0 500], '--k')
legend([h1, h2, hfill h3], 'HURDAT2 database', 'proposed relationship', '95% prediction band', 'landfall', 'location', 'northwest')
grid on; box on;
xlim([xstart xend])
ylim([0 250])
ylabel('radius of maximum winds [km]')
title({'Radius of maximum winds: observed in HURDAT2', 'versus predicted by proposed relationship', ['for hurricane ' fhurricane]})
grid on; box on;
datetick('x','mm/dd','keepticks')
print('-dpng','-r300', ['2_rmax.png']); close all


% R35
xtime   = [tcchoosen.time fliplr(tcchoosen.time)];
yvalue  = [dr35_mean-2*r35_RMSE fliplr(dr35_mean+2*r35_RMSE)]
id = ~isnan(yvalue);
close all; clc
Y = 29.0/2;   X = 20.5;  Margin = 0;             %# A4 paper size
xSize = X - 2*Margin;   ySize = Y - 2*Margin;        %# figure size on paper (width & height)
hFig = figure; hold on;
set(hFig, 'PaperUnits','centimeters')
set(hFig, 'PaperSize',[X Y])
set(hFig, 'PaperPosition',[0.5 0.5 xSize ySize])
set(hFig, 'PaperOrientation','portrait')
hfill = fill(xtime(id), yvalue(id), 'r'); set(hfill, 'facealpha',.5)
h1 = plot(tcchoosen.time, dr35_mean_HURDAT); set(h1, 'color', Cliner(2,:), 'linewidth', 2);
h2 = plot(tcchoosen.time, dr35_mean); set(h2, 'color', Cliner(1,:), 'linewidth', 2);
h3 = plot([landfall landfall], [0 500], '--k')
legend([h1, h2, hfill h3], 'HURDAT2 database', 'proposed relationship', '95% prediction band', 'landfall', 'location', 'northwest')
grid on; box on;
xlim([xstart xend])
ylabel('average difference radius of 35 knots [km]')
title({'Average radius of 35 knots: observed in HURDAT2' ,'versus predicted by proposed relationship', ['for hurricane ' fhurricane]})
grid on; box on;
datetick('x','mm/dd','keepticks')
print('-dpng','-r300', ['3_r35.png']); close all


end