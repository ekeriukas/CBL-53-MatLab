%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 
%% init
clear all; clc;close all;
addpath( "Functions","Nasa");
%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv
Runiv = 8.314;
[SpS,El]        = myload('Nasa\NasaThermalDatabase.mat',{'Diesel','O2','N2','CO2','H2O'});

%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 4.2;
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355;
CaIVC = -135;
CaEVO = 149;
CaEVC = -344;
CaSOI = -3.2;
% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams


%% Load data (if txt file)
FullName        = fullfile('Data','ExampleDataSet.txt');
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p               = reshape(dataIn(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
temp            = reshape(dataIn(:,3),[],Ncycles);


%% 1. Settings you can tune
gamma     = 1.35;     % constant gamma for a first result (see step 7 to make it vary with Ca)
sgOrder   = 3;        % Savitzky-Golay polynomial order
sgFrame   = 41;       % window in samples (must be ODD). 41 * 0.2 deg = 8.2 deg. Try 21 ... 61
pIntake   = [];       % [Pa] absolute intake pressure for pegging. Leave empty = no pegging
CaWinHR   = [-10 90]; % window used to integrate the heat release (just before combustion to ~end)

%% 2. One crank-angle vector, geometry and its derivative
Ca1  = Ca(:,1);                          % identical for every cycle, size (3600,1)
dCa  = Ca1(2) - Ca1(1);                  % 0.2 deg
V    = CylinderVolume(Ca1, Cyl);         % [m^3]
dVdCa = gradient(V, Ca1);                % [m^3/deg]  (analytic derivative is also possible)

%% 3. Cycle-averaged pressure (+ optional pegging)
pm = mean(p, 2);                         % ensemble mean over the 69 cycles [Pa]
if ~isempty(pIntake)
    [~, iBDC] = min(abs(Ca1 + 180));     % intake BDC
    pm = pm + (pIntake - pm(iBDC));      % shift so that p(BDC) = measured intake pressure
end

%% 4. Filter the pressure, then differentiate
% Averaging is not enough for a derivative, so smooth first (handbook Appendix 1).
ps    = sgolayfilt(pm, sgOrder, sgFrame);
dpdCa = gradient(ps, Ca1);               % [Pa/deg]

%% 5. aROHR  [J/deg]
%  aROHR = gamma/(gamma-1) * p * dV/dCa + 1/(gamma-1) * V * dp/dCa
aROHR = gamma/(gamma-1) .* ps .* dVdCa + 1/(gamma-1) .* V .* dpdCa;

%% 6. Cumulative heat release (aHR [J]) and CA10 / CA50 / CA90
iw   = Ca1 >= CaWinHR(1) & Ca1 <= CaWinHR(2);
Caw  = Ca1(iw);
aHR  = cumtrapz(Caw, aROHR(iw));         % [J]
aHR  = aHR - aHR(1);                     % zero at the start of the window
aHRtot = max(aHR);                       % total apparent heat released
xb   = aHR / aHRtot;                     % burned fraction 0..1

CAx = @(frac) interpAt(Caw, xb, frac);
CA10 = CAx(0.10);  CA50 = CAx(0.50);  CA90 = CAx(0.90);
fprintf('aHR_total = %.0f J | CA10 = %.1f | CA50 = %.1f | CA90 = %.1f deg\n', aHRtot, CA10, CA50, CA90);
fprintf('Ignition delay (CA10 - SOI) = %.1f deg = %.2f ms at 1500 rpm\n', ...
        CA10 - CaSOI, (CA10 - CaSOI)/(360*1500/60)*1e3);

%% 7. Plots (same layout as Figure 2.7 of the handbook)
f4 = figure(4); set(f4,'Position',[200 300 1200 400]);
subplot(1,2,1)
plot(Ca1, aROHR, 'LineWidth', 1.5); hold on
line([0 0], ylim, 'Color','k'); grid on
xlim([-45 135]); xlabel('Crank angle [deg]'); ylabel('aROHR [J/CA]');
subplot(1,2,2)
plot(Caw, aHR, 'LineWidth', 1.5); hold on
plot([CA10 CA50 CA90], aHRtot*[0.1 0.5 0.9], 'o', 'MarkerFaceColor','y');
grid on; xlabel('Crank angle [deg]'); ylabel('aHR [J]');

%% Local helper (put at the very END of the script file, or in its own file interpAt.m)
function ca = interpAt(Caw, xb, frac)
    i = find(xb >= frac, 1, 'first');            % first sample at/above the fraction
    ca = interp1(xb(i-1:i), Caw(i-1:i), frac);   % linear interpolation between the 2 samples
end



%% Plotting 

%% In-cylinder pressure
% All cycles in one plot
f1=figure(1);
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
subplot(1,2,1)
pp = plot(Ca,p/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')

% Filtered mean pressure in a separate plot
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
subplot(1,2,2)
pp = plot(Ca,ps/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('Filtered mean of all cycle pressures.')

%% pV-diagram
V = CylinderVolume(Ca(:,iselect),Cyl);
f3 = figure(2);
set(f3,'Position',[ 200 400 600 800]);              % Just a size I like. Your choice
subplot(2,1,1)
plot(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0 0.8]);ylim([0.5 50]);                      % Matter of taste
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram','(with new Volume function)'})
subplot(2,1,2)
loglog(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0.02 0.8]);ylim([0 50]);                      % Matter of taste
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram','(with new Volume function)'})
