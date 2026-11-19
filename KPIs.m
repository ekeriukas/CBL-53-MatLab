% Define static helper class for KPI functions
% (used to access multiple functions defined from outside of the script)
% Use functions like this: myEff = KPIs.efficiency(myPb, myLHV, mymdot_f);
classdef KPIs
methods (Static)

function result = efficiency(Pb, LHV, mdot_f)
    % Efficiency [-] as a function of
    % Pb       : break power [kW]
    % LHV      : lower heating value [MJ/kg]
    % mdot_f   : fuel mass flow [kg/s]
    result = (Pb)/(LHV*1e3*mdot_f);
end

function result = bsfc(Pb, mdot_f)
    % Break specific fuel consumption [g/kWhr] as a function of
    % Pb       : break power [kW]
    % mdot_f   : fuel mass flow [kg/s]
    result = mdot_f/Pb*(3.6e6);
end

function result = bsCO2(Pb, mdot_exh, Mexh, MCO2, XCO2)
    % Break specific CO2 emission [g/kWhr] as a function of
    % Pb       : break power [kW]
    % mdot_exh : exhaust mass flow [kg/s]
    % Mexh     : molar mass of exhaust [g/mol]
    % MCO2     : molar mass of CO2 [g/mol]
    % XCO2     : CO2 volume fraction in exhaust [vol %]
    YCO2 = (XCO2*1e-2)*MCO2/Mexh; % CO2 mass fraction in exhaust [-]
    mdot_CO2 = mdot_exh * YCO2;
    result = mdot_CO2/Pb*(3.6e6);
end

function result = bsNOx(Pb, mdot_exh, Mexh, MNOx, XNOx)
    % Break specific NOx emission [g/kWhr] as a function of
    % Pb       : break power [kW]
    % mdot_exh : exhaust mass flow [kg/s]
    % Mexh     : molar mass of exhaust [g/mol]
    % MNOx     : equivalent molar mass of NOx [g/mol]
    % XNOx     : NOx volume fraction in exhaust [vol ppm]
    YNOx = (XNOx*1e-6)*MNOx/Mexh; % NOx mass fraction in exhaust [-]
    mdot_NOx = mdot_exh * YNOx;
    result = mdot_NOx/Pb*(3.6e6);
end

end
end