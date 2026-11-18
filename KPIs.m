% Define static helper class for KPI functions, as apparently there is no
% other clean way to access multiple functions defined in a script from
% outside of the script.
% Use functions like this: myEff = KPIs.efficiency(myPb, myLHV, myMdot_f);
classdef KPIs
methods (Static)

function result = efficiency(Pb, LHV, Mdot_f)
    % Efficiency [-] as a function of
    % Pb : break power [kW]
    % LHV : lower heating value [MJ/kg]
    % Mdot_f : fuel mass flow [kg/s]
    result = (Pb*1e3)/(LHV*Mdot_f);
end

function result = bsfc(Pb, Mdot_f)
    % Break specific fuel consumption [g/kWhr] as a function of
    % Pb : break power [kW]
    % Mdot_f : fuel mass flow [kg/s]
    result = Mdot_f/Pb*(3.6e6);
end

function result = bsCO2(Pb, Mdot_CO2)
    % Break specific CO2 emission [g/kWhr] as a function of
    % Pb : break power [kW]
    % Mdot_f : CO2 mass flow [kg/s]
    result = Mdot_CO2/Pb*(3.6e6);
end

function result = bsNOx(Pb, Mdot_NOx)
    % Break specific NOx emission [g/kWhr] as a function of
    % Pb : break power [kW]
    % Mdot_f : NOx mass flow [kg/s]
    result = Mdot_NOx/Pb*(3.6e6);
end

end
end