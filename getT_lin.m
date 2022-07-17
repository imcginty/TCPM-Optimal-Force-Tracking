function [T_min T_max] = getT_lin()
% Returns range of temps the thermomechanical model is linearized about
T_min = 110; % [ºC] minimum of linear region
T_max = 160; % [ºC] maximum of linear region
end