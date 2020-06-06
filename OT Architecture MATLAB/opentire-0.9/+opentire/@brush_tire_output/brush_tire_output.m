classdef (Sealed = true) brush_tire_output < opentire.tire_output

properties
% Additional outputs for belt deflections, pressure, 
p = []; % pressure/force
e = []; % element tip deflections
b = []; % belt coordinates
end % properties

end % classdef

