classdef mf61 < tire
    % tire object with MF6.1 implementation
    
    properties
        
        % LATERAL_COEFFICIENTS
        PCY1                    % Shape factor Cfy for lateral forces
        PDY1                    % Lateral friction Muy
        PDY2                    % Variation of friction Muy with load
        PDY3                    % Variation of friction Muy with squared camber
        PEY1                    % Lateral curvature Efy at Fznom
        PEY2                    % Variation of curvature Efy with load
        PEY3                    % Zero order camber dependency of curvature Efy
        PEY4                    % Variation of curvature Efy with camber
        
        % ...etc
        % ...etc
        
    end
    
    methods
        
        % Class constructor
        function obj = mf61()
            
            % Assign the FITTYP
            obj.FITTYP = '61';
            
            % do more stuff
            
        end
        
        function output = special61method(input)
            % some special function specific to MF 6.1 implementation
            
        end
        
    end
end

