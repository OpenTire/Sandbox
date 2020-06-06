classdef tire < handle
    %TIRE represents a tire model
    
    properties
        
        % UNITS (with defaults)
        lengthUnit = 'meter'
        forceUnit = 'newton'
        angleUnit = 'radians'
        massUnit = 'kg'
        timeUnit = 'second'
        
        % MODEL
        FITTYP                  % Magic Formula Version number
        TYRESIDE
        LONGVL                  % Nominal speed
        VXLOW                   % Lower boundary of slip calculation
        ROAD_INCREMENT          % Increment in road sampling
        ROAD_DIRECTION          % Direction of travelled distance
        
        % DIMENSION
        UNLOADED_RADIUS         % Free tyre radius
        WIDTH                   % Nominal section width of the tyre
        ASPECT_RATIO            % Nominal aspect ratio
        RIM_RADIUS              % Nominal rim radius
        RIM_WIDTH               % Rim width
        
        % INERTIA
        MASS                    % Tyre Mass
        IXX                     % Tyre diametral moment of inertia
        IYY                     % Tyre polar moment of inertia
        BELT_MASS               % Belt mass
        BELT_IXX                % Belt diametral moment of inertia
        BELT_IYY                % Belt polar moment of inertia
        GRAVITY                 % Gravity acting on belt in Z direction
        
        %...etc
        %...etc
        
        
    end
    
    methods
        
        % Define the class constructor
        function obj=tire(varargin)
            
            % Check if no arguments are passed to the constructor
            if nargin==0
                
                % This allows the user to instantiate a tire
                % without specifying an input variable.
                return
                
            else
                
                % do something else with the input (auto run import TIR
                % function if TIR file is passed?)
                
            end
            
        end
        
        function out = ImportTIR(tirfile)
            % Fill a tire object with data from a TIR file
            
            
        end
        
        function WriteTIR(writelocation)
            % Write TIR file with tire object properties
            
        end
        
        % Compute functions
        function output = compute_steady(input, options)
            % compute_steady method
            % Computes steady-state forces, moments, and other tire_output values.
            % Some tire_input fields are optional, indicated by the value nan.
            
            % call steady state solver
            
        end
        
        
        function [output,state_dots,algebraic_loops] = compute_dynamic(time, input, road, states, algebraic_loops)
            % compute_dynamic method
            % algebraic_loops is [] or a vector of algebraic-loop values (outputs that are also inputs).
            % Integrators should call this function with values from previous time step.
            % Equilibrium solvers should solve for (output algebraic_loops - input_algebraic_loops) == 0.
            % states is [] or a vector of states whose time-derivatives are computed by the tire model.
            % state time derivatives are returned in state_dots.
            
            
            % call dynamic solver
            
        end
        
        
        
    end
    
    methods(Static)
        
        
        
    end
    
end

