% Work in progress - to share parameters across classes.

classdef SharedParams
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
    function nargout = parse_inputs(params, varargin)
            parser = inputParser();
            
            %Parameters used for f(t) and rho(t)
            parser.addParameter('startpulse', 200); 
            parser.addParameter('lengthpulse', 200);
            parser.addParameter('lengtht1', 10);

            parser.parse(varargin{:})
            params = parser.Results;
            params.t_0 = params.startpulse;
            params.t_1 = params.t_0 + params.lengtht1;
            params.t_2 = params.t_0 + params.lengthpulse;
    end
    end
    
end

% There may be other shared parameters? Check once code is working