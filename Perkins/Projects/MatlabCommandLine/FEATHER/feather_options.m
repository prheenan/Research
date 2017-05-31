classdef feather_options
   properties
        threshold
        tau
        surface_position_from_approach
   end
   methods
      function obj = feather_options(threshold,tau,...
                                     surface_position_from_approach)
            obj.threshold = threshold;
            obj.tau = tau;
            obj.surface_position_from_approach = ...
                surface_position_from_approach;
      end
   end
end
