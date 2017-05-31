classdef feather_options
   properties
        threshold
        tau
        base_path
   end
   methods
      function obj = feather_options(threshold,tau,base_path)
            %{
            constructor for object

            Args:
                threshold: probability threshold which feather uses (0,1)
                tau: fractional smoothing (0,1)         
                base_path: location of the python code
            Returns:
                constructed feather options object
            %} 
            obj.threshold = threshold;
            obj.tau = tau;
            obj.base_path = base_path;
      end
   end
end
