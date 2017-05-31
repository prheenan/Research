function []=feather_example()
    %{
    example which uses FEATHER.

    Path must be set properly 
    %}
    base = '/Users/patrickheenan/src_prh/Research/Perkins/Projects/';
    base_path = [base,'PythonCommandLine/FEATHER/'];
    % read the input file
    input_csv = 'example.csv';
    data = csvread(input_csv,2,0);
    % get the individual columns, for plotting purposes
    time = data(:,1);
    separation = data(:,2);
    force = data(:,3);
    % in this case, the instrument records  meta information we need 
    trigger_time = 0.382;
    dwell_time = 0.992;
    spring_constant = 6.67e-3;
    % get the force extension curve object to use
    obj = fec(time,separation,force,trigger_time,dwell_time,...
              spring_constant);
    % get the feather-specific options to use
    threshold = 1e-3;
    tau = 2e-2;
    opt = feather_options(threshold,tau,base_path);
    % get the predicted event locations
    indices = feather(obj,opt); 
    disp(indices)
    clf;
    hold all;
    plot(obj.time,obj.force*-1)
    for i=1:length(indices)
        plot(obj.time(indices(i)),obj.force(indices(i))*-1,'ro')
    end
        
end
