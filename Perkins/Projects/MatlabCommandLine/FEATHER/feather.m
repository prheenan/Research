function []=feather()
    base = '/Users/patrickheenan/src_prh/Research/Perkins/Projects/';
    base_path = [base,'PythonCommandLine/FEATHER/'];
    input_file = [base_path,'main_feather.py'];
    % read the input file
    [status,cmdout] =system('pwd');
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
    surface_position_from_approach = 0;
    opt = feather_options(threshold,tau,surface_position_from_approach);
    % save out the fec 
    struct_to_save = {};
    struct_to_save.time = obj.time;
    struct_to_save.separation = obj.separation;
    struct_to_save.force = obj.force;
    % XXX approach, retract velocities, spring constant
    matlab_file = [base_path,'tmp.mat'];
    output_file = [base_path,'out.csv'];
    save(matlab_file,'-struct','struct_to_save','-v7.3');
    command = ['//anaconda/bin/python2.7 ',input_file];
    command = append_numeric(command,...
                             '-spring_constant',obj.spring_constant);
    command = append_numeric(command,...
                             '-trigger_time',obj.trigger_time);   
    command = append_numeric(command,...
                             '-dwell_time',obj.dwell_time);    
    command = append_numeric(command,...
                             '-threshold',opt.threshold);      
    command = append_numeric(command,...
                             '-tau',opt.tau);    
    command = append_argument(command,...
                              '-file_input',matlab_file);                          
    command = append_argument(command,...
                              '-file_output',output_file);                             
    [status,cmdout] = system(command);
    disp(command);
    disp(cmdout);
end

function[output]=append_numeric(output,name,value)
    output = append_argument(output,name,num2str(value,'%15.3g'));
end

function[output]=append_argument(output,name,value)
    output = [output,' ',name,' ',value];
end