function cmd_output = cmd_history(cmd_name, cmd_parameters)
    %% Print command history to screen
    % Input:
    %   cmd_name: Command name;
    %   cmd_parameters: Row vactor containing parameters.
    % Output:
    %   cmd_output: String "'command'.mat".

    [~, cmd_parameters_length] = size(cmd_parameters);

    cmd_output = strcat(cmd_name, '( ');

    ii = 1;

    while ii < cmd_parameters_length
        cmd_output = strcat(cmd_output, num2str(cmd_parameters(ii)), ', ');
        ii = ii + 1;
    end

    cmd_output = strcat(cmd_output, num2str(cmd_parameters(ii)), ')');

    fprintf('\n\n\n%s', 'Command launched:');
    fprintf('\n%s\n\n', cmd_output);
    fprintf('%s\n', 'Launch time:');
    disp(datetime);

    cmd_output = strcat(cmd_output, '.mat');

end % End of cmd_history()
