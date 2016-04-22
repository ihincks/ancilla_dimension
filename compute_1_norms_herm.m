%%  compute_1_norms_herm
% This script is identical to compute_1_norms, but calls
% InducedSchattenModHerm instead of InducedSchattenMod. See compute_1_norms
% for an explanation of how these scripts work.

% The list of n values we wish to compute the norm for. The computation
% will run for m from n, ..., n^2
n_values = [2,3,4,5,6];
%n_values=[2];
% The number of initial guesses for computing the lower bounds
%N=2048;
%N = 1024;
N = 1000;
% The number of computations to do before saving
segments = 15*4;

% This sets the number of times that a single core can suffer SVD
% computation failure before the program quits. I'm arbitrarily setting
% this at 20. It seems unlikely that anything a single core can hit this
% many failures in a row, but if it does it should definitely quit.
max_attempts = 20;

save_file = 'norm_bounds_herm.mat';

% A stupid way of checking if the file already exists, if it does load the
% existing results, and if not create a new structure.
% The structure has the form:
%   results
%       - n# - for each n value, the results for each m value are stored in
%              a struct called n# (e.g. n6 for n=6).
%           - m# - similarly, for each m value for a given n, results are
%               stored as m#
dirlist = dir('.');
dirlist = {dirlist.name};
if sum(ismember(dirlist,save_file))
    results = load(save_file);
else
    results = struct();
end

% loop over values of n
for nidx=1:length(n_values)
    % the string 'n#'
    n = n_values(nidx);
    ns = ['n', num2str(n)];
    
    % check to see if the results struct already has an entry for this n
    % value
    if not(isfield(results, ns))
        results.(ns) = struct();
    end
    
    % for each n value, loop through the desired m values
    for m = n:(n^2)
        % before doing anything, check to see if there is already a stored
        % result for the m value. If the field already exists do nothing
        ms = ['m',num2str(m)];
        
        % If the field for the current value of n and m doesn't exist,
        % create and initialize it.
        if not(isfield(results.(ns),ms))
            results.(ns).(ms) = struct();
            results.(ns).(ms).val = 0;
            results.(ns).(ms).op = 0;
            results.(ns).(ms).runs = 0;
        end
        
        % If the currently stored number of runs factored into the
        % currently stored best value for the current value of n and m is
        % less than that specified at the beginning of the script, run more
        % initial guesses.
        if results.(ns).(ms).runs < N
            % run intial guesses until the number of runs exceeds those
            % specified at the beginning of the program
            while results.(ns).(ms).runs < N
                % compute how many runs are needed, and if it is more than
                % the value of segments, reduce it to that value
                runs_to_do = N - results.(ns).(ms).runs;
                if runs_to_do > segments
                    runs_to_do = segments;
                end
                
                % compute the new values
                values = cell(runs_to_do,1);
                ops = cell(runs_to_do,1);
                parfor k=1:runs_to_do
                    success =0;
                    attempts = 0;
                    while not(success)
                        % This error catching step is inserted as on occasion there
                        % is a problem with computing SVDs. It seems to happen very
                        % rarely (I only ever run into it when running many initial
                        % guesses), so the errors are caught, and the program
                        % only throws the error if it happens some threshold
                        % number of times in a row.
                        try
                            [nrm,X] = InducedSchattenNormModHerm(@(x)trash_transpose_map(x, n, m), @(x)trash_transpose_map_adjoint(x,n,m), n*n*m, 2*n*m);
                            values{k} = nrm;
                            ops{k} = X;
                            success = 1;
                        catch SVDE
                            attempts = attempts + 1;
                            %quit if the number of attempts hits max_attempts
                            if attempts == max_attempts
                                rethrow(SVDE)
                            end
                        end
                    end
                end
                
                % get the best value and its index
                value_list = zeros(runs_to_do,1);
                for k = 1:runs_to_do
                    value_list(k) = values{k}(1);
                end
                [y,i] = max(value_list);
                
                % If the found value is better than the stored one, replace
                % it
                if y > results.(ns).(ms).val
                    results.(ns).(ms).val = y;
                    results.(ns).(ms).op = ops{i};
                end
                
                %update the number of runs that have been factored into the
                %currently stored best value
                results.(ns).(ms).runs = results.(ns).(ms).runs + runs_to_do; 
                
                % save the results
                save(save_file, '-struct', 'results');
                display(['--- case n=', num2str(n), ', m=', num2str(m), ', total runs = ', num2str(results.(ns).(ms).runs),', at ', datestr(datetime('now'))]);
            end % end of while loop
            % display the time that the computation for the case n and m
            % have finished
            c = clock;
            display(['completed case n=', num2str(n), ', m=', num2str(m), ', at ', datestr(datetime('now'))]);
        else
            display(['case n=', num2str(n), ', m=', num2str(m), ' already done']);
        end
    end
end