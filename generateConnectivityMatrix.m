function [C, total_corr, sorted] = generateConnectivityMatrix(Ngroups, Ncells, p_same_group, p_diff_group)
    % Generate random group IDs for each cell
    group_vals = ceil(rand(1, Ncells) * Ngroups);

    % Initialize the connectivity matrix
    C = zeros(Ncells);
    % Calculate the difference between same-group and different-group probabilities
    p_diff = p_same_group - p_diff_group;
    % Loop through each cell to determine connections
    for i = 1:Ncells
        C(i, :) = rand(1, Ncells) < (p_diff_group + p_diff * (group_vals == group_vals(i)));
    end

    row_corr = corr(C');        % Row to row correlations
    col_corr = corr(C);         % Column to column correlations
    total_corr = col_corr + row_corr;   % Sum of column and row correlations

    threshold = 0.52;   

    groupid = 0; % This will count through the groups. Initialize it to zero                      
    grouped_matrices = zeros(1, Ncells);    % This will label the group of each cell
    
    for i = 1:Ncells
        if grouped_matrices(i) == 0
            groupid = groupid + 1; % Increment group ID for a new group
            remaining_list = i; % Initialize the remaining list with the current cell
            used_list = []; % Initialize the used list
            
            while ~isempty(remaining_list)
                cell = remaining_list(1); % Pop the first cell from the remaining list
                remaining_list(1) = []; % Remove the first cell from the remaining list
                
                groupcells = find(total_corr(cell, :) > threshold); % Find correlated cells
                grouped_matrices(groupcells) = groupid; % Assign group ID to these cells
                
                used_list = union(used_list, cell); % Add current cell to used list
                new_cells = setdiff(groupcells, used_list); % Find new cells to be processed
                remaining_list = union(remaining_list, new_cells); % Update remaining list
            end
        end
    end
    [~, neworder] = sort(grouped_matrices);

    sorted = C(neworder, neworder);


end
