function ratio_cell = motifRatios(groups, cells, p_same, p_diff)
    
    % STEP 1: Creating Initial Variables:
    
    % Parameters
    Ngroups = groups; % Number of distinct groups
    Ncells = cells; % Total number of cells to consider
    group_vals = ceil(rand(1, Ncells) * Ngroups); % Randomly assign group IDs
    
    p_same_group_values = p_same; % Array of same-group connection probabilities
    p_diff_group_values = p_diff; % Array of different-group connection probabilities
    
    % Ensure the input arrays are of the same length
    if length(p_same_group_values) ~= length(p_diff_group_values)
        error('The length of p_same_group_values and p_diff_group_values must be the same.');
    end
    
    % Generate all possible pairs of probabilities
    [p_same_group_list, p_diff_group_list] = meshgrid(p_same_group_values, p_diff_group_values);
    
    % Flatten the matrices into vectors
    p_same_group_list = p_same_group_list(:)';
    p_diff_group_list = p_diff_group_list(:)';
    
    % Number of pairs of probabilities generated
    num_pairs = length(p_same_group_list);

    % STEP 2: Make Conn Matrices

    % Initialize a cell array to hold the connectivity matrices
    connectivity_matrices = cell(1, num_pairs);
    sorted_conn_mats = cell(1, num_pairs);
    
    % Initialize a cell array to hold correlation matrices
    corr_list = cell(1, num_pairs);
    
    % Loop through the lists of probability values and generate connectivity
    % matrices
    for i = 1:num_pairs
        p_same_group = p_same_group_list(i);
        p_diff_group = p_diff_group_list(i);
        [connectivity_matrices{i}, corr_list{i}, sorted_conn_mats{i}] = generateConnectivityMatrix(Ngroups, Ncells, p_same_group, p_diff_group);
        
    end

    % STEP 3: Get ratios of motif classes and motifs

    motif_ratio_matrix = cell(1, num_pairs);
    class_ratio_matrix = cell(1, num_pairs);
    
    for i = 1:num_pairs
        [motif_expected, motif_empirical, class_expected, class_empirical] = tripletMotifs(connectivity_matrices{i});
        
        % All motifs involving three cells found in a connectivity matrix
        % ratioed with expected number of that motif
        class_ratio_list = class_empirical./class_expected;
        motif_ratio_list = motif_empirical./motif_expected;
    
        % Add list to cell matrix
        motif_ratio_matrix{i} = motif_ratio_list;
        class_ratio_matrix{i} = class_ratio_list;
    end

    % STEP 4: Analysis

    % Initialize an empty array to store the data
    numCells = length(class_ratio_matrix);
    numClasses = length(class_ratio_matrix{1});
    classdata = zeros(numCells, numClasses);
    
    numMotifs = length(motif_ratio_matrix{1});
    motifdata = zeros(numCells, numMotifs);
    
    % Extract the ratio vectors and store them in the array
    for i = 1:numCells
        classdata(i, :) = class_ratio_matrix{i};
        motifdata(i, :) = motif_ratio_matrix{i};
    end

    ratio_cell = cell(1,2);
    ratio_cell{1} = motifdata;
    ratio_cell{2} = classdata;