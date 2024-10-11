%% Final Project: 
% Analyzing different initial conditions in motif extraction to identify patterns in distribution of multiplicity

%% POD Summary

% Group Name: Motif Gang
% Group Members: Pranav Daggubati, Lyla Liu

%% Code
clear
close all

%% Creating Initial Variables:

% Parameters
Ngroups = 5;                                   % Number of distinct groups
Ncells = 200;                                   % Total number of cells to consider
group_vals = ceil(rand(1, Ncells) * Ngroups);   % Randomly assign group IDs

p_same_group_values = 0:0.1:1;                 % Array of same-group connection probabilities
p_diff_group_values = 0:0.1:1;                 % Array of different-group connection probabilities

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

%% Creating Random Initial Variables Connectivity Matrices

% Initialize a cell array to hold the connectivity matrices
connectivity_matrices = cell(1, num_pairs);
sorted_conn_mats = cell(1, num_pairs);

% Initialize a cell array to hold correlation matrices
corr_list = cell(1, num_pairs);

% Select a threshold for same-group membership, based partly on the
% distribution observed in the histogram. This can be adjusted. A higher
% threshold results in more distinct groups being extracted.
threshold = 0.5;

% Loop through the lists of probability values and generate connectivity matrices
for i = 1:num_pairs
    p_same_group = p_same_group_list(i);
    p_diff_group = p_diff_group_list(i);
    [C, total_corr, sorted] = generateConnectivityMatrix(Ngroups, Ncells, p_same_group, p_diff_group);
    
    connectivity_matrices{i} = C;
    corr_list{i} = total_corr;   % Adding corr matrix of each initial condition to 
    sorted_conn_mats{i} = sorted;
end

motif_ratio_matrix = cell(1, num_pairs);

for i = 1:num_pairs
    [unidirectProb, bidirectProb, numConnections, numBiConnections, numTriples, expectedTriplets] = tripletMotifs(connectivity_matrices{i});
    
    % All motifs involving three cells found in a connectivity matrix
    % ratioed with expected number of that motif
    motif_ratio_list = numTriples./expectedTriplets;
    
    % Add list to cell matrix
    motif_ratio_matrix{i} = motif_ratio_list;
end


%% Test Plots

% Plot initial conditions as a scatter plot
figure;
scatter(p_same_group_list, p_diff_group_list);
title('Same vs Diff Group Connection Prob')
xlabel('Same P');
ylabel('Diff P');

% Display the first connectivity matrix as an example
figure;
subplot(1,2,1);
imagesc(connectivity_matrices{50});
subplot(1,2,2);
histogram(corr_list{50});