function [motif_expected, motif_empirical, class_expected, class_empirical] = tripletMotifs(conn_matrix)
    % motif_tricounts counts numbers of triplet connections vs chance
    %
    %   [expected, empirical] = tripletMotifs(conn_matrix)
    %
    %   The function requires as inputs:
    %   conn_matrix is the input square connectivity matrix
    %
    %   It returns the following set of statistics:
    %   uniP is mean connection probability (i.e. fraction of connections that are present)
    %   biP is empirical probability of a bidirectional connection (i.e. the fraction of cell-pairs connected bidirectionally)
    %   Nconns is total number of connections
    %   Nbi is double the total number of bidirectionally connected pairs (i.e. the number of connections)
    %   Ntriples is the actual number of each motifs for triplets of cells
    %   N3expect is the expected numbers of each motif for the connection patterns between three different neurons

    % Ensure the matrix is square
    [N1, N2] = size(conn_matrix);
    if N1 ~= N2
        error('Error! The input connectivity matrix is not square.');
    end

    % Total number of connections
    Nconns = nnz(conn_matrix);
    % Mean connection probability
    uniP = Nconns / (N1 * N2);

    % Number of self-connections
    Nself = sum(diag(conn_matrix) > 0);
    selfP = Nself / N1;
    
    % Count bidirectional connections
    biConnMatrix = conn_matrix & conn_matrix';
    Nbi = sum(biConnMatrix(:)) / 2;
    biP = 2 * Nbi / Nconns;

    % Number of potential pairs
    Npairs = N1 * (N1 - 1) / 2;
    biProb = Nbi / Npairs;
    uniProb = (Nconns - Nself - 2 * Nbi) / Npairs /2;
    zeroProb = 1 - biProb - uniProb;

    % Expected number of motifs
    motif_expected = zeros(1, 16);
    motif_expected(1) = zeroProb^3;
    motif_expected(2) = 6 * zeroProb^2 * uniProb;
    motif_expected(3) = 3 * zeroProb^2 * biProb;
    motif_expected(4) = 3 * zeroProb * uniProb^2;
    motif_expected(5) = 3 * zeroProb * uniProb^2;
    motif_expected(6) = 6 * zeroProb * uniProb^2;
    motif_expected(7) = 6 * zeroProb * uniProb * biProb;
    motif_expected(8) = 6 * zeroProb * uniProb * biProb;
    motif_expected(9) = 3 * zeroProb * biProb^2;
    motif_expected(10) = 6 * uniProb^3;
    motif_expected(11) = 2 * uniProb^3;
    motif_expected(12) = 3 * uniProb^2 * biProb;
    motif_expected(13) = 6 * uniProb^2 * biProb;
    motif_expected(14) = 3 * uniProb^2 * biProb;
    motif_expected(15) = 6 * uniProb * biProb^2;
    motif_expected(16) = biProb^3;
    motif_expected = motif_expected * N1 * (N1 - 1) * (N1 - 2) / 6;
    
    % Initialize Ntriples
    motif_empirical = zeros(1, 16);
    
    for i = 1:N1
        for j = i + 1:N1
            type1 = conn_matrix(i,j) + 2*conn_matrix(j,i);
            for k = j + 1:N1
                type2 = conn_matrix(j,k) + 2*conn_matrix(k,j);
                type3 = conn_matrix(k,i) + 2*conn_matrix(i,k);
                type = [type1, type2, type3];
                
                % Number of connections in triplet
                % Counts if type > 0 which means any two nodes connected
                % anywhere
                num = sum(type > 0);
    
                % Num2 conditional to check number of connections that are
                % bidirectional
                num2 = sum(type == 3);
    
                % Motif numbering is based on the multiplicity defined
                % above in the N3expect array.
                % Motifs are defined out of order below to optimize conditionals. 
    
                % No Connections in Triplet
                if num == 0
                    % Motif 1: No connections in the triplet
                    motif = 1;
    
                % One Connection in Triplet
                elseif num == 1
                    % Motif 2: One unidirectional connection
                    % Motif 3: One bidirectional connection
                    motif = 2 + (num2 == 1);    % Boolean here to check if it is a 
                                                % bidirectional connection
                
                % Two Connections in Triplet
                elseif num == 2
                    if num2 == 2 % Two bidirecitonals

                        % Motif 9: L-shaped both bidirectional connections
                        motif = 9;

                    elseif num2 == 1 % One bidirectional
                        % Motifs 7-8: L-shaped mixed connections
                            % Motif 7: One unidirectional connection to one of 
                            % the bidirectionally connected cells
                            % Motif 8: One unidirectional connection away from 
                            % one of bidirectionally connected cells

                        if sum(type) == 4  % must be a "3" and a "1"
                            type2 = min(type, 2);
                            if mod(type2(2) - type2(1), 3) == 1
                                motif = 7;
                            else
                                motif = 8;
                            end
                        else  % sum(type) = 5, must be a "3" and a "2"
                            type2 = max(type, 1);
                            if mod(type2(1) - type2(2), 3) == 1
                                motif = 7;
                            else
                                motif = 8;
                            end
                        end

                    else % No bidirectional in two connections
                        % Motifs 4-6: L-shaped unidirectional connections
                            % Motif 4: Two unidirectional connections with same start point
                            % Motif 5: Two unidirectional connections with same endpoint
                            % Motif 6: Two unidirectional connections in a chain
                        
                        % Which pairs are connected and directionality
                        directionality = find(type>0);
                        if( type(directionality(1)) == type(directionality(2)) ) 
                            motif = 6;
                        elseif ( mod(type(2)-type(1),3) == 1 ) 
                            motif = 5;
                        else
                            motif = 4;
                        end
                    end
    
                % Three Connections in Triplet
                else
                    if num2 == 3
                        % Motif 16: All bidirectional connections
                        motif = 16;
                    elseif num2 == 2
                        % Motif 15: Two bidirectional connections and one unidirectional connection
                        motif = 15;
                    elseif num2 == 1
                        % Motif 12-14: One bidirectional connection and two unidirectional connections
                            % Motif 12: Both unidirectional conenctions point toward the bidirectionally connected pair
                            % Motif 13: One unidirectional points away, the other toward the bidirectionally connected pair
                            % Motif 14: Both unidirectional connections point away from bidirectionally connected pair  
                        if sum(type) == 6
                           if ( mod(type(2)-type(1),3) == 1 )
                               motif = 14;  
                           else
                               motif = 12;  
                           end                                        
                        else
                            motif = 13; 
                        end
                    else
                        % Motif 10-11: Three unidirectional connections
                            % Motif 10: Two connections in chain with third going from start node to end node
                            % Motif 11: Connections in circle
                        motif = 10 + (mod(sum(type), 3) == 0);
                    end

                end
                motif_empirical(motif) = motif_empirical(motif) + 1;
            end
        end
    end
    
    % Expected Values for 0 connections, purely uni connections, mixed
    % connections, and purely bi connections
    expZero = motif_expected(1);
    expUni = sum(motif_expected([2, 4, 5, 6, 10, 11]));
    expMixed = sum(motif_expected([7, 8, 12, 13, 14, 15, 16]));
    expBi = sum(motif_expected([3, 9, 16]));


    % Empirical Values for 0 connections, purely uni connections, mixed
    % connections, and purely bi connections
    empZero = motif_empirical(1);
    empUni = sum(motif_empirical([2, 4, 5, 6, 10, 11]));
    empMixed = sum(motif_empirical([7, 8, 12, 13, 14, 15, 16]));
    empBi = sum(motif_empirical([3, 9, 16]));

    class_expected = [expZero, expUni, expMixed, expBi];
    class_empirical = [empZero, empUni, empMixed, empBi];
end
