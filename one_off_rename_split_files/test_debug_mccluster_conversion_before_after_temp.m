iTrial = 4
fn = 'myPoints';
% Map all points to respective trial
                tids = clusterToTrialId(clusterData.(fn));
                % Find the indices that belong to the current trial
                ids = find(tids == iTrial);
                %ids(11:end) = [];
                ids
                
                % Get the points that belong to the current trial
                tmp = clusterData.(fn)(ids)
                % Convert from unsplit values to split values
                trialPoints = unsplitToSplitIndex(tmp)
                