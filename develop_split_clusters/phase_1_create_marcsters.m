close all
clear all
clc

expFiles = dir(fullfile('T:\projects\object_task_2021\recordings\nonmoved_controls', '**', 'experiment_description.json'));

for i = 1:length(expFiles)
    %unsplitParentFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls\Young\CMG-young_OTNM_final';
    unsplitParentFolder = expFiles(i).folder;
    process_mouse(unsplitParentFolder);
end

function process_mouse(unsplitParentFolder)
fprintf('Processing (%s)\n', unsplitParentFolder);

splitParentFolder = strrep(unsplitParentFolder, 'nonmoved_controls', 'nonmoved_controls_renamed');

% Get all of the regular cluster files. Do not keep autosave.clusters or
% anything else that is not like TT#.clusters.
clusterFiles = dir(fullfile(unsplitParentFolder, 'hab', '*.clusters'));
keep = zeros(1, length(clusterFiles));
for i = 1:length(clusterFiles)
    if regexp(clusterFiles(i).name, '^(TT)\d(.clusters)$')
        keep(i) = 1;
    else
        keep(i) = 0;
    end
end
clusterFiles(~keep) = [];
numClusterFiles = length(clusterFiles);

%%
for iClusterFile = 1:numClusterFiles
    clusterFilename = fullfile(clusterFiles(iClusterFile).folder, clusterFiles(iClusterFile).name);
    tetrodeName = sprintf('TT%s', clusterFiles(iClusterFile).name(3));
    
    nlxUnsplitFilename = fullfile(unsplitParentFolder, 'hab', sprintf('%s.ntt', tetrodeName));
    
    sessionNames = {'s1', 's2', 's3', 's4'};
    numSessions = length(sessionNames);
    numSpikesSplit = 0;
    trialLengths = zeros(1, numSessions);
    trialBounds = zeros(numSessions, 2);
    for iSession = 1:numSessions
        nlxFilenameSplit = fullfile(splitParentFolder, sessionNames{iSession}, sprintf('%s.ntt', tetrodeName));
        %fprintf('Processing %s\n', nlxFilenameSplit);
        
        % Load the associated NTT file
        [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nlxFilenameSplit, [1 1 1 1 1], 1, 1, [] );
        
        numSpikesSplit = numSpikesSplit + length(Timestamps);
        
        trialLengths(iSession) = length(Timestamps);
        % We need to keep the trial indices start and end so we can split
        % the cluster data later
        if iSession == 1
            trialBounds(iSession,1) = 1;
        else
            trialBounds(iSession,1) = sum(trialLengths(1:iSession-1)) + 1;
        end
        trialBounds(iSession,2) = sum(trialLengths(1:iSession));
    end
    
    % Load the associated NTT file
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nlxUnsplitFilename, [1 1 1 1 1], 1, 1, [] );
    
    % Load the unsplit cluster data
    [numClusters, clusterData, clusterColours, allPoints] =  load_cluster_data(clusterFilename);
    
    numSpikesUnsplit = length(Timestamps);
    maxUnsplitIndex = max(allPoints);
    
    
    if numSpikesUnsplit ~= numSpikesSplit || maxUnsplitIndex > numSpikesSplit
        fprintf('Total spikes before split (%d)\n', numSpikesUnsplit);
        fprintf('Total spikes after split (%d)\n', numSpikes);
        fprintf('Max index in clusters before split (%d)\n', maxUnsplitIndex);
    end
    
    %%
    clusterToTrialId = zeros(1, trialBounds(end,2));
    for iTrial = 1:size(trialBounds,1)
        clusterToTrialId(trialBounds(iTrial,1):trialBounds(iTrial,2)) = iTrial;
    end
    
    unsplitToSplitIndex = zeros(1, trialBounds(end,2));
    for iTrial = 1:size(trialBounds,1)
        unsplitToSplitIndex(trialBounds(iTrial,1):trialBounds(iTrial,2)) = 1:trialLengths(iTrial);
    end
    
    %%
    % Load the data in the unsplit cluster file
    data = load(clusterFilename, '-mat');
    numClusters = length(data.MClust_Clusters);
    clusterColours = data.MClust_Colors;
    
    MCC = {};
    trialClusterData = cell(numClusters, numSessions);
    for iCluster=1:numClusters
        clusterData = data.MClust_Clusters{iCluster};
        
        % Convert the unsplit indices into local trial indices
        %clusterData.myPoints = unsplitToSplitIndex(clusterData.myPoints);
        %clusterData.myOrigPoints = unsplitToSplitIndex(clusterData.myOrigPoints);
        %clusterData.ForbiddenPoints = unsplitToSplitIndex(clusterData.ForbiddenPoints);
        
        fieldNames = {'myPoints', 'myOrigPoints', 'ForbiddenPoints'};
        
        % Each cluster is split into 4 sessions
        
        for iTrial = 1:numSessions
            trialClusterData{iCluster, iTrial} = clusterData;
            for iField = 1:length(fieldNames)
                fn = fieldNames{iField};
                % Map all points to respective trial
                tids = clusterToTrialId(clusterData.(fn));
                % Find the indices that belong to the current trial
                ids = find(tids == iTrial);
                % Get the points that belong to the current trial
                trialClusterData{iCluster, iTrial}.(fn) = clusterData.(fn)(ids);
                % Convert from unsplit values to split values
                trialClusterData{iCluster, iTrial}.(fn) = unsplitToSplitIndex(trialClusterData{iCluster, iTrial}.(fn));
                
                % Reshape to same idea as MClust (#points x 1)
                trialClusterData{iCluster, iTrial}.(fn) = reshape(trialClusterData{iCluster, iTrial}.(fn), length(trialClusterData{iCluster, iTrial}.(fn)), 1);
                
                % THIS WAS A HACK
                %trialClusterData{iCluster, iTrial}.(fn) = trialClusterData{iCluster, iTrial}.(fn) - 1;
            end
        end
    end
    
    outputFolder = fullfile(unsplitParentFolder, 'marc_clusters');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
     
    newFilename = fullfile(outputFolder, sprintf('%s.marcsters', tetrodeName));
    fprintf('Saving file: %s\n', newFilename);
    if isfile(newFilename)
        delete(newFilename);
    end
    
    save(newFilename, 'trialClusterData', 'clusterColours', 'numClusters', 'numSessions', 'tetrodeName');
end % iCluster
end % function

%%

%%



%%
function [numClusters, clusterData, clusterColours, allPoints] =  load_cluster_data(clusterFilename)

data = load(clusterFilename, '-mat');
numClusters = length(data.MClust_Clusters);
clusterColours = data.MClust_Colors;

allPoints = [];
for iCluster=1:numClusters
    clusterData{iCluster} = data.MClust_Clusters{iCluster};
    
    allPoints = [allPoints; clusterData{iCluster}.myOrigPoints];
end
allPoints = sort(allPoints);
end % function