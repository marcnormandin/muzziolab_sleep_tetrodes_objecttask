close all
clear all
clc

%unsplitParentFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls\Young\CMG-young_OTNM_final';

expFiles = dir(fullfile('T:\projects\object_task_2021\recordings\nonmoved_controls', '**', 'experiment_description.json'));

for i = 1:length(expFiles)
    %unsplitParentFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls\Young\CMG-young_OTNM_final';
    unsplitParentFolder = expFiles(i).folder;
    process_mouse(unsplitParentFolder);
end

function process_mouse(unsplitParentFolder)

    marcstersFolder = fullfile(unsplitParentFolder, 'marc_clusters');

    splitParentFolder = strrep(unsplitParentFolder, 'nonmoved_controls', 'nonmoved_controls_renamed');

    sessionNames = {'s1', 's2', 's3', 's4'};

    % Get list of the marcsters
    marcsters = dir(fullfile(marcstersFolder, '*.marcsters'));
    numMarcsters = length(marcsters);

    for iMarcster = 1:numMarcsters
        mfn = fullfile(marcsters(iMarcster).folder, marcsters(iMarcster).name);

        marcster = load(mfn, '-mat');

        tetrodeName = marcster.tetrodeName;
        numClusters = marcster.numClusters;
        numSessions = marcster.numSessions;
        numTrials = numSessions;

        for iTrial = 1:numSessions
            %MCC = struct('MClust_Clusters', {}, 'MClust_Colors', []);
            MClust_Colors = marcster.clusterColours;
            MClust_Clusters = {};
            for iCluster = 1:numClusters
                MClust_Clusters{iCluster} = mccluster(marcster.trialClusterData{iCluster, iTrial});
            end

            outputFolder = fullfile(splitParentFolder, sessionNames{iTrial});
            ofn = fullfile(outputFolder, sprintf('%s.clusters', tetrodeName));

            if isfile(ofn)
                delete(ofn)
            end
            fprintf('Saving file (%s)\n', ofn);
            save(ofn, 'MClust_Clusters', 'MClust_Colors', '-mat');

            % Create the FD directory if it doesnt already exist
            fdDirectory = fullfile(outputFolder, 'FD');
            if ~exist(fdDirectory, 'dir')
                mkdir(fdDirectory);
            end
        end
    end
end % function
