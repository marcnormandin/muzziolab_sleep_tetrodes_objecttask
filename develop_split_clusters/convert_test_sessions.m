% This script was written as one-off to convert only the test session of
% the object task clusters into muzzis so that the custom cluster quality
% code could be run on the clusters.
% 2020-02-09

close all
clear all
clc

expFiles = dir(fullfile('T:\projects\object_task_2021\recordings\nonmoved_controls_renamed', '**', 'experiment_description.json'));

for i = 1:length(expFiles)
    unsplitParentFolder = expFiles(i).folder;
    process_mouse(unsplitParentFolder);
end

function process_mouse(unsplitParentFolder)
    testFolder = fullfile(unsplitParentFolder, 's5');
    
    % Get all of the regular cluster files. Do not keep autosave.clusters or
    % anything else that is not like TT#.clusters.
    clusterFiles = dir(fullfile(testFolder, '*.clusters'));
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

    for iFile = 1:numClusterFiles
        mclustClustersFilename = fullfile(clusterFiles(iFile).folder, clusterFiles(iFile).name);
        outputMuzzisFilename = strrep(mclustClustersFilename, '.clusters', '.muzzis');

        ml_nlx_mclust_convert_clusters_to_muzzis(mclustClustersFilename, outputMuzzisFilename)
    end
end

