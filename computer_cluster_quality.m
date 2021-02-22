close all
clear all
clc

startingFolder = pwd;

parentFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls_renamed';
outputFolder = strrep(parentFolder, 'recordings', 'analysis');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%%
% Search all subfolders of the parent folder for tfiles, and then CD into
% each folder in turn and run the cluster quality on it.
tfiles = dir(fullfile(parentFolder, '**', '*.t'));

trialFolders = unique({tfiles.folder});
problemFolders = {};
for iFolder = 1:length(trialFolders)
    fprintf('Processing %d of %d\n', iFolder, length(trialFolders));
    cd(trialFolders{iFolder})
    
    % Remove any pre-existing cluster quality files
    prevFiles = dir(fullfile(trialFolders{iFolder}, 'ML*-ClusterQual.mat'));
    for iFile = 1:length(prevFiles)
        delete(fullfile(prevFiles(iFile).folder, prevFiles(iFile).name));
    end
    
    % Remove any pre-existing cluster "wv" files
    prevFiles = dir(fullfile(trialFolders{iFolder}, '*-wv.mat'));
    for iFile = 1:length(prevFiles)
        delete(fullfile(prevFiles(iFile).folder, prevFiles(iFile).name));
    end
    
    try
        Create_CQ_File
    catch e
        problemFolders{end+1} = trialFolders{iFolder};
    end
end

%%
cqFiles = dir(fullfile(parentFolder, '**', 'ML*-ClusterQual.mat'));
results = [];
for i = 1:length(cqFiles)
    path = cqFiles(i).folder;
    name = cqFiles(i).name;

    s = split(path, filesep);
    session = s{end};
    mouse = s{end-1};
    group = s{end-2};
    experiment = s{end-3};
    
    s = split(name, '-');
    tetrode = [s{2} '_' s{3}];
    
    cq = load(fullfile(path, name));
    x = cq.CluSep;
    x.experiment = experiment;
    x.group = group;
    x.mouse = mouse;
    x.session = session;
    x.tetrode = tetrode;
    
    if isempty(results)
        results = x;
    else
        results(end+1) = x;
    end
end

%%

cd(startingFolder);

T = struct2table(results);

writetable(T, fullfile(outputFolder, 'cluster_analysis_NEW.xlsx'));

%%
%writetable(T, fullfile('Y:\\marc', '64bit_20200205.xlsx'));
