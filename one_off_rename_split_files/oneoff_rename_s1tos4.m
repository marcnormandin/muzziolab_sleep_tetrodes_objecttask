parentFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls_split';

expFiles = dir(fullfile(parentFolder, '**', 'experiment_description.json'));

for iExp = 1:length(expFiles)
    sessionParentFolder = expFiles(iExp).folder;
    fprintf('Processing %s\n', sessionParentFolder);

    sessionNames = {'s1', 's2', 's3', 's4'}; % s5 is already single
    for iSession = 1:length(sessionNames)
        sessionFolder = fullfile(sessionParentFolder, sessionNames{iSession});
        fprintf('\t %s\n', sessionFolder);

        files = dir(sessionFolder);
        files([files.isdir] == 1) = [];

        % First load all of the filename information to make sure that it
        % all matches
        ext = {};
        main = {};
        t1 = {};
        t2 = {};
        for iFile = 1:length(files)
            originalName = files(iFile).name;
            fprintf('\t\t%s\n', originalName);

            tmp = split(originalName, '.');
            s1 = tmp{1};
            ext{iFile} = tmp{2};
            tmp = split(s1, '_');
            main{iFile} = tmp{1};
            t1{iFile} = tmp{2};
            t2{iFile} = tmp{3};
            %fprintf('%s %s %s <--> %s\n', main{iFile}, t1{iFile}, t2{iFile}, s2);
        end
        
        if length(unique(t1)) ~= 1 || length(unique(t2)) ~= 1
            disp(t1)
            disp(t2)
            error('Check folder: %s\n', sessionFolder);
        end

        for iFile = 1:length(files)
            originalName = files(iFile).name;
            newName = sprintf('%s.%s', main{iFile}, ext{iFile});

            inputFolder = sessionFolder;
            outputFolder = strrep(sessionFolder, 'nonmoved_controls_split', 'nonmoved_controls_renamed');
            if ~exist(outputFolder, 'dir')
                mkdir(outputFolder);
            end
            fprintf('%s\n', outputFolder);
            fprintf('%s -> %s\n', originalName, newName);

            copyfile(fullfile(inputFolder, originalName), fullfile(outputFolder, newName));
        end
    end
end