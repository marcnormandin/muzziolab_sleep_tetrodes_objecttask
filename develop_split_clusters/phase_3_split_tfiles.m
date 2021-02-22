close all
clear all
clc

% This code splits the tfiles that were written for the unsplit (4 trials)
% data to separate t-files for each trial

%unsplitParentFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls\Young\CMG-young_OTNM_final';

expFiles = dir(fullfile('T:\projects\object_task_2021\recordings\nonmoved_controls', '**', 'experiment_description.json'));

for i = 1:length(expFiles)
    unsplitParentFolder = expFiles(i).folder;
    process_mouse(unsplitParentFolder);
end

function process_mouse(unsplitParentFolder)
    splitParentFolder = strrep(unsplitParentFolder, 'nonmoved_controls', 'nonmoved_controls_renamed');

    trialNames = {'s1', 's2', 's3', 's4'}; % Dont process s5 because that is the test session
    numTrials = length(trialNames);

    % For each trial read in the start and stop timestamps
    trialBounds_mus = zeros(numTrials, 2); % start, stop
    for iTrial = 1:numTrials
       trialFolder = fullfile(splitParentFolder, trialNames{iTrial});
       nvtFilename = fullfile(trialFolder, 'VT1.nvt');
       if ~isfile(nvtFilename)
           error('Expected nvt file does not exist (%s).', nvtFilename);
       end

       [timestamps_mus, ~, ~, ~, ~, ~, ~] = Nlx2MatVT(nvtFilename, [1 1 1 1 1 1], 1, 1, [] );

       trialBounds_mus(iTrial,1) = timestamps_mus(1);
       trialBounds_mus(iTrial,2) = timestamps_mus(end); % Hopefully this works even though another 512 values are after the last timestamp
    end
    
    % Read the original timestamps
    unsplitNvtFilename = fullfile(unsplitParentFolder, 'hab', 'VT1.nvt');
    [unsplitNlxNvtTimeStamps_mus, ~, ~, ~, ~, ~, ~] = Nlx2MatVT(unsplitNvtFilename, [1 1 1 1 1 1], 1, 1, [] );

    % Tfiles
    unsplitTFileFolder = fullfile(unsplitParentFolder, 'hab');
    tFileList = dir(fullfile(unsplitTFileFolder, 'TT*.t'));
    numTFiles = length(tFileList);
    
    if numTFiles < 1
        warning('There are no tfiles to process for (%s). Skipping.\n', unsplitParentFolder);
        return
    end

    % For each t-file
    for iT = 1:numTFiles

        tFilename = tFileList(iT).name;
        tFilenameFull = fullfile(unsplitTFileFolder, tFilename);

        % Get the original header so that we can replicate it.
        [header] = ml_nlx_mclust_load_spikes_header(tFilenameFull);

        % Read the original time timestamps
        %unsplitSpikeTimes_mus = ml_nlx_mclust_load_spikes_64bit(tFilenameFull) .* 10^6; % Convert to mus
        
        [unsplitSpikeTimes_mus, is32BitValid, is64BitValid] = read_spikes_mus(unsplitNlxNvtTimeStamps_mus, tFilenameFull);

        splitSpikeTimes_mus = cell(numTrials,1);
        for iTrial = 1:numTrials
           inds = find( unsplitSpikeTimes_mus >= trialBounds_mus(iTrial,1) & unsplitSpikeTimes_mus <= trialBounds_mus(iTrial,2) );
           splitSpikeTimes_mus{iTrial} = unsplitSpikeTimes_mus(inds);
        end

        % Check that the number of total spikes will remain the same
        numSpikesBeforeSplit = length(unsplitSpikeTimes_mus);
        numSpikesAfterSplit = 0;
        for iTrial = 1:numTrials
            numSpikesAfterSplit = numSpikesAfterSplit + length(splitSpikeTimes_mus{iTrial});
        end
        if numSpikesBeforeSplit ~= numSpikesAfterSplit
            warning('Spikes dont match for (%s). Num spikes before split (%d). Num spikes after split (%d).', tFilenameFull, numSpikesBeforeSplit, numSpikesAfterSplit);
        end

        % Save the respective spike times for each trial
        if ~is32BitValid && ~is64BitValid
            warning('The original file is not valid for 32 or 64 bits (%s). Skipping.\n', tFilenameFull);
        else
            for iTrial = 1:numTrials
                trialFolder = fullfile(splitParentFolder, trialNames{iTrial});
                outputTFullFilename = fullfile(trialFolder, tFilename);

                % Open a file to write the trial spikes times
                if isfile(outputTFullFilename)
                    delete(outputTFullFilename);
                end
                fp = fopen(outputTFullFilename, 'wb', 'b');
                if (fp == -1)
                 error(['Could not open file"' fn '".']);
                end
                %fprintf('Saving data to %s ... ', outputTFullFilename);
                for iH = 1:length(header)
                   fwrite(fp, sprintf('%s\n', header{iH}));
                end

                % Save as 64 bit regardless
                if is32BitValid || is64BitValid
%                     fwrite(fp, splitSpikeTimes_mus{iTrial} ./ 100, 'uint32');
%                 elseif is64BitValid
                    fwrite(fp, splitSpikeTimes_mus{iTrial} ./ 100, 'uint64');
                else
                    fprintf('Spikes are not 32 or 64 bit. Skipping.\n');
                end
                fclose(fp);
                %fprintf('done!\n');
            end % iTrial
        end % if 32 or 64 bit
    end % iT

end % function 


function [spikeTimes_mus, is32BitValid, is64BitValid] = read_spikes_mus(nlxNvtTimeStamps_mus, tFilename)
    % Load it as 32 bit
    spikeTimes_mclust = ml_nlx_mclust_load_spikes_32bit(tFilename);
    spikeTimes_mus_32bit = spikeTimes_mclust .* 10^6;
    is32BitValid = ml_nlx_mclust_spiketimes_are_valid(nlxNvtTimeStamps_mus, spikeTimes_mus_32bit, false);

    % Load it as 64 bit
    spikeTimes_mclust = ml_nlx_mclust_load_spikes_64bit(tFilename);
    spikeTimes_mus_64bit = spikeTimes_mclust .* 10^6;
    is64BitValid = ml_nlx_mclust_spiketimes_are_valid(nlxNvtTimeStamps_mus, spikeTimes_mus_64bit, false);

    if is32BitValid && ~is64BitValid
        spikeTimes_mus = spikeTimes_mus_32bit;
    elseif ~is32BitValid && is64BitValid
        spikeTimes_mus = spikeTimes_mus_64bit;
    elseif ~is32BitValid && ~is64BitValid
        spikeTimes_mus = [];
        warning('Spikes times are invalid!');
    end
end % function