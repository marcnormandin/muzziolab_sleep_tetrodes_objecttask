close all
clear all
clc

recordingsParentDirectory = 'T:\projects\object_task_2021\recordings\nonmoved_controls_renamed';
analysisParentDirectory = strrep(recordingsParentDirectory, 'recordings', 'analysis');
if ~exist(analysisParentDirectory, 'dir')
    mkdir(analysisParentDirectory);
end

% Find all muzzis
muzzis = dir(fullfile(recordingsParentDirectory, '**', '*.muzzis'));

T = [];

for iMuzzi = 1:length(muzzis)
    fprintf('Processing muzzis (%d) of (%d)\n', iMuzzi, length(muzzis));
    
    tetrodeDataFolder = muzzis(iMuzzi).folder;
    [filepath, name, ext] = fileparts(muzzis(iMuzzi).name);
    tetrodeName = name;
    nttFilename = fullfile( tetrodeDataFolder, sprintf('%s.ntt', tetrodeName) );
    clusterFilename = fullfile( tetrodeDataFolder, sprintf('%s.muzzis', tetrodeName) );
    
    [clusterQuality] = compute_cluster_quality_all(nttFilename, clusterFilename);
    
    s = split(tetrodeDataFolder, filesep);
    sessionName = s{end};
    mouseName = s{end-1};
    groupName = s{end-2};
    experimentName = s{end-3};
    
    numClusters = length(clusterQuality);
    for iCluster = 1:numClusters
       clusterName = sprintf('%s_%d', tetrodeName, iCluster);
       
       clusterQuality(iCluster).experiment = experimentName;
       clusterQuality(iCluster).group = groupName;
       clusterQuality(iCluster).mouse = mouseName;
       clusterQuality(iCluster).session = sessionName;
       clusterQuality(iCluster).tetrode = tetrodeName;
       clusterQuality(iCluster).cluster = clusterName;
    end
    
    fprintf('Appending to table\n');
    T = [T; struct2table(clusterQuality)];
end % iMuzzi

% tetrodeDataFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls_renamed\Young\CMG-young_OTNM_final\s1';
% tetrodeName = 'TT3';

%% Save the results
outputFolder = analysisParentDirectory;
outputFilename = fullfile(outputFolder, 'cluster_analysis_muzzio_lab_nonmoved_20200208.xlsx');

if isfile(outputFilename)
    delete(outputFilename);
end

writetable(T, outputFilename)

%%






function [clusterQuality] = compute_cluster_quality_all(nttFilename, clusterFilename)
    [filepath1, name1, ext1] = fileparts(nttFilename);
    [filepath2, name2, ext2] = fileparts(clusterFilename);
    
    if strcmp(name1, name2) ~= 1
        error('%s does not match %s', name1, name2);
    end
    
    tetrodeName = name1; % or name2

    % Load the spike waveform data
    [Timestamps_mus, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nttFilename, [1 1 1 1 1], 1, 1, [] );

    clusterData = load( clusterFilename, '-mat' )

    numClusters = length(clusterData.Muzzio_Clusters)
    clusterInfo = struct('clusterName', [], 'spikeIndices', [], 'numSpikes', []);

    for iCluster = 1:numClusters
        clusterInfo(iCluster).clusterName = clusterData.Muzzio_Clusters{iCluster}.name;
        clusterInfo(iCluster).spikeIndices = clusterData.Muzzio_Clusters{iCluster}.myPoints;
        clusterInfo(iCluster).numSpikes = length(clusterInfo(iCluster).spikeIndices);
    end

    % Find shorted wires
    iChannelShorted = find(sum(sum(abs(Samples),3),1) < 1);
    if ~isempty(iChannelShorted)
        warning('%s channel %d is shorted', tetrodeName, iChannelShorted)
    end

    numOrigWires = size(Samples,2);
    wireIds = 1:numOrigWires;
    
    % Remove the shorted data since it is zero
    if ~isempty(iChannelShorted)
        Samples(:,iChannelShorted,:) = [];
        wireIds = setdiff(wireIds, iChannelShorted);
    end

    numSpikeSamples = size(Samples, 1);
    numWires = size(Samples,2);
    numSpikes = size(Samples,3);

    % Compute the feature matrix for ALL spikes (not just those of a cluster)
    featureMatrix = compute_feature_matrix(Samples);

    clusterQuality = struct('tetrodeName', [], 'clusterName', [], 'clusterNumber', [], 'isolationDistance', [], 'lRatio', []);
    for iCluster = 1:numClusters
        ci = clusterInfo(iCluster);

        fprintf('processing %s cluster (%d: %s) ... ', tetrodeName, iCluster, ci.clusterName);
        [cq] = compute_cluster_quality(featureMatrix, ci);
        wire_snr = compute_snr(Samples, ci);
        fprintf('done!\n');

        clusterQuality(iCluster).tetrodeName = tetrodeName;
        clusterQuality(iCluster).clusterName = ci.clusterName;
        clusterQuality(iCluster).clusterNumber = iCluster;
        clusterQuality(iCluster).isolationDistance = cq.isolationDistance;
        clusterQuality(iCluster).lRatio = cq.lRatio;
        clusterQuality(iCluster).numSpikes = ci.numSpikes;
        %clusterQuality(iCluster).wireIds = wireIds;
        
        % default of 0
        for iWire = 1:numOrigWires
            clusterQuality(iCluster).(sprintf('SNR_%d', iWire)) = 0;
        end
        
        for iWire = 1:length(wire_snr)
            clusterQuality(iCluster).(sprintf('SNR_%d', wireIds(iWire))) = wire_snr(iWire);
        end
    end
end % function



    
%%


% iCluster = 1;
% 
% ci = clusterInfo(iCluster);
% 
% wv = Samples(:,:, ci.spikeIndices);
% plot_spike_waveforms(ci, wv)

%%
function [clusterSnr] = compute_snr(Samples, ci)
    numWires = size(Samples,2);
    %% Compute the SNR (not exact match...)
    wv = Samples(:,:, ci.spikeIndices);
    clusterAvg = mean(wv, 3);
    clusterStd = std(wv, 0, 3);

    noiseAvg = mean(mean(Samples, 3), 1);
    noiseStd = mean(std(Samples, 0, 3), 1);

    clusterSNR = zeros(1, numWires);

    magicNumber = sqrt(2);

    for iWire = 1:numWires
        % For each wire, find the maximum absolute value of the average
        % waveform
        [~, clusterPeakIndex] = max( abs(clusterAvg(:,iWire)) );
        peakValue = clusterAvg(clusterPeakIndex, iWire);

        %x = abs(peakValue) - sign(peakValue) .* noiseAvg(clusterPeakIndex,iWire);

        clusterSnr(iWire) = abs( peakValue - noiseAvg(iWire) ) ./ noiseStd(iWire) * magicNumber;
    end % iWire
end % function

%%
function [clusterQuality] = compute_cluster_quality(featureMatrix, ci)
    clusterSpikeIndices = ci.spikeIndices;
    clusterSpikeIndices(clusterSpikeIndices > size(featureMatrix,1)) = [];
    
    %numSpikes = size(featureMatrix,1);
    df = size(featureMatrix,2); % number of features
    
    if length(clusterSpikeIndices) <= df
        warning('Number of cluster spikes (%d) is below minimum of (%d). Cannot compute mahalanobis distance. Returning nan.', length(clusterSpikeIndices), df);
        clusterQuality.lRatio = nan;
        clusterQuality.isolationDistance = nan;
        return
    end
    
    m = mahal(featureMatrix,featureMatrix(clusterSpikeIndices,:));
    
    % Remove the few nan values that might exist
    badi = isnan(m);
    if ~isempty(badi)
        validIndices = setdiff(1:length(m), badi);
        m(badi) = [];

        % Now we have to map each cluster spike index to its location in
        % the new list
        newClusterSpikeIndices = [];
        for iClusterSpike = 1:length(clusterSpikeIndices)
           ind = find(validIndices == clusterSpikeIndices(iClusterSpike));
           if ~isempty(ind)
               k = length(newClusterSpikeIndices) + 1;
               newClusterSpikeIndices(k) = ind;
           end
        end
        
        clusterSpikeIndices = newClusterSpikeIndices;
    end
    
    numSpikes = length(m);
    nClusterSpikes = length(clusterSpikeIndices);
    nNotClusterSpikes = numSpikes - nClusterSpikes;
    
    %mClusterSpikes = m(ci.spikeIndices);
    mNotClusterSpikes = m(setdiff(1:numSpikes, clusterSpikeIndices));
    

    % For cluster C, containing nc spikes, the isolation distance is defined as the squared Mahalanobis distance of the nc-th closest non-c spike to the center of C
    % Larger values indicate greater isolation
    % Not normalized. Clusters with lots of spikes will have a higher isolation
    % distance.
    if nNotClusterSpikes >= nClusterSpikes && nClusterSpikes <= numSpikes / 2
        tmp = sort(mNotClusterSpikes); % sort ascending
        isolationDistance = tmp(nClusterSpikes);
    else
        isolationDistance = nan;
    end

    pNotClusterSpikes = 1 - chi2cdf(mNotClusterSpikes, df);
    
    if any(isnan(pNotClusterSpikes))
        error('found a nan');
    end
    
    %pClusterSpikes = 1 - chi2cdf(mClusterSpikes, df);

    %pMeanNotClusterSpikes = mean(pNotClusterSpikes);
    %pMeanClusterSpikes = mean(pClusterSpikes);

    % Low L-Ratio means lots of empty space between cluster and other spikes
    lRatio = sum( pNotClusterSpikes ) / nClusterSpikes;
    
    if isnan(lRatio)
        error('Found a nan in L-ratio!');
    end
    


    
    clusterQuality.isolationDistance = isolationDistance;
    clusterQuality.lRatio = lRatio;
end % function

function [spikeFeatures] = compute_feature_matrix(Samples)
    numSpikeSamples = size(Samples, 1);
    numWires = size(Samples,2);
    numSpikes = size(Samples,3);

    numFeatures = 2; % energy and PC1
    spikeFeatures = nan(numSpikes, numFeatures * numWires);

    % Calculate energy of each spikes on each wire
    for iWire = 1:numWires
        x = squeeze(Samples(:,iWire,:));

        % Compute energy of each spikes per wire and normalize by the number of samples
        spikeFeatures(:, iWire) = sqrt(sum(x .* x, 1)) ./ numSpikeSamples;
    end

    % Calculate the first principal component of the energy-normalized spikes
    for iWire = 1:numWires
        x = squeeze(Samples(:,iWire,:));
        e = sqrt(sum(x .* x, 1));
        E = repmat(e, numSpikeSamples, 1);
        x = x ./ E; % normalize each spike component by the corresponding spike energy

        [pc,wpc] = pca(x');

        % Get the first principal component coefficients for each spike
        spikeFeatures(:,iWire + numWires) = wpc(:,1);
    end

end % function


function plot_spike_waveforms(clusterInfo, wv)
    
    numWires = size(wv, 2); % 4
    numSpikes = size(wv,3);
    spikeLength = size(wv, 1);
    
    figure('name', clusterInfo.clusterName)
    for iWire = 1:numWires
        ax(iWire) = subplot(2,numWires, iWire);
        for iSpike = 1:numSpikes
           plot( 1:spikeLength, wv(:, iWire, iSpike) );
           hold on
        end
    end % iWire
    
    wvAvg = mean(wv, 3);
    for iWire = 1:numWires
        ax(iWire + numWires) = subplot(2,numWires, iWire+numWires);
        plot( 1:spikeLength, wvAvg(:, iWire), 'r-', 'linewidth', 4 );
    end % iWire
    
    linkaxes(ax, 'xy')
    
end % function
