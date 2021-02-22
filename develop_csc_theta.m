close all
clear all
clc

recFolder = 'T:\projects\object_task_2021\recordings\main\Old_Controls\R53c\s1';
%recFolder = 'T:\projects\object_task_2021\recordings\main\Young_SD\R93sd\s1';
%recFolder = 'T:\projects\object_task_2021\recordings\main\Old_Controls\R76c\s1';

cscFilename = fullfile(recFolder, 'CSC2.ncs');
nvtFilename = fullfile(recFolder, 'VT1.nvt');
tFilename = fullfile(recFolder, 'TT1_1.t');
numBits = -1;

[nlxNvtTimeStamps_mus, X, Y, Angles, Targets, Points, Header] = ...
           Nlx2MatVT(nvtFilename, [1 1 1 1 1 1], 1, 1, [] );



[Timestamps_mus, ChannelNumbers, SampleFrequencies,NumberOfValidSamples, Samples, Header] = ...
    Nlx2MatCSC(cscFilename,[1 1 1 1 1], 1, 1, [] );

[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ...
           Nlx2MatSpike('test.nst', [1 1 1 1 1], 1, 1, [] );

ts = zeros(1, numel(Samples));

dt = median(diff(Timestamps_mus));

k = 1;
for i = 1:numel(Timestamps_mus)
   
   for j = 1:size(Samples,1) % 1:512
       ts(k) = Timestamps_mus(i) + dt / (size(Samples,1)+1) * j;
       k = k + 1;
   end
end


ts_ms = ts ./ 10^3;
ts_origin_ms = ts_ms(1);

ts_ms = ts_ms - ts_origin_ms; % zero


spikeTimes_mus = ml_nlx_load_mclust_spikes_as_mus(nlxNvtTimeStamps_mus, tFilename, numBits);

spikeTimes_ms = spikeTimes_mus ./ 10^3;
spikeTimes_ms = spikeTimes_ms - ts_origin_ms;







fs = 1 ./ (ts_ms(2) - ts_ms(1)) * 1000;

cscr = reshape(Samples, 1, numel(Samples));

figure
plot(ts_ms, cscr)
title('CSC Raw')


cscf = ml_ephys_remove_powerline_noise(cscr, fs);
figure
plot(ts_ms, cscf)
title('Filtered')


% detrend
W_secs = 0.1;
W_points = W_secs * fs;
cscz = cscf;
a = 1;
b = a + W_points;
while 1
    x = cscz(a:b);
    y = detrend(x,1);
    cscz(a:b) = y;
    
    a = a + W_points;
    b = b + W_points;
    if a > numel(cscz)
        break;
    end
    if b > numel(cscz)
        b = numel(cscz);
    end
end

% figure
% plot(ts_ms, csc, 'k-')
% hold on
% plot(ts_ms, cscz, 'r-')
% grid on

figure
plot(ts_ms, cscz, 'k-');
hold on
plot(spikeTimes_ms, zeros(1, numel(spikeTimes_ms)), 'ro', 'markerfacecolor', 'r', 'markersize', 10)
