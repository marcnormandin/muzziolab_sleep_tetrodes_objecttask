close all

M = mltetrodeplacemap.meanFiringRateMapSmoothed;
O = mltetrodeplacemap.dwellTimeMapTrue;
O(O > 0) = 1;
O(O < 1) = 0;
O = 1 - O;
%O(end+1,:) = nan;
%O(:, end+1) = nan;

W = ones(size(O,1), size(O,2), 3);

figure
subplot(2,2,1)
mltetrodeplacemap.plot()
subplot(2,2,2)
%pcolor(flipud(O))
imagesc(O)
axis equal
subplot(2,2,3)
pcolor(M);
shading flat
hold on
h = imagesc(W)
hold off
set(h, 'AlphaData', O)
set(gca, 'ydir', 'reverse')
axis equal off

subplot(2,2,4)
mltetrodeplacemap.plot_path_with_spikes()
