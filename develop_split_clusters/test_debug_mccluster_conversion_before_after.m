close all
clear all
clc

origFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls\Young\CMG-young_OTNM_final\hab';
origData = load(fullfile(origFolder, 'TT2.clusters'), '-mat')


iCluster = 3;

fn = 'myPoints';
x1 = origData.MClust_Clusters{iCluster}.(fn);
nx = length(x1)

pts = [];
for iTrial = 1:4
    newFolder = fullfile('T:\projects\object_task_2021\recordings\nonmoved_controls_renamed\Young\CMG-young_OTNM_final', sprintf('s%d', iTrial));
    newData = load(fullfile(newFolder, 'TT2.clusters'), '-mat');
   z =  newData.MClust_Clusters{iCluster}.(fn);
   pts = [pts; z];
end
np = length(pts)
