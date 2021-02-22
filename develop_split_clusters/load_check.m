close all
clear all
clc

unsplitParentFolder = 'T:\projects\object_task_2021\recordings\nonmoved_controls\Young\CMG-young_OTNM_final';
splitParentFolder = strrep(unsplitParentFolder, 'nonmoved_controls', 'nonmoved_controls_renamed');

nlxFilenameUnsplit = fullfile(unsplitParentFolder, 'hab', 'TT3.ntt');
nlxFilenameSplit = fullfile(splitParentFolder, 's1', 'TT3.ntt');


% Load the split data
[Timestamps, ScNumbers, CellNumbers, Features1, Samples1, Header1] = Nlx2MatSpike(nlxFilenameSplit, [1 1 1 1 1], 1, 1, [] );

