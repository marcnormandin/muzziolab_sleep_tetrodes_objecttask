% I wrote this script because I was given two separate datasets. The first
% has only the t files and nvt file used for the object task placemaps and
% correlations, and the second has the raw data needed for the cluster
% metrics. Since the metrics should be done on the same data, I wrote this
% script to copy the required data into the correct folders so the metrics
% are run on the actual clusters used. Everything should have been given to
% me as one dataset....

close all
clear all
clc

DONT RUN THIS AGAIN


fromParent = 'T:\ObjectTask_SleepStudy\datasets\Object Task';
toParent = 'T:\projects\object_task_2021\recordings\main';

Tfrom = make_data_table(fromParent);
Tto = make_data_table(toParent);

Tcommon = intersect(Tfrom, Tto)

for iRec = 1:size(Tcommon,1)
   rec = Tcommon(iRec,:);
   
   fromFolder = fullfile(fromParent, rec.group{1}, rec.animal{1}, rec.session{1});
   toFolder = fullfile(toParent, rec.group{1}, rec.animal{1}, rec.session{1});
   
   fileTypes = {'*.txt', '*.mclust', '*.clusters', '*.ncs', '*.nev', '*.ntt', '*.fd'};
   dirTypes = {'FD'};
   
   for iFileType = 1:length(fileTypes)
      ft = fileTypes{iFileType};
      
      files = dir( fullfile(fromFolder, ft) );
      if ~isempty(files)
          for iFile = 1:length(files)
              copyfile( fullfile(fromFolder, files(iFile).name), fullfile(toFolder, files(iFile).name) );
          end
      end
   end
   
   for iDir = 1:length(dirTypes)
      ft = dirTypes{iDir};
      
      if exist(fullfile(fromFolder, ft), 'dir')
        copyfile( fullfile(fromFolder, ft), fullfile(toFolder, ft) );
      end
   end
end

function T = make_data_table(parentFolder)
    matches = dir(fullfile(parentFolder, '**', '*.t'));
    folders = unique({matches.folder});
    t = [];
    for i = 1:length(folders)
        k = length(t) + 1;
        s = split(folders(i), filesep);
        t(k).group = s{end-2};
        t(k).animal = s{end-1};
        t(k).session = s{end};
    end
    T = struct2table(t);
end