function MCC = make_mccluster(clusterData)
fieldNames = fields(clusterData);

% Hack to split a cluster file
% MCC.name = '';
% MCC.xdimNames = {}; 
% MCC.ydimNames = {}; 
% MCC.xdimSources = {};
% MCC.ydimSources = {};
% MCC.cx = {}; 
% MCC.cy = {};
% MCC.AddFlag = [];
% MCC.recalc = -1;
% MCC.myPoints = [];
% MCC.myOrigPoints = [];
% MCC.ForbiddenPoints = [];

for iField = 1:length(fieldNames)
   fn = fieldNames{iField};
   MCC.(fn) = clusterData.(fn);
end

MCC = class(MCC, 'mccluster');