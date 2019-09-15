dbstop if error
clear; 
close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to aggregate majoy hyperspectral vegetation indices for every
% image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc
dbstop if error

dataPath  = 'T:\Box2\Drone Flight Data and Reference Files\Flight Data - All Sites\CLMB STND 2019 Flight Data\100081_2019_06_11_17_57_06\';

hdrPath   = strrep(dataPath, 'T:\Box2\Drone Flight Data and Reference Files\Flight Data - All Sites', 'T:\Results\AnalysisDroneData\ReflectanceCube\ReadableHDR');
hyperPath = strrep(dataPath, 'T:\Box2\Drone Flight Data and Reference Files\Flight Data - All Sites', 'T:\Results\AnalysisDroneData\ReflectanceCube\MATdataCube');
hyperPath = [hyperPath '\56\'];
gtPath    = 'T:\Results\AnalysisDroneData\grounTruth\CLMB STND 2019 Flight Data\100081_2019_06_11_17_57_06\gt_processed\';
INDICESpath = strrep(hyperPath, 'MATdataCube', 'indices');


list = dir([gtPath, '*.mat']);

% get the correct order of the files
fileIdx = [];
for ii = 1:length(list)
    tempFile = list(ii).name;
    fileIdx  = [fileIdx str2double(tempFile(isstrprop(tempFile, 'digit')))];
end
[~, idx] = sort(fileIdx);
list = list(idx);

% list_indices = {'aci', 'ari', 'cari', 'ci_red_edge', 'cri', 'evi', 'mari', ...
%     'mcari', 'mtci', 'ndvi', 'pri', 'psnd', 'rgri', 'rvsi', 'sipi', 'sr', 'vari', ...
%     'vi_green', 'wbi'};

list_indices = {'aci', 'ari', 'cari', 'ci_red_edge', 'evi', 'mari', ...
    'mcari', 'mtci', 'ndvi', 'pri', 'rgri', 'rvsi', 'sipi', 'sr', 'vari', ...
    'vi_green', 'wbi'};

result_overall_sw  = zeros(length(list), length(list_indices));
result_overall_nsw = zeros(length(list), length(list_indices));
result_overall     = {};
pixel_count        = zeros(length(list), 6);
for iIDX = 1:length(list_indices)
    result_overall{iIDX} = zeros(length(list), 6);
end

for iFile = 1:length(list)
    gtName          = list(iFile).name;
    load(fullfile(gtPath, gtName)) %gt
    gt = gt_final;
    gt_map = zeros(size(gt));
    gt_map(find(gt>0)) = 1;
    for iClass = 1:6
        x = find(gt == iClass);
        if ~isempty(x)
            pixel_count(iFile, iClass) = length(x);
        end
    end
    
    for iIDX = 1:length(list_indices)
        result_overall{iIDX} = zeros(length(list), 6);
        idxName   = strrep(gtName, 'ground_truth', 'raw');
        idxName   = strrep(idxName, '.mat', ['_', list_indices{iIDX}, '.mat']);
        loaded    = load(fullfile(INDICESpath, idxName));
        data      = getfield(loaded, list_indices{iIDX});
        sw        = data.*gt_map;
        
        for iClass = 1:6
            x = find(gt == iClass);
            if ~isempty(x)
                result_overall{iIDX}(iFile, iClass) = mean(sw(x));
            end
        end
        
        temp      = find(sw == 0);
        sw(temp)  = [];
        cond = isinf(sw) + isnan(sw);
        sw(cond == 1) = [];
        
        nsw       = data.*(1-gt_map);
        temp      = find(nsw == 0);
        nsw(temp) = [];
        cond = isinf(nsw) + isnan(nsw);
        nsw(cond == 1) = [];
        
        result_overall_sw(iFile, iIDX) = mean(sw);
        result_overall_nsw(iFile, iIDX) = mean(nsw);
    end   
end
