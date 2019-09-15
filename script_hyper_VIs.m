dbstop if error
clear; 
close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to calculate majoy hyperspectral vegetation indices
% Reference: Hyperspectral Remote Sensing of Vegetation (Second Edition, Volumne II)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataPath  = 'T:\Box2\Drone Flight Data and Reference Files\Flight Data - All Sites\CLMB GWAS 2019 Flight Data\100083_2019_06_25_15_59_59\';
hdrPath   = strrep(dataPath, 'T:\Box2\Drone Flight Data and Reference Files\Flight Data - All Sites', 'T:\AnalysisDroneData\ReflectanceCube\ReadableHDR');
hyperPath = strrep(dataPath, 'T:\Box2\Drone Flight Data and Reference Files\Flight Data - All Sites', 'T:\AnalysisDroneData\ReflectanceCube\MATdataCube');
% hyperPath = [hyperPath '\56\'];

%%%% Get the list of all files 
list = dir([hyperPath, 'raw_*.mat']);
    % get the correct order of the files
fileIdx = [];
for ii = 1:length(list)
    tempFile = list(ii).name;
    fileIdx  = [fileIdx str2double(tempFile(isstrprop(tempFile, 'digit')))];
end
[~, idx] = sort(fileIdx);
list = list(idx);

% load flags of wavelengths
% load('flagGoodWvlen.mat')

% definition of indices result save path. If directory not exist, create one
indexPath = strrep(dataPath, 'T:\Box2\Drone Flight Data and Reference Files\Flight Data - All Sites', 'T:\AnalysisDroneData\ReflectanceCube\indices');
if ~exist(indexPath, 'dir')
    mkdir(indexPath)
end

% load wavelength infomation
load([hdrPath, list(1).name])% parameters, wavelength

%% 
% list of interested wavelengths which are used in calculating the vege indices
list_wv  = [445, 500, 510, 531, 550, 570, 650, 670, 675, 680, 681.25, 700, 708.75, 714, 733, 752, 753.75, 800, 900, 970];
idx = {};
for i = 1:length(list_wv)
    wv = list_wv(i);
%     [~, idx(i)] = min(abs(wavelength-wv));
    [~,temp]=sort(abs(wavelength-wv), 'ascend');
%     idx{i} = temp(1:1);
    idx{i} = [temp(1)-2 temp(1)-1 temp(1) temp(1)+1 temp(1)+2];
end

%  Starting and ending indices of wavelengths of specific colors
% green
[~, idx_green(1)] = min(abs(wavelength-540));
[~, idx_green(2)] = min(abs(wavelength-560));

% red-edge
[~, idx_red(1)] = min(abs(wavelength-690));
[~, idx_red(2)] = min(abs(wavelength-710));

% blue
[~, idx_blue(1)] = min(abs(wavelength-450));
[~, idx_blue(2)] = min(abs(wavelength-490));

% nir
[~, idx_nir(1)] = min(abs(wavelength-760));
[~, idx_nir(2)] = min(abs(wavelength-800));

h      = ones(3,3)/8;
h(2,2) = 0;
%% loop over all images
for iFile = 1:length(list)
    fileName  = list(iFile).name;
    fileName_ = strrep(fileName, '_rd_rf.mat', '');
    load([hyperPath, fileName]) % data
    R_nir    = mean(data(:,:,idx_nir(1):idx_nir(2)), 3);
    R_red    = mean(data(:,:,idx_red(1):idx_red(2)), 3);
    R_green  = mean(data(:,:,idx_green(1):idx_green(2)), 3);
    R_blue   = mean(data(:,:,idx_blue(1):idx_blue(2)), 3);

	R_wv   = [];
    R_wv_  = [];
    for i = 1:length(list_wv)
        R_wv(:,:,i)  = mean(data(:,:,idx{i}), 3);
%         R_wv_(:,:,i) = imfilter(R_wv(:,:,i), h);
    end
    
    % getting rid of 0's to avoid the problem of "deviding by 0" when calculating indices
    % use imfilter?
    
    % NIR
    index      = find(R_nir == 0);
    while ~isempty(index)
        R_nir_   = imfilter(R_nir,h);
        R_nir(index) = R_nir_(index);
        index        = find(R_nir == 0);
    end
    
    % RED
	index      = find(R_red == 0);    
    while ~isempty(index)
        R_red_   = imfilter(R_red,h);
        R_red(index) = R_red_(index);
        index        = find(R_red == 0);
    end

    % GREEN
	index      = find(R_green == 0);    
    while ~isempty(index)
        R_green_ = imfilter(R_green, h);
        R_green(index) = R_green_(index);
        index = find(R_green == 0);
    end
    
    % BLUE
	index      = find(R_blue == 0); 
    while ~isempty(index)
        R_blue_  = imfilter(R_blue, h);  
        R_blue(index) = R_blue_(index);
        index = find(R_blue == 0);
    end 
    
    % other specific wavelengths
    for i = 1:length(list_wv)
        temp = R_wv(:,:,i);
        index = find(temp == 0);        
        while ~isempty(index)
            temp_ = imfilter(temp, h);
            temp(index) = temp_(index);
            index = find(temp == 0);
        end
        R_wv(:,:,i) = temp;
    end
    
    % calculation of indices 
    % 1DL_DGVI
    
    % 1DZ_DGVI
    
    % ACI
    aci = R_green./R_nir;
    save(fullfile(indexPath, [fileName_,'_aci.mat']), 'aci')
    figure, imagesc(aci), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_aci.png']), 'png')
    
    % ARI
    ari = (1./R_wv(:,:,list_wv == 550))-(1./R_wv(:,:,list_wv == 700));
	save(fullfile(indexPath, [fileName_, '_ari.mat']), 'ari')
    figure, imagesc(ari), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_ari.png']), 'png')
    
    % ARVI
    
    % CAI
    
    % CARI
    cari = (R_wv(:,:,list_wv == 700)-R_wv(:,:,list_wv == 670)) - 0.2*(R_wv(:,:,list_wv == 700)-R_wv(:,:,list_wv == 550));
    save(fullfile(indexPath, [fileName_, '_cari.mat']), 'cari')
    figure, imagesc(cari), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_cari.png']), 'png')
    
    % CI_red_edge
    ci_red_edge = R_nir./R_red-1;
	save(fullfile(indexPath, [fileName_,'_ci_red_edge.mat']), 'ci_red_edge')
    figure, imagesc(ci_red_edge), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_ci_red_edge.png']), 'png')
    
    % CRI
    cri = [];
    cri(:,:,1) = (1./R_wv(:,:,list_wv == 510))-(1./R_wv(:,:,list_wv == 550));
    cri(:,:,2) = (1./R_wv(:,:,list_wv == 510))-(1./R_wv(:,:,list_wv == 700));
	save(fullfile(indexPath, [fileName_, '_cri.mat']), 'cri')
    figure, subplot(1,2,1), imagesc(cri(:,:,1)), axis image
    subplot(1,2,2), imagesc(cri(:,:,2)), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_cri.png']), 'png')
    
    % EVI
    evi = 2.5*(R_nir-R_red)./(R_nir+6*R_red-7.5*R_blue+1);
	save(fullfile(indexPath, [fileName_, '_evi.mat']), 'evi')
    figure, imagesc(evi), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_evi.png']), 'png')
    
    % MARI
    mari = (1./R_wv(:,:,list_wv == 550)-1./R_wv(:,:,list_wv == 700)).*R_nir;
    save(fullfile(indexPath, [fileName_, '_mari.mat']), 'mari')
    figure, imagesc(mari), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_mari.png']), 'png')
    
    % MCARI
    mcari = ((R_wv(:,:,list_wv == 700)-R_wv(:,:,list_wv == 670))-0.2*(R_wv(:,:,list_wv == 700)-R_wv(:,:,list_wv == 550))).*(R_wv(:,:,list_wv == 700)./R_wv(:,:,list_wv == 670));
    save(fullfile(indexPath, [fileName_, '_mcari.mat']), 'mcari')
    figure, imagesc(mcari), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_mcari.png']), 'png')
    
    % MSI (SWIR not available)

    % MTCI
    mtci = (R_wv(:,:,list_wv == 753.75)-R_wv(:,:,list_wv == 708.75))./(R_wv(:,:,list_wv == 708.75)-R_wv(:,:,list_wv == 681.25));
	save(fullfile(indexPath, [fileName_, '_mtci.mat']), 'mtci')
    figure, imagesc(mtci), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_mtci.png']), 'png')
    
    % NDII (SWIR not available)
    
    % NDLI (wavelengths 1754, 1680 not available)
    
    % NDNI (wavelengths 1510, 1680 not available)
    
    % NDVI
    ndvi = (R_nir-R_red)./(R_nir+R_red);
	save(fullfile(indexPath, [fileName_, '_ndvi.mat']), 'ndvi')
    figure, imagesc(ndvi), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_ndvi.png']), 'png')
    
    % NDWI (wavelength 1240 not available)
    
    % PRI
    pri = (R_wv(:,:,list_wv == 531)-R_wv(:,:,list_wv == 570))./(R_wv(:,:,list_wv == 531)+R_wv(:,:,list_wv == 570));
	save(fullfile(indexPath, [fileName_, '_pri.mat']), 'pri')
    figure, imagesc(pri), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_pri.png']), 'png')
    
    % PSND
    psnd = [];
    psnd(:,:,1) = (R_wv(:,:,list_wv == 800)-R_wv(:,:,list_wv == 675))./(R_wv(:,:,list_wv == 800)+R_wv(:,:,list_wv == 675)); % chl_a
    psnd(:,:,2) = (R_wv(:,:,list_wv == 800)-R_wv(:,:,list_wv == 650))./(R_wv(:,:,list_wv == 800)+R_wv(:,:,list_wv == 650)); % chl_b
    psnd(:,:,3)   = (R_wv(:,:,list_wv == 800)-R_wv(:,:,list_wv == 500))./(R_wv(:,:,list_wv == 800)+R_wv(:,:,list_wv == 500)); % car
	save(fullfile(indexPath, [fileName_, '_psnd.mat']), 'psnd')
   
    figure, subplot(1,3,1), imagesc(psnd(:,:,1)), axis image
    subplot(1,3,2), imagesc(psnd(:,:,2)), axis image
    subplot(1,3,3), imagesc(psnd(:,:,3)), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_psnd.png']), 'png')
    
    % REP (slope?)
    
    % RGRI
    rgri = R_red./R_green;
	save(fullfile(indexPath, [fileName_, '_rgri.mat']), 'rgri')
    figure, imagesc(rgri), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_rgri.png']), 'png')
    
    % RVSI
    rvsi = (R_wv(:,:,list_wv == 714)+R_wv(:,:,list_wv == 752))/2-R_wv(:,:,list_wv == 733);
    save(fullfile(indexPath, [fileName_, '_rvsi.mat']), 'rvsi')
    figure, imagesc(rvsi), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_rvsi.png']), 'png')
    
    % SAVI (L?)
    
    % SIPI
    sipi = (R_wv(:,:,list_wv == 800)-R_wv(:,:,list_wv == 445))./(R_wv(:,:,list_wv == 800)-R_wv(:,:,list_wv == 680));
	save(fullfile(indexPath, [fileName_, '_sipi.mat']), 'sipi')
    figure, imagesc(sipi), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_sipi.png']), 'png')
    
    % SR
    sr = R_wv(:,:,list_wv == 800)./R_wv(:,:,list_wv == 675);
	save(fullfile(indexPath, [fileName_, '_sr.mat']), 'sr')
    figure, imagesc(sr), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_sr.png']), 'png')
    
    % TSAVI (a, b?)
    
    % VARI
    vari = (R_green-R_red)./(R_green+R_red-R_blue);
    save(fullfile(indexPath, [fileName_, '_vari.mat']), 'vari')
    figure, imagesc(vari), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_vari.png']), 'png')
    
    % VI_green
    vi_green = (R_green-R_red)./(R_green+R_red); 
    save(fullfile(indexPath, [fileName_, '_vi_green.mat']), 'vi_green')
    figure, imagesc(vi_green), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_vi_green.png']), 'png')
    
    % WDVI (a?)
%     wdvi = R_nir
    
    % WBI
    wbi = R_wv(:,:,list_wv == 900)./R_wv(:,:,list_wv == 970);
    save(fullfile(indexPath, [fileName_, '_wbi.mat']), 'wbi')
    figure, imagesc(wbi), axis image
    saveas(gcf, fullfile(indexPath, [fileName_, '_wbi.png']), 'png')
    
    close all
end


