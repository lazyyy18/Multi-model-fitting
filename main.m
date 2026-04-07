clc , clear , close all;
vlf;
addpath('./model_specific/')
addpath('./data/adelaidermf/H')
addpath('./data/adelaidermf/F')
datafile = 'F';
if datafile == 'F'
    matFiles = dir(fullfile('.\data\adelaidermf\F','*.mat'));
elseif datafile == 'H'
    matFiles = dir(fullfile('.\data\adelaidermf\H','*.mat'));
end
numfiles = length(matFiles);
min_error = [];
t =1;
average_time = 0;
nan = 0;
average_error = 0;
recall = 0;
precision = 0;
f1_score = 0;
average_time = 0;
average_error = 0;
for iter = 1 : numfiles
    % rng(42);
    total_time = 0;
    total_error = 0;
    close all
    residual_result = [];
    ttime = [];
    
    label_result = [];

    % % 讀資料
    disp(['running seq: ', matFiles(iter).name])
    load(matFiles(iter).name);
    AA = matFiles(iter).name(1:end-4);
    [~, base] = fileparts(matFiles(iter).name);   % -> '1_1.5' 或 '1'
    numModels = max(label) - min(label);
     data1 = data'; 
    [data, ia, ic] = unique(data1, "rows");
    data = data';
    numPoints=[];
    for i=min(label):max(label)
        numPoints = [ numPoints , sum(label==i)];
    end
    disp(['Num Points(outliersFirst): ', num2str(numPoints)])
    disp(['numModels: ',num2str(numModels)])

    N = sum(numPoints);
    X = data(1:3,:);
    Y = data(4:6,:);
    source = X(1:3,:)';
    target = Y(1:3,:)';
    if datafile == 'H'
        [time,est_label,model] = GH_plane_coosac(source, target, numModels, img1, img2);
    else
        [time,est_label,model] = GH_motion_coosac(source, target, numModels, img1, img2);
    end
    new_estlabel = est_label(ic);
    
    total_time = total_time + time;
    [mean_error,index] = missclass(new_estlabel, label);
    new_elabel = zeros(length(new_estlabel),1);
    for i=1:max(new_estlabel)
        new_elabel(new_estlabel == index(i+1)) = i;
    end
    new_estlabel = new_elabel;
    average_error = average_error + mean_error;
    average_time = average_time + time;
    fprintf('time ： %.2f\n', total_time);
    fprintf('error ： %.6f\n', mean_error);
end
min_error = [min_error, average_error];
fprintf('average time ： %.2f\n', average_time/ numfiles);
fprintf('mean_error ： %.6f\n', average_error/numfiles);

