clear; clc;

% Define parameter ranges
eta_range = [0,1,3,5];
K_range = [50,100,150,200];
delta_range = [0.001,0.01,0.1,0,1];
beta_range = [0.01,0.1,1];
r_range = [0,1,3,5];
% eta_range = 3;
% K_range = 150;
% delta_range = 0.001;
% beta_range = 0.01;
% r_range = 5;

% Open a file to save results
fid = fopen('新数据测试.txt', 'w');
path_x13 = 'Protein1.csv';
path_x14 = 'RNA5001.csv'; 
y = 'truth_10X10k.csv';  
% Read the CSV files into MATLAB
x13 = readmatrix(path_x13);
x14 = readmatrix(path_x14);  
gt=readmatrix(y);  
save('x13.mat', 'x13');
save('x14.mat', 'x14');
save('y.mat', 'y');
% Loop through parameter combinations
for eta = eta_range
    for K = K_range
        for delta = delta_range
            for beta = beta_range
                for r = r_range
                    try
                        rand('seed', 900);
%                         load('SMAGE.mat');
%                         for i = 1:size(P, 2)
%                             X{i} = double(P{i}');
%                         end
%                         gt = double(real_label);
%                         load('pbmc_3kY.mat');
%                         gt=data';
                        path_x13 = 'Protein1.csv';
                        path_x14 = 'RNA5001.csv'; 
                        y = 'truth_10X10k.csv';  
                        % Read the CSV files into MATLAB
                        x13 = readmatrix(path_x13);
                        x14 = readmatrix(path_x14);  
                        gt=readmatrix(y);  
                        X={x13,x14};
%                         for i=1:size(X,2)
%                             X{i}=NormalizeData(X{i});
%                         end
                        option.eta = eta;
                        option.K = K;
                        option.delta = delta;
                        option.beta = beta;
                        option.r = r;
                        option.Vnum = size(X, 2);  % Number of views
                        option.alpha = ones(option.Vnum, 1) / option.Vnum;
                        option.numClust = size(unique(gt), 1);  % Number of classes
                        option.threshold = 1e-1;  % Threshold
                        option.lambda = 1e-4;  % Regularization parameter
                        option.max_iter = 20;  % Max iterations
                        option.N = size(X{1}, 2);  % Data dimension

                        % Call model
                        [Y, Q, C, G, U] = PLCMF(X, gt, option);
                        numClust = size(unique(gt), 1); 
                        clusterLabels= kmeans(G', numClust);
                        permutedLabels = bestMap(gt, clusterLabels);
                        accuracy = sum(permutedLabels == gt) / length(gt);
                        nmi = compute_NMI(gt, permutedLabels);
                        ami = AMI(gt, permutedLabels);  % Calculate AMI
                        ari = ARI(gt, max(gt), permutedLabels, max(permutedLabels));  % Calculate ARI
% 
%                         for e = 1:size(G, 2)
%                             v = G(:, e);
%                             ma = max(v);
%                             [s, t] = find(v == ma);
%                             l(e) = s;
%                         end
%                         l = l';  % Labels obtained by DRjCC
%                         ll = gt;
%                         [newl] = bestMap(ll, l);  % Permute label of l
%                         nmi = compute_NMI(ll, newl);  % Calculate NMI
%                         ami = AMI(ll, newl);  % Calculate AMI
%                         ari = ARI(ll, max(ll), newl, max(newl));  % Calculate ARI
%                         pre_label = newl;
%                         if ~isempty(ll)
%                             exact = find(pre_label == ll);
%                             accuracy = length(exact) / length(newl);  % Calculate accuracy
%                         else
%                             accuracy = [];
%                         end
                        result = [accuracy, nmi, ami, ari]
                        
                        % Save result to the file
%                         fprintf(fid, 'eta: %d, K: %d, delta: %.2f, beta: %.2f, r: %d, results: [%s]\n', ...
%                             eta, K, delta, beta, r, num2str(result));
                        fprintf(fid, 'Parameters: eta: %d, K: %d, delta: %.2f, beta: %.2f, r: %d | Results: Accuracy: %.4f, NMI: %.4f, AMI: %.4f, ARI: %.4f\n', ...
                           eta, K, delta, beta, r, result(1), result(2), result(3), result(4));

                    catch ME
                        % Skip this parameter set and print the error message
                        fprintf('Error with parameters eta: %d, K: %d, delta: %.2f, beta: %.2f, r: %d. Skipping...\n', ...
                            eta, K, delta, beta, r);
                        fprintf('Error message: %s\n', ME.message);
                    end
                end
            end
        end
    end
end

% Close the results file
fclose(fid);
