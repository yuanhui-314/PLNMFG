clc;clear;
% 初始化参数
eta_range = 0;
K_range = 40;
delta_range = 0.001;
beta_range = 0.01;
r_range = 5;
alpha1_range = 1;

for eta = eta_range
    for K = K_range
        for delta = delta_range
            for beta = beta_range
                for r = r_range
                    for alpha1 = alpha1_range
                        try
                            rand('seed', 900);
                            load('anno.mat');
                            for i = 1:size(P, 2)
                                X{i} = double(P{i}');
                            end
                            gt = double(real_label);
                            option.eta = eta;
                            option.K = K;
                            option.delta = delta;
                            option.beta = beta;
                            option.r = r;
                            option.Vnum = size(X, 2);  
                            option.alpha = ones(option.Vnum, 1) / option.Vnum;
                            option.numClust = size(unique(gt), 1);  
                            option.threshold = 1e-1; 
                            option.lambda = 1e-4; 
                            option.max_iter = 20;  
                            option.N = size(X{1}, 2);  
                            option.alpha1 = alpha1; 

                            % Call model
                            [Y, Q, C, G, U] = PLNMFG(X, gt, option);
                            numClust = size(unique(gt), 1); 
                            clusterLabels = kmeans(G', numClust);
                            permutedLabels = bestMap(gt, clusterLabels);
                            accuracy = sum(permutedLabels == gt) / length(gt);
                            nmi = compute_NMI(gt, permutedLabels);
                            ami = AMI(gt, permutedLabels);  
                            ari = ARI(gt, max(gt), permutedLabels, max(permutedLabels)); 
                            result = [accuracy, nmi, ami, ari]                        

                        end
                    end
                end
            end
        end
    end
end

