
% Step 1: 加载多组学数据
% RNA 数据 (行: 样本, 列: 特征)
rna_data = csvread('RNA5001.csv');
% 蛋白质数据
protein_data = csvread('Protein1.csv');
rna_data = rna_data';
protein_data = protein_data';
size(rna_data)
size(protein_data)

% Step 2: 数据整合 (按列合并)
% 假设多组学数据有相同的样本顺序
multiomics_data = [rna_data, protein_data]; 

% Step 4: 运行 UMAP
% 参数解释: n_neighbors 控制局部区域大小, min_dist 控制点间距离
[reduction, umap_params] = run_umap(multiomics_data, ...
    'n_neighbors', 15, 'min_dist', 0.1);
reduction
% Step 5: 聚类分析 (使用 k-means 聚类)
num_clusters = 8;
% 假设reduction是UMAP降维后的数据
if size(reduction, 1) > num_clusters
    % 确保数据集行数大于簇数
    cluster_labels = kmeans(reduction, num_clusters);
else
    disp('数据集的样本数小于簇数，请减少簇的数量');
end

% Step 6: UMAP 可视化
figure;
scatter(reduction(:, 1), reduction(:, 2), 30, cluster_labels, 'filled');
colormap(jet(num_clusters)); % 设置颜色映射
title('UMAP Clustering of Multi-Omics Data');
xlabel('UMAP Dimension 1');
ylabel('UMAP Dimension 2');
colorbar;

