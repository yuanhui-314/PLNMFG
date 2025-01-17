function   [Y,Q,C,G,U]=PLNMFG(X,gt,option)
%%----------------Initialize-------------------
numClust=option.numClust ;
K=option.K;
threshold=option.threshold;
delta=option.delta;
lambda=option.lambda;         
r=option.r;
max_iter=option.max_iter;
Vnum=option.Vnum;
N=option.N ;                            
alpha=option.alpha; 
alpha1=option.alpha1; 
beta=option.beta; 
eta=option.eta;
alpha_r=alpha.^r;
Jlast=99999;
IsConverge = 0;
iter = 1;
U=cell(size(X,1),size(X,2));
Y=cell(size(X,1),size(X,2));
V1=zeros(K,K);
V2=zeros(K,N);
V=rand(K,N);
J1=ones(Vnum,1);
J2=ones(Vnum,1);
Ja=ones(Vnum,1);
maxDnorm = Inf;  % 初始化最大范数
tic;
% Initialize arrays to store loss and accuracy
loss_values = zeros(max_iter, 1);
ACC_values = zeros(max_iter, 1);
for i=1:Vnum 
    Q{i}=rand(numClust,K);
    U{i}=rand(size(X{i},1),K);
    Y{i}=rand(numClust,N);
    YY=litekmeans(X{i}',numClust,'MaxIter', 100);
    Y{i}=ToM(YY,size(Y{i},1),size(Y{i},2));
end 
for i=1:Vnum
    window{i}=(X{i}==0); 
    dim{i}=size(X{i});
    rand('seed',23);
    S0{i}=rand(dim{i}).*window{i};
end
C=abs(rand(K,numClust));
G=abs(rand(numClust,N));
Norm = 2;
NormV = 1;
[C,G] = NormalizeUV(C, G', NormV, Norm);G=G';

%%------------------Update---------------------------
loss_values = zeros(max_iter, 1); 
iter_indices = 1:max_iter;       
while (IsConverge == 0&&iter<max_iter+1) 
    %---------UpdateQi---------
     for i=1:Vnum
        Q{i}=Y{i}*V'/(V*V');
     end
     %---------UpdateUi---------
     for i=1:Vnum
        U{i}=(X{i}+S0{i})*V'/(V*V');
     end     
    %---------UpdateV---------
    V1=zeros(K,K);
    V2=zeros(K,N);
%     alpha_r(i)=1;
     for i=1:Vnum
        V1=V1+(alpha_r(i)*U{i}'*U{i}+alpha_r(i)*delta*Q{i}'*Q{i});
        V2=V2+(alpha_r(i)*U{i}'*(X{i}+S0{i})+alpha_r(i)*delta*Q{i}'*Y{i});
     end
     V=(V1+(beta+lambda)*eye(K))\(V2+beta*C*G);   
     options = [];
     option.Metric = 'Cosine';
     options.NeighborMode = 'KNN';
     options.k =5;%5 nearest neighbors
     options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
     A=constructA(V',options);
     Dbar = diag(sum(double(A),2));
     L=Dbar-A;
     clear options;

     %---------UpdateC---------
     C=V*G'/(G*G');
     C = max(C, 0);

     %---------UpdateG---------
	G1 = beta*C'*V+alpha1*G*A;
	G2 = beta*C'*C*G + alpha1*G*Dbar;
	G = G.*(G1)./(G2);
    G = max(G,0);
    [C,G] = NormalizeUV(C, G', NormV, Norm);
    G=G';


     %---------UpdateS---------
     for i=1:Vnum
        temp=U{i}*V-X{i};
        miu=diag(sum(X{i})/median(sum(X{i})));
        Skl=SoftThreshold(temp,ones(dim{i})*(eta*miu));
        S{i}=Skl.*window{i};
        S0{i}=S{i};
     end

    %---------Calculate J---------
        %% Calculate Loss
    Jcurrent = 0; 
    for i = 1:Vnum
        term2 = 0;
        for j = 1:min(size(S{i}, 2), numel(miu))
            term2 = term2 + miu(j) * norm(S{i}(:, j), 1);
        end
        J1(i)=alpha_r(i)*norm(X{i}+S{i}-U{i}*V,'fro')^2;
        J2(i)=alpha_r(i)*norm(Y{i}-Q{i}*V,'fro')^2;
        Ja(i)=J1(i)+delta*J2(i)+alpha_r(i)*eta*term2; 
        term1 = norm(X{i} + S{i} - U{i} * V, 'fro')^2;
        term3 = norm(Y{i} - Q{i} * V, 'fro')^2;
        Jcurrent = Jcurrent + (alpha_r(i) * (term1 + eta * term2 + delta * term3));
    end
    % Regularization term
    Jcurrent = Jcurrent + beta * norm(V - C * G, 'fro')^2+alpha1*trace(G*L*G');
    loss_values(iter) = Jcurrent; 
    disp(Jcurrent);

    %---------Calculate alpha---------   
    alpha=(Ja.^(1/(1-r)))/sum(Ja.^(1/(1-r)));
    alpha_r=alpha.^r;
    
%     ---------Iscoverage---------   
    J(iter)=Jcurrent;
    if (abs(Jlast - Jcurrent)) < threshold
            IsConverge=1;
    end

    if isnan(maxDnorm)
        fprintf('maxDnorm为NaN，停止执行：迭代 #%d\n', iter);
        result = NaN(1, 6); 
        return;
    end
    maxDnorm = Jcurrent;
    iter = iter + 1;
end
end
