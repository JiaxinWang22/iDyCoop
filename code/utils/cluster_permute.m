function cpinfo = cluster_permute(data1, data2, dpdt, thrP, numP, twoside, numC, Conn)
% 
% Reference: https://ww2.mathworks.cn/matlabcentral/fileexchange/71737-permutest
% Modified by Jiaxin Wang at SANP lab https://mylab.bnu.edu.cn/
% 
% cluster_permutation conducts permutation test with parallel computing
% INPUTS:
% data1: 2D/3D matrix with testing trials at the LAST dimension
% data2: scalar (for one-sample test) or 2D/3D matrix with testing trials at the LAST dimension
% dpdt: if true, data1 and data2 should have the same size, or data2 be a scalar
%       if false, other dimensions but the first one should be the same size
%       if data1 and data2 have the same size (or data2 be a scalar), default is true
%       if data1 and data2 have the same size but the first dimension, default is false
% thrP: scalar, p-value threshold for cluster significance, default is 0.05
% numP: integer, number of permutations, default is 10000
% side: true (to find negative clusters) or false, default is true
% numC: integer, number of clusters (sorted by summed t-values) to determine
%       significant (set to inf to include all clusters), default is inf
% Conn: cluster connectivity (when using function bwconncomp), default is 4
%
% OUTPUTS:
%
% Examples:
% [~,~,~,~] = permutest(rand(3,4,5), 0);
% [~,~,~,~] = permutest(rand(3,4,5), rand(3,4,5), true);
% [~,~,~,~] = permutest(rand(3,4,5), rand(6,4,5), false);
% [~,~,~,~] = permutest(rand(3,4,5), rand(6,4,5), false, 0.01, 5000, true, inf, 8);

%% Check inputs
checknan = [any(isnan(data1),'all') any(isnan(data2),'all')];
if any(checknan)
    warning('There is NaN in input data, will be ignored.')
end

if nargin < 1
    error('Not enough input arguments');
end
if nargin < 2 || isempty(data2)
    data2 = 0;
end

if ndims(data1) == 3
    [data1x, data1y, ntrials1] = size(data1);
elseif ismatrix(data1)
    [data1x, ntrials1] = size(data1);
    data1y = 1;
else
    error('Input data1 needs to be 2D or 3D');
end
if ndims(data2) == 3
    [data2x, data2y, ntrials2] = size(data2);
elseif isscalar(data2)
    data2x = data1x;
    data2y = data1y;
    ntrials2 = ntrials1;
elseif ismatrix(data2)
    [data2x, ntrials2] = size(data2);
    data2y = 1;
else
    error('Input data2 needs to be 2D or 3D or a scalar');
end

if nargin < 3 || isempty(dpdt)
    if isequal(size(data1),size(data2)) || isscalar(data2)
        dpdt = true;
        ntrials = ntrials1;
    elseif isequal(data1x,data2x) && isequal(data1y,data2y)
        dpdt = false;
    end
else
    if dpdt
        if ~isequal(size(data1),size(data2)) && ~isscalar(data2)
            error('Size of all dimensions should be identical for two dependent samples');
        end
        ntrials = ntrials1;
    else
        if ~isequal(data1x,data2x) || ~isequal(data1y,data2y)
            error('Size of all dimensions but the last one should be identical for two independent samples');
        end
    end
end
datax = data1x;
datay = data1y;

if nargin < 4 || isempty(thrP)
    thrP = 0.05;
end
if nargin < 5 || isempty(numP)
    numP = 10^4;
end
if nargin < 6 || isempty(twoside)
    twoside = true;
end
if nargin < 7 || isempty(numC)
    numC = inf;
end
if nargin < 8 || isempty(Conn)
    Conn = 4;
end

if dpdt
    maxnumP = 2^ntrials;
    if numP > maxnumP
        warning('Only %d permutations are possible. Using this value instead of %d.', maxnumP, numP);
        numP = maxnumP;
    end
else
    warning('off','MATLAB:nchoosek:LargeCoefficient');
    maxnumP = nchoosek(ntrials1 + ntrials2, ntrials1);
    warning('on','MATLAB:nchoosek:LargeCoefficient')
    if numP > maxnumP
        warning('Only %d permutations are possible. Using this value instead of %d', maxnumP, numP);
        numP = maxnumP;
    end
end

%% Initialization
clusters = cell(1);
pvalue = ones(1);
if twoside
    if dpdt
        tThreshold = abs(tinv(thrP/2, ntrials - 1));
    else
        tThreshold = abs(tinv(thrP/2, ntrials1 + ntrials2 - 2));
    end
else
    if dpdt
        tThreshold = abs(tinv(thrP, ntrials - 1));
    else
        tThreshold = abs(tinv(thrP, ntrials1 + ntrials2 - 2));
    end
end

%% PRODUCE PERMUTATION VECTORS
if numP < maxnumP / 1000
    if dpdt
        permutation_vectors = round(rand(ntrials,numP)) * 2 - 1;
    else
        permutation_vectors = ones(ntrials1+ntrials2,numP);
        for p = 1:numP
            idx = randperm(ntrials1+ntrials2,ntrials2);
            permutation_vectors(idx,p) = 2;
        end
    end
else
    if dpdt
        rndB = dec2bin(randperm(2^ntrials,numP) - 1);
        nBits = size(rndB,2);
        if nBits < ntrials
            rndB(:,(ntrials-nBits+1):ntrials) = rndB;
            rndB(:,1:ntrials-nBits) = '0';
        end
        permutation_vectors = (rndB == '1')' * 2 - 1;
    else
        permutation_vectors = ones(ntrials1 + ntrials2,numP);
        idx_matrix = nchoosek(1:(ntrials1 + ntrials2), ntrials2);
        idx_matrix = idx_matrix(randperm(size(idx_matrix, 1),numP),:)';
        for p = 1:numP
            permutation_vectors(idx_matrix(:,p),p) = 2;
        end
    end
end

%% RUN PRIMARY TTEST
if ~isscalar(data2)    
    if dpdt && datay == 1
        t_value_vector = modttest(data1'-data2',0);
        cpinfo.imgrdiff = mean(data1,2,'omitnan') - mean(data2,2,'omitnan');
        cpinfo.imgtdiff = t_value_vector;
    end
    
    if dpdt && datay > 1        
        t_value_vector = zeros(datax,datay);
        for ii = 1:datay
            t_value_vector(:,ii) = modttest(squeeze(data1(:,ii,:)-data2(:,ii,:))',0);
        end
        cpinfo.imgrdiff = mean(data1,3,'omitnan') - mean(data2,3,'omitnan');
        cpinfo.imgtdiff = t_value_vector;
    end
    
    if ~dpdt && datay == 1
        t_value_vector = modttest2(data1',data2');
        cpinfo.imgrdiff = mean(data1,2,'omitnan') - mean(data2,2,'omitnan');
        cpinfo.imgtdiff = t_value_vector;
    end
    
    if ~dpdt && datay > 1
        t_value_vector = zeros(datax,datay);
        for ii = 1:datay
            t_value_vector(:,ii) = modttest2(squeeze(data1(:,ii,:))',squeeze(data2(:,ii,:))');
        end
        cpinfo.imgrdiff = mean(data1,3,'omitnan') - mean(data2,3,'omitnan');
        cpinfo.imgtdiff = t_value_vector;
    end
    
else
    if dpdt && datay == 1
        t_value_vector = modttest(data1',data2);
        cpinfo.imgrzero = mean(data1-data2,2,'omitnan');
        cpinfo.imgtzero = t_value_vector;
    end
    
    if dpdt && datay > 1
        t_value_vector = zeros(datax,datay);
        for ii = 1:datay
            t_value_vector(:,ii) = modttest(squeeze(data1(:,ii,:))',data2);
        end
        cpinfo.imgrzero = mean(data1-data2,3,'omitnan');
        cpinfo.imgtzero = t_value_vector;
    end
end

% Find the above-threshold clusters:
connectivity = Conn;
CC = bwconncomp(t_value_vector > tThreshold, connectivity);
cMapPrimary = zeros(size(t_value_vector));
tSumPrimary = zeros(CC.NumObjects,1);
for i = 1:CC.NumObjects
    cMapPrimary(CC.PixelIdxList{i}) = i;
    tSumPrimary(i) = sum(t_value_vector(CC.PixelIdxList{i}));
end
if twoside
    n = CC.NumObjects;
    CC = bwconncomp(t_value_vector < -tThreshold, connectivity);
    for i = 1:CC.NumObjects
        cMapPrimary(CC.PixelIdxList{i}) = n + i;
        tSumPrimary(n+i) = sum(t_value_vector(CC.PixelIdxList{i}));
    end
end
[~,tSumIdx] = sort(abs(tSumPrimary),'descend');
tSumPrimary = tSumPrimary(tSumIdx);

%% RUN PERMUTATIONS
permutation_distribution = zeros(numP,1);
if ~isscalar(data2)
    if dpdt && datay == 1
        tempdata1 = data1';
        tempdata2 = data2';
        parfor p = 1:numP
            D = bsxfun(@times,tempdata1-tempdata2,permutation_vectors(:,p));
            t_value_vector = modttest(D,0);
            permutation_distribution(p) = TmaxC(t_value_vector, tThreshold, connectivity, twoside);
        end
    end
    
    if dpdt && datay > 1
        tempdata1 = permute(data1,[3 1 2]);
        tempdata2 = permute(data2,[3 1 2]);
        parfor p = 1:numP
            t_value_vector = zeros(datax,datay);
            for ii = 1:datay
                D = tempdata1(:,:,ii)-tempdata2(:,:,ii);
                D = bsxfun(@times,D,permutation_vectors(:,p));
                t_value_vector(:,ii) = modttest(D,0);
            end
            permutation_distribution(p) = TmaxC(t_value_vector, tThreshold, connectivity, twoside);
        end
    end
    
    if ~dpdt && datay == 1
        all_trials = cat(ndims(data1), data1, data2);
        parfor p = 1:numP
            tempdata1 = all_trials(:,permutation_vectors(:,p)==1);
            tempdata2 = all_trials(:,permutation_vectors(:,p)==2);
            t_value_vector = modttest2(tempdata1',tempdata2');
            permutation_distribution(p) = TmaxC(t_value_vector, tThreshold, connectivity, twoside);
        end
    end
    
    if ~dpdt && datay > 1
        all_trials = cat(ndims(data1), data1, data2);
        parfor p = 1:numP
            t_value_vector = zeros(datax,datay);
            tempdata1 = all_trials(:,:,permutation_vectors(:,p)==1);
            tempdata2 = all_trials(:,:,permutation_vectors(:,p)==2);
            for ii = 1:datay
                t_value_vector(:,ii) = modttest2(squeeze(tempdata1(:,ii,:))',squeeze(tempdata2(:,ii,:))');
            end
            permutation_distribution(p) = TmaxC(t_value_vector, tThreshold, connectivity, twoside);
        end
    end
else
    if dpdt && datay == 1
        tempdata1 = data1';
        parfor p = 1:numP
            D = bsxfun(@times,tempdata1,permutation_vectors(:,p));
            t_value_vector = modttest(D,data2);
            permutation_distribution(p) = TmaxC(t_value_vector, tThreshold, connectivity, twoside);
        end
    end
    
    if dpdt && datay > 1
        tempdata1 = permute(data1,[3 1 2]);
        parfor p = 1:numP
            t_value_vector = zeros(datax,datay);
            for ii = 1:datay
                D = tempdata1(:,:,ii);
                D = bsxfun(@times,D,permutation_vectors(:,p));
                t_value_vector(:,ii) = modttest(D,data2);
            end
            permutation_distribution(p) = TmaxC(t_value_vector, tThreshold, connectivity, twoside);
        end
    end
    
end
%% DETERIMNE SIGNIFICANCE
for clustIdx = 1:min(numC,length(tSumPrimary))
    if twoside
        ii = sum(abs(permutation_distribution) >= abs(tSumPrimary(clustIdx)));
    else
        ii = sum(permutation_distribution >= tSumPrimary(clustIdx));
    end
    
    clusters{clustIdx,1} = find(cMapPrimary == tSumIdx(clustIdx));
    clusters{clustIdx,1} = clusters{clustIdx,1}(:);
    pvalue(clustIdx,1) = (ii+1) / (numP+1);
end
cpinfo.cluster = clusters;
cpinfo.pvalue = pvalue;
cpinfo.permutation_distribution = permutation_distribution;
cpinfo.t_sums = tSumPrimary;
% disp('Permutation finished.')

end

%% Functions
function t = modttest(x,m)
if nargin < 2 || isempty(m)
    m = 0;
end
samplesize = size(x,1);
sqrtn = sqrt(samplesize);
xmean = sum(x,1,'omitnan')/samplesize;
xc = bsxfun(@minus,x,xmean);
xstd = sqrt(sum(conj(xc).*xc,1,'omitnan')/(samplesize-1));
xdiff = xmean - m;

% Check for rounding issues causing spurious differences (from ttest function)
fix = (xdiff ~= 0) & (abs(xdiff) < 100*sqrtn.*max(eps(xmean), eps(m)));
if any(fix(:))
    constvalue = min(x,[],1);
    fix = fix & all(x == constvalue, 1);
end
if any(fix(:))
    xdiff(fix) = 0;
    xstd(fix) = 0;
end

ser = xstd ./ sqrtn;
t = xdiff ./ ser;
end

function t = modttest2(x1,x2)

n1 = size(x1,1);
n2 = size(x2,1);

xmean1 = sum(x1,1,'omitnan')/n1;
xmean2 = sum(x2,1,'omitnan')/n2;

xc = bsxfun(@minus,x1,xmean1);
xstd1 = sqrt(sum(conj(xc).*xc,1,'omitnan')/(n1-1));
xc = bsxfun(@minus,x2,xmean2);
xstd2 = sqrt(sum(conj(xc).*xc,1,'omitnan')/(n2-1));

sx1x2 = sqrt(xstd1.^2/n1 + xstd2.^2/n2);

t = (xmean1 - xmean2) ./ sx1x2;

end

function pd = TmaxC(t_value, tThr, connect, side)

CC = bwconncomp(t_value > tThr, connect);
tSum = zeros(CC.NumObjects,1);
for i = 1:CC.NumObjects
    tSum(i) = sum(t_value(CC.PixelIdxList{i}));
end
if side
    n = CC.NumObjects;
    CC = bwconncomp(t_value < -tThr, connect);
    for i = 1:CC.NumObjects
        tSum(n+i) = sum(t_value(CC.PixelIdxList{i}));
    end
end

if isempty(tSum)
    pd = 0;
else
    [~,idx] = max(abs(tSum));
    pd = tSum(idx);
end
end
