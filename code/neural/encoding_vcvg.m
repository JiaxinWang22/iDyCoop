% This code calculates partial correlation between Power and VC/VG
% p values are determined by a permutation method for each channel

permuteN = 5000;
freqBand = 30:110; % frequency of interest
temp_NeuData = tempNeu;  % neural data, channel x freq x time
temp_BehData = tempBeh;  % behavioral data, variable x time

temp_Beta = cell(numel(temp_NeuData), size(temp_BehData{1,1},1));
temp_PVal = cell(numel(temp_NeuData), size(temp_BehData{1,1},1));
temp_FVal = cell(numel(temp_NeuData), size(temp_BehData{1,1},1));
temp_FVal_Permute = cell(numel(tempNeuData), size(tempBehData{1,1},1));

for vari = 1:size(temp_BehData{1,1},1) % VC and VG
    for subi = 1:numel(temp_NeuData) % subject
        tempBehW = temp_BehData{subi,1};
        tempNeuW = squeeze(mean(temp_NeuData{subi,1}(:,freqBand,:), 2, 'omitnan'));

        % size error handling when there is only 1 channel
        if size(tempNeuW,1) == size(tempBehW,2)
            tempNeuW = tempNeuW';
        end

        % stat variables
        v = nan(size(tempNeuW,1), 1); 
        p = nan(size(tempNeuW,1), 1); 
        F = nan(size(tempNeuW,1), 1); 
        F_permute = nan(permuteN, size(tempNeuW,1)); 
        for chani = 1:size(tempNeuW,1) % channel
            tempX = tempNeuW(chani,:)';             % neural time series
            tempY = tempBehW(vari,:)';              % behavioral time series
            tempZ = tempBehW(setdiff(1:2,vari),:)'; % behavioral control variable
            tempValid = ~isnan(tempX) & ~isnan(tempY) & ~isnan(tempZ);

            if sum(tempValid) >= 2 % ensure there are more than 3 time points
                tempX = tempX(tempValid,:);
                tempY = tempY(tempValid,:);
                tempZ = tempZ(tempValid,:);
                tempX = znscore(tempX,1);
                tempY = znscore(tempY,1);
                tempZ = znscore(tempZ,1);

                % calculate partial correlation, same as using partialcorr
                b1 = regress(tempX, [ones(numel(tempZ),1) tempZ]);
                residualX = tempX - b1(2) * tempZ;
                b2 = regress(tempY, [ones(numel(tempZ),1) tempZ]);
                residualY = tempY - b2(2) * tempZ;
                residualX = znscore(residualX,1);
                residualY = znscore(residualY,1);
                [b,~,~,~,stats] = regress(residualX, [ones(numel(residualY),1) residualY]);
                v(chani,1) = b(2);
                F(chani,1) = stats(2);
                p(chani,1) = stats(3);

                % shuffle time series and calculate F distribution
                parfor permutei = 1:permuteN
                    % shuffle neural time series
                    tempXShuffle = tempX(randperm(numel(tempX)));

                    b1 = regress(tempXShuffle, [ones(numel(tempZ),1) tempZ]);
                    residualX = tempXShuffle - b1(2) * tempZ;
                    b2 = regress(tempY, [ones(numel(tempZ),1) tempZ]);
                    residualY = tempY - b2(2) * tempZ;
                    residualX = znscore(residualX,1);
                    residualY = znscore(residualY,1);
                    [~,~,~,~,stats] = regress(residualX, [ones(numel(residualY),1) residualY]);
                    F_permute(permutei,chani) = stats(2);
                end
            else
                v(chani,1) = nan;
                F(chani,1) = nan;
                p(chani,1) = nan;
            end
        end
        temp_Beta{subi,vari} = v;
        temp_PVal{subi,vari} = p;
        temp_FVal{subi,vari} = F;
        temp_FVal_Permute{subi,vari} = F_permute;
    end
end

% calculate p value
temp_PVal_Permute = cell(size(temp_FVal_Permute));
for vari = 1:size(tempStdFVal_Permute,2)
    for subi = 1:size(tempStdFVal_Permute,1)
        F = temp_FVal_Permute{subi,vari};
        p = nan(size(temp_FVal_Permute{subi,vari},2), 1);

        % calculate num of permute value that are larger than the true value 
        for chani = 1:size(temp_FVal_Permute{subi,vari},2)
            if ~all(isnan(F(:,chani,1)))
                % two-sided test
                counts = sum(abs(F(:,chani,1)) >= abs(temp_FVal{subi,vari}(chani,1)));
                p(chani, 1) = (counts + 1) / (permuteN + 1);
            end
        end
        temp_PVal_Permute{subi,vari} = p;
    end
end

% FDR correction across all channels
p_value = nan(size(cell2mat(temp_PVal_Permute(:,1)),1), size(temp_PVal_Permute,2));
for vari = 1:size(tempStdFVal_Permute,2)
    p_value(:,vari) = cell2mat(temp_PVal_Permute(:,vari));
end
pval_grpFDR = [];
[~,~,~,pval_grpFDR(:,1)] = fdr_bh(p_value(:,1));
[~,~,~,pval_grpFDR(:,2)] = fdr_bh(p_value(:,2));

% remove outliers of beta value
tempOutlier = [];
tempBeta = nan(size(cell2mat(temp_Beta(:,1)),1), size(temp_Beta,2));
for vari = 1:size(temp_Beta,2)
    tempCol = cell2mat(temp_Beta(:,vari));
    tempBeta(:,vari) = tempCol;
    tempOutlier{vari,1} = find(isoutlier(tempBeta(:,vari),'mean','ThresholdFactor',3));
    tempBeta(tempOutlier{vari,1}, vari) = nan;
end
encoding_beta = atanh(tempBeta); % Fisher-z transform

%% z-score data with omitnan
function outputData = znscore(inputData,dim)
tempMean = mean(inputData,dim,'omitnan');
tempSTD = std(inputData,[],dim,'omitnan');
outputData = (inputData - tempMean)./tempSTD;
end
