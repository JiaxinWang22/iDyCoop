% This code calculates inter-brain correlation between teammates
% also calculates partial correlation after regressing out VC and VG

%% correlation power (self) ~ power (teammate)
freqBand = 30:150; % frequency of interest
temp_NeuData = neu_data; % neural data for each dyad, channel x freq x time
temp_Beta = cell(size(temp_NeuData,1),1);
for iPair = 1:size(temp_NeuData,1) % dyads
    tempNeuW1 = squeeze(mean(temp_NeuData{iPair,1}(:,freqBand,:),2,'omitnan')); % player1
    tempNeuW2 = squeeze(mean(temp_NeuData{iPair,2}(:,freqBand,:),2,'omitnan')); % player2

    % size error handling when there is only 1 channel
    if size(tempNeuW1,2) == size(tempNeuW2,1)
        tempNeuW2 = tempNeuW2';
    elseif size(tempNeuW2,1) == size(tempNeuW1,2)
        tempNeuW1 = tempNeuW1';
    end
    % calculate correlation
    tempV = corr(tempNeuW1', tempNeuW2', 'rows', 'pairwise', 'type', 'Pearson');
    temp_Beta{iPair,1} = tempV(:);
end
temp_corr = cell2mat(temp_Beta);
% remove outliers
temp_corr(any(isoutlier(temp_corr,'mean','ThresholdFactor',3),2),:) = nan;
ins_corr = atanh(temp_corr); % Fisher-z transform

%% partial correlation power (self) ~ power (teammate) + GOM + COM
freqBand = 30:150; % frequency of interest
temp_NeuData = neu_data; % neural data for each dyad, channel x freq x time
temp_BehData = beh_data; % behavioral data for each dyad, variable x time
temp_Beta = cell(size(temp_NeuData));
for iPair = 1:size(temp_NeuData,1) % dyads
    tempBehW = temp_BehData{iPair,1};
    tempNeuW1 = squeeze(mean(temp_NeuData{iPair,1}(:,freqBand,:),2,'omitnan'));
    tempNeuW2 = squeeze(mean(temp_NeuData{iPair,2}(:,freqBand,:),2,'omitnan'));

    % size error handling when there is only 1 channel
    if size(tempNeuW1,2) == size(tempNeuW2,1)
        tempNeuW2 = tempNeuW2';
    elseif size(tempNeuW2,1) == size(tempNeuW1,2)
        tempNeuW1 = tempNeuW1';
    end

    corr = nan(size(tempNeuW1,1)*size(tempNeuW2,1), 1);
    for chani = 1:size(tempNeuW1,1)
        temp_corr = nan(size(tempNeuW2,1), 1);
        for chanj = 1:size(tempNeuW2,1)
            tempX = tempNeuW2(chanj,:)';
            tempY = tempNeuW1(chani,:)';
            tempZ = tempBehW';
            tempValid = ~isnan(tempX) & ~isnan(tempY) & ~isnan(tempZ(:,1));

            % calculate partial correlation
            if sum(tempValid) >= 2
                tempX = tempX(tempValid,:);
                tempY = tempY(tempValid,:);
                tempZ = tempZ(tempValid,:);
                tempX = znscore(tempX,1);
                tempY = znscore(tempY,1);
                tempZ = znscore(tempZ,1);
                temp_corr(chanj,:) = partialcorr(tempX, tempY, tempZ);
            else
                temp_corr(chanj,:) = nan;
            end
        end
        corr(chani*chanj-chanj+1:chani*chanj,:) = temp_corr;
    end
    temp_Beta{iPair,1} = corr;
end
temp_partial_corr = cell2mat([temp_Beta(:,1); temp_Beta(:,2)]); % concatenate subjects
temp_partial_corr(any(isoutlier(temp_partial_corr,'mean','ThresholdFactor',3),2),:) = nan;
ins_partial_corr = atanh(temp_partial_corr); % Fisher-z transform

%% z-score data with omitnan
function outputData = znscore(inputData,dim)
tempMean = mean(inputData,dim,'omitnan');
tempSTD = std(inputData,[],dim,'omitnan');
outputData = (inputData - tempMean)./tempSTD;
end
