% build LME model for each behavioral variable
% predictors: behavioral variables
% output: trial win/lose
% data is NxV matrix, with N as trials, V as variables
tempBehVar = znscore(data,1);

% make into a table for model fitting
TrialNum = categorical(tempYWin(:,2));
SubPairNum = categorical(tempYWin(:,3));
tbl1 = table(tempYWin(:,1), TrialNum, SubPairNum, 'VariableNames', {'WinLose','TrialNum','DyadNum'});
tbl2 = array2table(tempBehVar, 'VariableNames', {'MRatio','IRatio','resetN','VC','VG','VCDiff','VGDiff','AC','AG','Vel','Acc'});
tbl = [tbl1 tbl2];
varNames = tbl.Properties.VariableNames(4:end)';

% stats variables
temp_Beta = nan(size(varNames,1),1);
temp_CI = nan(size(varNames,1),2);
temp_PValue = nan(size(varNames,1),1);
temp_TValue = nan(size(varNames,1),1);
temp_AdjustedR = nan(size(varNames,1),1);
for vari = 1:length(varNames)
    % fit the model
    tempVarName = varNames{vari};
    tempFormula = ['WinLose ~ ' tempVarName ' + (1|DyadNum) + (1|DyadNum:TrialNum)'];
    tempFit = fitglme(tbl,tempFormula,'Distribution','Binomial','Link','logit','FitMethod','MPL');

    % stats
    temp_Beta(vari,1) = tempFit.Coefficients.Estimate(2:end);
    temp_CI(vari,:) = [tempFit.Coefficients.Lower(2:end) tempFit.Coefficients.Upper(2:end)];
    temp_TValue(vari,1) = tempFit.Coefficients.tStat(2:end);
    temp_PValue(vari,1) = tempFit.Coefficients.pValue(2:end);
    temp_AdjustedR(vari,1) = tempFit.Rsquared.Adjusted;
end
% FWE correction for p values
pvalue_FWE = temp_PValue * length(temp_PValue);

%% z-score data with omitnan
function outputData = znscore(inputData,dim)
tempMean = mean(inputData,dim,'omitnan');
tempSTD = std(inputData,[],dim,'omitnan');
outputData = (inputData - tempMean)./tempSTD;
end
