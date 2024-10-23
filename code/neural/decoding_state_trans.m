% decode state transition from power data
useTimeEndList = -99:4:100; % time points around transition points to conduct decoding
freqBand = 30:150; % power frequency bands
CVFold = 10; % 10-fold cross-validation

%% decoding
accuracy = nan(numel(useTimeEndList), CVFold);
accuracy_ch = nan(numel(useTimeEndList), CVFold); % for chance level

for time_i = 1:numel(useTimeEndList)
    temp_tran_data = tran_data{time_i,1};
    temp_nontran_data = nontran_data{time_i,1};

    % tran_data size is epoch x channel x freq x time
    temp_subdata_class1 = cellfun(@(c) permute(mean(c(:,:,freqBand,:),[2 3],'omitnan'), [1,4,2,3]), temp_tran_data, 'UniformOutput', 0);
    temp_subdata_class2 = cellfun(@(c) permute(mean(c(:,:,freqBand,:),[2 3],'omitnan'), [1,4,2,3]), temp_nontran_data, 'UniformOutput', 0);

    temp_data_class1 = [];
    temp_data_class2 = [];
    for subi = 1:size(temp_subdata_class1,1)
        % normalize within subject
        tempSubList = normalize([temp_subdata_class1{subi,1}; temp_subdata_class2{subi,1}], 1);

        % concatenate across subjects
        temp_data_class1 = [temp_data_class1; tempSubList(1:end/2,:)];
        temp_data_class2 = [temp_data_class2; tempSubList(end/2+1:end,:)];
    end

    temp_acc = nan(CVFold,1);
    temp_acc_ch = nan(CVFold,1);
    temp_label = [ones(size(temp_data_class1,1),1); -ones(size(temp_data_class1,1),1)];
    temp_input = [temp_data_class1; temp_data_class2];
    temp_CV = cvpartition(temp_label,'KFold',CVFold);

    for CVi = 1:CVFold
        temp_train_data  = temp_input(training(temp_CV, CVi),:);
        temp_test_data   = temp_input(test(temp_CV, CVi),:);
        temp_train_label = temp_label(training(temp_CV, CVi));
        temp_test_label  = temp_label(test(temp_CV, CVi));

        tempFitInfo = fitcsvm(normalize(temp_train_data,1), temp_train_label,...
            'KernelFunction', 'linear',...
            'CacheSize', 'maximal',...
            'Standardize', false, ...
            'KernelScale','auto');
        tempTestPredict = predict(tempFitInfo, normalize(temp_test_data,1));
        temp_acc(CVi,1) = sum(tempTestPredict == temp_test_label)./numel(temp_test_label);

        % chance level classifier: random shuffle training lables
        tempFitInfo_ch = fitcsvm(normalize(temp_train_data,1), temp_train_label(randperm(numel(temp_train_label)),:),...
            'KernelFunction','linear',...
            'CacheSize', 'maximal',...
            'Standardize', false, ...
            'KernelScale','auto');
        tempTestPredict = predict(tempFitInfo_ch, normalize(temp_test_data,1));
        temp_acc_ch(CVi,1) = sum(tempTestPredict == temp_test_label)./numel(temp_test_label);
    end

    accuracy(time_i, CVRi, :) = temp_acc;
    accuracy_ch(time_i, CVRi, :) = temp_acc_ch;
end
save('decoding_data.mat', 'accuracy', 'accuracy_ch', 'time_i');
