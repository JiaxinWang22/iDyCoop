function epoch = find_epoch(inputData, inputType, outputType)
% extract epochs with consective time-stamps in a vector
% inputType:  numeric   - integers
%             binary    - 0 or 1 
% outputType: valuecell - extract epochs in cells
%             valuemat  - extract the start and end element of epochs
%             indexmat  - extract the idx of the start and end element of epochs
% NOTE:       binary input can only have indexmat output 
%             inputData should be a vector
% Written by Jiaxin Wang at SANP lab https://mylab.bnu.edu.cn/

if isempty(inputData)
    error('input data is empty.')
elseif ~isvector(inputData)
    error('input data should be a row or column vector.')
end
inputData = double(inputData(:));

switch inputType
    case 'numeric'
        B = [0; find(diff(inputData) ~= 1); numel(inputData)];
    case 'binary'
        B = [0; find(diff(inputData) ~= 0); numel(inputData)];
        if strcmp(outputType,'valuemat') || strcmp(outputType,'valuecell')
            outputType = 'indexmat';
            warning('output can only be ''indexmat'' for binary input.');
        end
    otherwise
        error('Method not supported.')
end

switch outputType
    case 'valuecell'
        epoch = cell(length(B) - 1, 1);
        for i = 1:length(B) - 1
            epoch{i,1} = inputData((B(i)+1) : B(i+1));
        end
    case 'valuemat'
        epoch = [inputData(B(1:end-1) + 1), inputData(B(2:end))];
    case 'indexmat'
        epoch = [B(1:end-1) + 1, B(2:end)];
        if strcmp(inputType,'binary')
            if inputData(1) == 1
                epoch(2:2:end,:) = [];
            else
                epoch(1:2:end,:) = [];
            end
        end

    otherwise
        error('Output format not supported.')
end

end