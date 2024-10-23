function pcor = power_corr_sw(input1, input2, slidw, overlap, method)
% power_corr_sw calculates power correlation using a sliding window with overlap using parallel computing (parfor)
% pcor = power_corr_sw(input1, input2, slidw, overlap, method)
%
% INPUT
% input1 and input2 should be 3D CHANNEL x FREQUENCY x TIME matrices with same size of frequency and time
% slidw & overlap should be time length in SAMPLES (depending on your sampling rate)
% method supported: 'pearson' or 'spearman' (default is 'pearson')
%
% OUTPUT
% corspectrm is the power correlation matrix
% 
% NOTICE
% The output time length is shorter than the input time length, so you can add padding to the input data
%
% Written by Jiaxin Wang at SANP lab https://mylab.bnu.edu.cn/

% checking input
if ~exist('method', 'var')
    method = 'pearson';
end

nchan1 = size(input1, 1);
nchan2 = size(input2, 1);

nfreq  = size(input1, 2);

ntime1  = size(input1, 3);
ntime2  = size(input2, 3);
if ntime1 < ntime2
    ntime = ntime1;
    input2 = input2(:,:,1:ntime);
    fprintf('Cutting data to the same length: input1 length %d', ntime);
elseif ntime1 > ntime2
    ntime = ntime2;
    input1 = input1(:,:,1:ntime);
    fprintf('Cutting data to the same length: input2 length %d', ntime);
end
    
input1 = permute(input1,[3 2 1]);
input2 = permute(input2,[3 2 1]);

% initializing output matrices
timepoint = slidw/2 : slidw - overlap : ntime - slidw/2;
timnumber = length(timepoint);
corspectrm = ones(nchan1*nchan2, nfreq, length(timepoint));

parfor freq = 1:nfreq
    temp1 = squeeze(input1(:, freq, :));
    temp2 = squeeze(input2(:, freq, :));
    v = zeros(nchan1*nchan2, timnumber);
    timepoint = slidw/2 : slidw - overlap : ntime - slidw/2;

    for tim = 1:timnumber
        timewindow = timepoint(tim) - slidw/2 + 1 : timepoint(tim) + slidw/2;
        r = corr(temp1(timewindow,:), temp2(timewindow,:), 'Type', method);
        v(:, tim) = r(:); 
    end
    corspectrm(:, freq, :) = v;
end

% output
pcor = corspectrm;
