% This code extracts M and I state epochs 
%
% distance: Nx1 cell array of Euclidean distance between teammates
% time_sr: Nx1 cell array of time stamps for each trial (seconds)
% reset_point: Nx1 cell array of time points where reset happened

runN = size(distance,1); % total session numbers

%% M state
M_crt = 2; % 2s time window for M state
temp_Mstate = cell(runN,1);
for iRun = 1:runN % session
    for iTeam = 1:2 % team
        for iTrial = 1:10 % trial

            % find line-up period, this doesn't include the reset points
            tempLineUp = find(distance{iRun,1}{iTrial,iTeam}(:,3) < 140);
            temp_lineup_timewd = find_epoch(tempLineUp, 'numeric', 'valuemat');

            if ~isempty(temp_lineup_timewd)
                temp_exclude = [];
                for pointer = 1:size(temp_lineup_timewd,1)
                    tempbeg = temp_lineup_timewd(pointer,1);
                    tempend = temp_lineup_timewd(pointer,2);

                    % exclude those epochs with less than 2s connection
                    if diff(time_sr{iRun,1}{iTrial,1}([tempbeg tempend])) < M_crt
                        temp_exclude = [temp_exclude; pointer];
                    end
                end
                temp_Mstate{iRun,1}{iTrial,iTeam}(temp_exclude,:) = [];
            end
        end
    end
end
timewin_Mstate = temp_Mstate;

%% I state epochs
temp_Mstate = timewin_Mstate;
temp_Iperiod = cell(runN,1);
temp_epoch_Istate = cell(runN,1);
for iRun = 1:runN
    for iTeam = 1:2
        for iTrial = 1:10

            if ~isempty(temp_Mstate{iRun,1}{iTrial,iTeam})
                % concatenate M states
                temp_Mperiod = [];
                for pointer = 1:size(temp_Mstate{iRun,1}{iTrial,iTeam},1)
                    tempmovebeg = temp_Mstate{iRun,1}{iTrial,iTeam}(pointer,1);
                    tempmoveend = temp_Mstate{iRun,1}{iTrial,iTeam}(pointer,2);
                    temp_Mperiod = [temp_Mperiod; (tempmovebeg:tempmoveend)'];
                end
                % time periods other than M states are I states 
                temp_Iperiod{iRun,1}{iTrial,iTeam} = setdiff(1:size(distance{iRun,1}{iTrial,1},1), temp_Mperiod)';
            else
                temp_Iperiod{iRun,1}{iTrial,iTeam} = (1:size(distance{iRun,1}{iTrial,1},1))';
            end

            % this includes reset points
            temp_epoch_Istate{iRun,1}{iTrial,iTeam} = find_epoch(temp_Iperiod{iRun,1}{iTrial,iTeam}, 'numeric', 'valuemat');
        end
    end
end
epoch_Istate = temp_epoch_Istate;

%% I state epochs with same length as M states for temporal comparison
M_crt = 2; % same as M crt for comparison
temp_epoch_Istate = epoch_Istate;
for iRun = 1:runN
    for iTeam = 1:2
        for iTrial = 1:10
            
            temp_exclude = [];
            for pointer = 1:size(temp_epoch_Istate{iRun,1}{iTrial,iTeam},1)
                tempbeg = temp_epoch_Istate{iRun,1}{iTrial,iTeam}(pointer,1);
                tempend = temp_epoch_Istate{iRun,1}{iTrial,iTeam}(pointer,2);

                % do not include the reset points
                if ismember(tempbeg, reset_point{iRun,1}{iTrial,iTeam}) && (tempbeg < tempend)
                    tempbeg = tempbeg + 1;
                    temp_epoch_Istate{iRun,1}{iTrial,iTeam}(pointer,1) = temp_epoch_Istate{iRun,1}{iTrial,iTeam}(pointer,1) + 1;
                end

                if diff(time_sr{iRun,1}{iTrial,1}([tempbeg tempend])) < M_crt
                    temp_exclude = [temp_exclude; pointer];
                end
            end
            temp_epoch_Istate{iRun,1}{iTrial,iTeam}(temp_exclude,:) = [];
        end
    end
end
timewin_Istate = temp_epoch_Istate;
