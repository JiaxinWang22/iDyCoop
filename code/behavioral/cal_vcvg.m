% This code calculates VC and VG (differences) across teammates
% 
% iVelocity, iPosition are Nx1 cell arrays with N rows representing N sessions
% each cell contains a 10x4 cell matrix, with rows as trials and columns as players  

%% velocity and position theta angle
tempVelocity = iVelocity100;
tempPosition = iPosition100;
runN = size(iVelocity100,1);
tempCosTheta = cell(runN,1);
for iRun = 1:runN
    for iTeam = 1:2
        for iTrial = 1:10
            % vector pointing from player 1 to 2
            temp1to2 = tempPosition{iRun,1}{iTrial,2*iTeam} - tempPosition{iRun,1}{iTrial,2*iTeam-1};
            temp1to2 = temp1to2(1:end-1,:);  % because we are using velocity
            % if two dots are at the exact same place, temp1to2 will be [0 0], we give NaN to them, but not 0
            
            tempVel1 = tempVelocity{iRun,1}{iTrial,2*iTeam-1}(:,1:2);
            tempVel2 = tempVelocity{iRun,1}{iTrial,2*iTeam}(:,1:2);
            
            % if any one of the Velocity is zero, CosTheta will be NaN, it's ok
            tempCosTheta{iRun,1}{iTrial,2*iTeam-1} = dot(tempVel1,temp1to2,2)./(vecnorm(tempVel1,2,2).*vecnorm(temp1to2,2,2));
            tempCosTheta{iRun,1}{iTrial,2*iTeam} = dot(tempVel2,-temp1to2,2)./(vecnorm(tempVel2,2,2).*vecnorm(-temp1to2,2,2));
            
            % some error handling
            if any(abs(tempCosTheta{iRun,1}{iTrial,2*iTeam-1}) > 1)
                tempCosTheta{iRun,1}{iTrial,2*iTeam-1}(tempCosTheta{iRun,1}{iTrial,2*iTeam-1} > 1) = 1;
                tempCosTheta{iRun,1}{iTrial,2*iTeam-1}(tempCosTheta{iRun,1}{iTrial,2*iTeam-1} < -1) = -1;
            end
            if any(abs(tempCosTheta{iRun,1}{iTrial,2*iTeam}) > 1)
                tempCosTheta{iRun,1}{iTrial,2*iTeam}(tempCosTheta{iRun,1}{iTrial,2*iTeam} > 1) = 1;
                tempCosTheta{iRun,1}{iTrial,2*iTeam}(tempCosTheta{iRun,1}{iTrial,2*iTeam} < -1) = -1;
            end
        end
    end
end
% save cosTheta if needed

%% velocity pointing to finishing line
% extract directly from velocity
tempVelocity = iVelocity100;
tempVel_G = cell(runN,1);
for iRun = 1:runN
    for iTeam = 1:2
        for iTrial = 1:10
            tempVel_G{iRun,1}{iTrial,2*iTeam-1} = tempVelocity{iRun,1}{iTrial,2*iTeam-1}(:,1);
            tempVel_G{iRun,1}{iTrial,2*iTeam} = tempVelocity{iRun,1}{iTrial,2*iTeam}(:,1);
            
            % replace points with temp1to2 = 0 with real x-axis velocity
            temp1to2 = tempPosition{iRun,1}{iTrial,2*iTeam} - tempPosition{iRun,1}{iTrial,2*iTeam-1};
            temp1to2 = temp1to2(1:end-1,:);  % because we are using velocity
            tempOverlap = (temp1to2(:,1) == 0 & temp1to2(:,2) == 0);
            tempVel_G{iRun,1}{iTrial,2*iTeam-1}(tempOverlap) = tempVelocity{iRun,1}{iTrial,2*iTeam-1}(tempOverlap,1);
            tempVel_G{iRun,1}{iTrial,2*iTeam}(tempOverlap) = tempVelocity{iRun,1}{iTrial,2*iTeam}(tempOverlap,1);
        end
    end
end
iVel_G = tempVel_G;

%% velocity pointing to teammate
tempVelocity = iVelocity100;
tempVel_C = cell(runN,1);
for iRun = 1:runN
    for iTeam = 1:2
        for iTrial = 1:10
            tempVel_C{iRun,1}{iTrial,2*iTeam-1} = tempVelocity{iRun,1}{iTrial,2*iTeam-1}(:,3).*tempCosTheta{iRun,1}{iTrial,2*iTeam-1};
            tempVel_C{iRun,1}{iTrial,2*iTeam} = tempVelocity{iRun,1}{iTrial,2*iTeam}(:,3).*tempCosTheta{iRun,1}{iTrial,2*iTeam};
            
            % replace NaN values with 0, because NaN is caused by 0 velocity
            tempVel_C{iRun,1}{iTrial,2*iTeam-1}(isnan(tempVel_C{iRun,1}{iTrial,2*iTeam-1})) = 0;
            tempVel_C{iRun,1}{iTrial,2*iTeam}(isnan(tempVel_C{iRun,1}{iTrial,2*iTeam})) = 0;
            
            % replace points with temp1to2 = 0 with 0, assume no velocity to mate
            temp1to2 = tempPosition{iRun,1}{iTrial,2*iTeam} - tempPosition{iRun,1}{iTrial,2*iTeam-1};
            temp1to2 = temp1to2(1:end-1,:);  % because we are using velocity
            tempOverlap = (temp1to2(:,1) == 0 & temp1to2(:,2) == 0);
            tempVel_C{iRun,1}{iTrial,2*iTeam-1}(tempOverlap) = 0;
            tempVel_C{iRun,1}{iTrial,2*iTeam}(tempOverlap) = 0;          
        end
    end
end
iVel_C = tempVel_C;

%% sum and diff of velocity
tempVel_G = iVel_G;
tempVel_C = iVel_C;
tempVel_GPair = cell(runN,1);
tempVel_CPair = cell(runN,1);
for iRun = 1:runN
    for iTeam = 1:2
        for iTrial = 1:10
            tempVel_GPair{iRun,1}{iTrial,iTeam}(:,1) = sum([tempVel_G{iRun,1}{iTrial,2*iTeam-1}, tempVel_G{iRun,1}{iTrial,2*iTeam}],2);
            tempVel_GPair{iRun,1}{iTrial,iTeam}(:,2) = abs(tempVel_G{iRun,1}{iTrial,2*iTeam-1} - tempVel_G{iRun,1}{iTrial,2*iTeam});
            
            tempVel_CPair{iRun,1}{iTrial,iTeam}(:,1) = sum([tempVel_C{iRun,1}{iTrial,2*iTeam-1}, tempVel_C{iRun,1}{iTrial,2*iTeam}],2);
            tempVel_CPair{iRun,1}{iTrial,iTeam}(:,2) = abs(tempVel_C{iRun,1}{iTrial,2*iTeam-1} - tempVel_C{iRun,1}{iTrial,2*iTeam});
        end
    end
end
iVel_GPair = tempVel_GPair;
iVel_CPair = tempVel_CPair;
