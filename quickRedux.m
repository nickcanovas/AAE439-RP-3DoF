clear all;
close all

testid = 1;
saveData =true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE THESE THINGS
dataFileNameMB = ['motor_1.tdms']; %ENSURE .tdms IS AT END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LFmatFilename = sprintf('Test_%d_Data',testid);

        LFMB = reduceTDMS(dataFileNameMB,1,nan);
 
   % Package and Save
    if saveData% && ZeroSave
        fprintf('Saving Low Frequency Data...\n')
        save([pwd,'\',LFmatFilename],'LFMB')
        fprintf('Low Frequency Data Saved.\n\n')
    end

allData.LFMB = LFMB;

t = LFMB.time.Value;
F = LFMB.cc_lc_02.Value;

plot(t, F);