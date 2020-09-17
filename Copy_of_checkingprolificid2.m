clear, clc

% % % load table with prolific id
load('Masterlabvancedtable.mat')

% % % stimulus sets
load(strcat('BKKmirrorss 2.mat')); % BKKmirrors
load(strcat('BKKoriginals 2.mat')); % BKKorigs

% % % hard-coding tempvec
temp_lims = [100 300]; 

S = struct();
cnt = 1

for b = 1:size(Masterlabvancedtable1,1)
    
    %get rid of the quotes
    bb = split(Masterlabvancedtable1.crowdsourcing_subj_id{b},'""');
    
    bbb = bb{1};
    
    %S.prolificid
    
    c = strfind(Masterlabvancedtable1.crowdsourcing_subj_id,bbb);
    
    %x = [];
    x = b;
    
    for d = b+1:length(c)
        if ~ isempty(c{d})
            x = [x d];
        end
    end
    
    if length(x) == 1
        
    else % party here 
        megax(cnt,:) = x;
        S(cnt).prolificid=char(bbb);
        S(cnt).prolificidstr=convertCharsToStrings(bbb); %make a string for easy labeling
        S(cnt).stims1 = getNRKstimorder(megax(cnt,1));
        S(cnt).stims2 = getNRKstimorder(megax(cnt,2));
        
       save('test.mat','S')
    
    % % % session 1 
    datadir_1 = ['/Users/nicoleahuizinga/Desktop/TEMPITAPS_april/data' filesep bbb filesep 'session_1/'];
    for t = 1:1
        wavname = strcat(datadir_1,'blockNr_3_taskNr_2_trialNr_',num2str(t-1),'_tapping.wav');
        [audio, sf] = audioread(wavname);
        
        if [S(cnt).stims1(t,2)] == 1
           BKK_curr_trial = BKKorigs{t};
        else
           BKK_curr_trial = BKKmirrors{t};
        end
        
        [taptimes_raw, taptimes_zeroed, rhythm_onsets, trial_thresh] = mollyread (audio, sf, BKK_curr_trial);
        S(cnt).taptimes_raw_sesh1{t} = taptimes_raw;
        S(cnt).taptimes_zeroed_sesh1{t} = taptimes_zeroed;
        S(cnt).stimulus_sesh1{t} = rhythm_onsets; 
        S(cnt).trialthresh_sesh1{t} = trial_thresh;
        
    end
    
    for m = 1:2
        % get SMT trials 1&2
        wavcomfy = strcat(datadir_1,'blockNr_1_taskNr_2_trialNr_0_SMT_comfortable',num2str(m),'.wav');
        [audio, sf] = audioread(wavcomfy);
        
        [ctap_times_raw, ctrial_thresh, cints, cmean, cstd] = molly2read(audio, sf);
        S(cnt).comfytaptimes_raw_sesh1_{m} = ctap_times_raw; % raw comfy tap times
        S(cnt).comfytaptimes_thresh_sesh1_{m} = ctrial_thresh; % comfy trial thresholds
        S(cnt).comfyints_sesh1_{m} = cints;
        S(cnt).comfyintsmean_sesh1_{m} = cmean;
        S(cnt).comfyintsstd_sesh1_{m} = cstd;
        
        
        % get slow forced tapping trials 1&2
        wavslow = strcat(datadir_1,'blockNr_1_taskNr_2_trialNr_0_SMT_slow',num2str(m),'.wav');
        [audio, sf] = audioread(wavslow);
        
        [ctap_times_raw, ctrial_thresh, cints, cmean, cstd] = molly2read(audio, sf);
        S(cnt).slowtaptimes_raw_sesh1_{m} = ctap_times_raw;
        S(cnt).slowtaptimes_thresh_sesh1_{m} = ctrial_thresh;
        S(cnt).slowints_sesh1_{m} = cints;
        S(cnt).slowintsmean_sesh1_{m} = cmean;
        S(cnt).slowintsstd_sesh1_{m} = cstd;
        
        % get fast forced tapping trials 1&2
        wavfast = strcat(datadir_1,'blockNr_1_taskNr_2_trialNr_0_SMT_fast',num2str(m),'.wav');
        [audio, sf] = audioread(wavfast);
        
        [ctap_times_raw, ctrial_thresh, cints, cmean, cstd] = molly2read(audio, sf);
        S(cnt).fasttaptimes_raw_sesh1_{m} = ctap_times_raw;
        S(cnt).fasttaptimes_thresh_sesh1_{m} = ctrial_thresh;
        S(cnt).fastints_sesh1_{m} = cints;
        S(cnt).fastintsmean_sesh1_{m} = cmean;
        S(cnt).fastintsstd_sesh1_{m} = cstd;
        
    end
    
    for m = 3:3
        
        % % % load SMT trial 3 wav and retrieve tap times
        wavcomfy = strcat(datadir_1,'blockNr_4_taskNr_1_trialNr_0_SMT_comfortable3.wav');
        [audio, sf] = audioread(wavcomfy);
        
        [ctap_times_raw, ctrial_thresh, cints, cmean, cstd] = molly2read(audio, sf);
        S(cnt).comfytaptimes_raw_sesh1_{m} = ctap_times_raw;
        S(cnt).comfytaptimes_thresh_sesh1_{m} = ctrial_thresh;
        S(cnt).comfyints_sesh1_{m} = cints;
        S(cnt).comfyintsmean_sesh1_{m} = cmean;
        S(cnt).comfyintsstd_sesh1_{m} = cstd;
        S(cnt).comfyintsmean_cell2mat = cell2mat(S(cnt).comfyintsmean_sesh1_);
        %S(cnt).comfyintsmean_cell2matfinal = 
        S(cnt).comfyintsmean_final_sesh1 = mean(S(cnt).comfyintsmean_cell2mat);
        S(cnt).comfyintsstd_final_sesh1 = std(S(cnt).comfyintsmean_cell2mat);
        
         % get slow forced tapping trials 3
        wavslow = strcat(datadir_1,'blockNr_4_taskNr_1_trialNr_0_SMT_slow3.wav');
        [audio, sf] = audioread(wavslow);
        
        [ctap_times_raw, ctrial_thresh, cints, cmean, cstd] = molly2read(audio, sf);
        S(cnt).slowtaptimes_raw_sesh1_{m} = ctap_times_raw;
        S(cnt).slowtaptimes_thresh_sesh1_{m} = ctrial_thresh;
        S(cnt).slowints_sesh1_{m} = cints;
        S(cnt).slowintsmean_sesh1_{m} = cmean;
        S(cnt).slowintsstd_sesh1_{m} = cstd;
        S(cnt).slowintsmean_cell2mat = cell2mat(S(cnt).slowintsmean_sesh1_);
        S(cnt).slowintsmean_final_sesh1 = mean(S(cnt).slowintsmean_cell2mat);
        S(cnt).slowintsstd_final_sesh1 = std(S(cnt).slowintsmean_cell2mat);
        
        
        % get fast forced tapping trials 3
        wavfast = strcat(datadir_1,'blockNr_4_taskNr_1_trialNr_0_SMT_fast3.wav');
        [audio, sf] = audioread(wavfast);
        
        [ctap_times_raw, ctrial_thresh, cints, cmean, cstd] = molly2read(audio, sf);
        S(cnt).fasttaptimes_raw_sesh1_{m} = ctap_times_raw;
        S(cnt).fasttaptimes_thresh_sesh1_{m} = ctrial_thresh;
        S(cnt).fastints_sesh1_{m} = cints;
        S(cnt).fastintsmean_sesh1_{m} = cmean;
        S(cnt).fastintsstd_sesh1_{m} = cstd;
        S(cnt).fastintsmean_cell2mat = cell2mat(S(cnt).fastintsmean_sesh1_);
        S(cnt).fastintsmean_final_sesh1 = mean(S(cnt).fastintsmean_cell2mat);
        S(cnt).fastintsstd_final_sesh1 = std(S(cnt).fastintsmean_cell2mat);
        
        % get difference between mean 
        
    end
   
    % % % session 2 
    datadir_2 = ['/Users/nicoleahuizinga/Desktop/TEMPITAPS_april/data' filesep bbb filesep 'session_2/'];
    
    for t = 1:1 %specify the trials
        wavname = strcat(datadir_2,'blockNr_3_taskNr_2_trialNr_',num2str(t-1),'_tapping.wav');
        [audio, sf] = audioread(wavname);
        
        if [S(cnt).stims1(t,2)] == 1
           BKK_curr_trial = BKKorigs{t};
        else
           BKK_curr_trial = BKKmirrors{t};
        end
        
        [taptimes_raw, taptimes_zeroed, rhythm_onsets, trial_thresh] = mollyread (audio, sf, BKK_curr_trial);
        S(cnt).taptimes_raw_sesh2{t} = taptimes_raw;
        S(cnt).taptimes_zeroed_sesh2{t} = taptimes_zeroed;
        S(cnt).stimulus_sesh2{t} = rhythm_onsets; 
        S(cnt).trialthresh_sesh2{t} = trial_thresh;
    end
    
    
    
    cnt = cnt+1; % done with this guy
    end
end

%% analysis part  ideas 

% tappies = S(1).taptimes_raw_sesh1{1,1};
% stimz = S(1).stimulus_sesh1{1,1};

%tempvec for mirror vs original



%% fun
function [stimorder] = getNRKstimorder (x)
matfilename = sprintf('/Users/nicoleahuizinga/Desktop/TEMPITAPS_april/NRKstimorders/NKRstimorder_%d.mat', x);
load (matfilename);
stimorder = NKstimorder;
end

function [taptimes_raw, taptimes_zeroed, rhythm_onsets, trial_thresh] = mollyread (audio, sf, BKK_curr_trial) % todo: find BKK & put here 

y = BKK_curr_trial;
data = audio;

t  = 1/sf:1/sf:length(data)/sf;

% get the stimulus ------------------------------------------------
% y = BKKmirrors{kk};
y = y(1:end-sf); % this chops off exactly 1 s, which is the noise
tstim = 1/sf:1/sf:length(y)/sf;
%         figure,  plot(tstim,y)

% % get rid of babies (noise floor) -----------------------------------------
bthresh = .01; % this threshold is manual for now
y(y < bthresh) = 0; % set everything below threshold to 0
%             figure, plot(tstim,y) % plot the nice clean signal

% first derivative --------------------------------------------------------
deriv = diff(y); % get the first derivative
Td = tstim(1:end-1); % get a new time vector
%         figure, plot(Td,deriv) % plot the derivative

% find onsets based on derivative % --------------------------------------
Thresh = bthresh; % why not just keep the threshold chosen above?
all_onsets = find(deriv > Thresh); % find all onsets (this will be wayyyyyy too many)
tones = []; % = taps
distance = diff(all_onsets); % what is the distance between each "onset"?
for ii = 1:length(all_onsets)-1 % loop through all of them
    if ii == 1 % always keep the first one
        tones = [tones all_onsets(ii)];
    else
        if distance(ii-1) > .05*sf % this is another manual threshold, here I don't take onsets that are closer than 100 ms to each other (people can't tap that fast)
            tones = [tones all_onsets(ii)];
        end
    end
end

% make 1s & 0s out of the red circles (Td(tones))
Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
Tdlgth(tones) = 1;
%         figure, stem(Td,Tdlgth)

rhythm_onsets = Td(tones);

rhythm_onsets_zeroes = rhythm_onsets - rhythm_onsets(1);

nints = length(diff(rhythm_onsets));

%% x (wavfile  data)

x = data(:,1);

% this is a new little interactive module for you to set the threshold
% on a trial-by-trial basis. Not fully automated, but neither is my
% brain right now.
ok = 0;
while ok == 0
    %titlename = S.prolificid;
    figure,  plot(t,x)
    %title(sprintf('Trial %s [titlename]);
    
    % % get rid of babies (noise floor) -----------------------------------------
    bthresh = input('Enter threshold to try:  '); % this threshold is manual for now
    x(x < bthresh) = 0; % set everything below threshold to 0
    %     figure, plot(t,x) % plot the nice clean signal
    
    % first derivative --------------------------------------------------------
    deriv = diff(x); % get the first derivative
    Td = t(1:end-1); % get a new time vector
    figure, plot(Td,deriv) % plot the derivative
    
    check = input('Is it ok (y=1|n=0)?  ');
  
    trial_thresh = bthresh; % save the single-trial threshold
    
    close all
    ok = check;
end

% trial_thresh(kk) = bthresh; % save the single-trial threshold

% find onsets based on derivative % --------------------------------------
Thresh = bthresh; % why not just keep the threshold chosen above?
all_onsets = find(deriv > Thresh); % find all onsets (this will be wayyyyyy too many)
tones = []; % = taps
distance = diff(all_onsets); % what is the distance between each "onset"?
for ii = 1:length(all_onsets)-1 % loop through all of them
    if ii == 1 % always keep the first one
        tones = [tones all_onsets(ii)];
    else
        if distance(ii-1) > .1*sf % this is another manual threshold, here I don't take onsets that are closer than 100 ms to each other (people can't tap that fast)
            tones = [tones all_onsets(ii)];
        end
    end
end


% make 1s & 0s out of the red circles (Td(tones))
Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
Tdlgth(tones) = 1;

taptimes_raw = Td(tones);

taptimes_zeroed = taptimes_raw - rhythm_onsets(1); % now have to assume rhythm play is perfect until new info

end

function [ctaptimes_raw, ctrial_thresh, cints, cmean, cstd] = molly2read (audio, sf) % todo: find BKK & put here 

%y = BKK_curr_trial;
data = audio;

t  = 1/sf:1/sf:length(data)/sf;

% % get the stimulus ------------------------------------------------
% % y = BKKmirrors{kk};
% y = y(1:end-sf); % this chops off exactly 1 s, which is the noise
% tstim = 1/sf:1/sf:length(y)/sf;
% %         figure,  plot(tstim,y)

% % % get rid of babies (noise floor) -----------------------------------------
% bthresh = .01; % this threshold is manual for now
% y(y < bthresh) = 0; % set everything below threshold to 0
% %             figure, plot(tstim,y) % plot the nice clean signal

% % first derivative --------------------------------------------------------
% deriv = diff(y); % get the first derivative
% Td = tstim(1:end-1); % get a new time vector
% %         figure, plot(Td,deriv) % plot the derivative

% % find onsets based on derivative % --------------------------------------
% Thresh = bthresh; % why not just keep the threshold chosen above?
% all_onsets = find(deriv > Thresh); % find all onsets (this will be wayyyyyy too many)
% tones = []; % = taps
% distance = diff(all_onsets); % what is the distance between each "onset"?
% for ii = 1:length(all_onsets)-1 % loop through all of them
%     if ii == 1 % always keep the first one
%         tones = [tones all_onsets(ii)];
%     else
%         if distance(ii-1) > .05*sf % this is another manual threshold, here I don't take onsets that are closer than 100 ms to each other (people can't tap that fast)
%             tones = [tones all_onsets(ii)];
%         end
%     end
% end

% % make 1s & 0s out of the red circles (Td(tones))
% Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
% Tdlgth(tones) = 1;
% %         figure, stem(Td,Tdlgth)

% rhythm_onsets = Td(tones);
% 
% rhythm_onsets_zeroes = rhythm_onsets - rhythm_onsets(1);
% 
% nints = length(diff(rhythm_onsets));

%% z (wavfile data)

z = data(:,1);

% this is a new little interactive module for you to set the threshold
% on a trial-by-trial basis. Not fully automated, but neither is my
% brain right now.
ok = 0;
while ok == 0
    figure,  plot(t,z)
    
    % % get rid of babies (noise floor) -----------------------------------------
    bthresh = input('Enter threshold to try:  '); % this threshold is manual for now
    z(z < bthresh) = 0; % set everything below threshold to 0
    %     figure, plot(t,x) % plot the nice clean signal
    
    % first derivative --------------------------------------------------------
    deriv = diff(z); % get the first derivative
    Td = t(1:end-1); % get a new time vector
    figure, plot(Td,deriv) % plot the derivative
    
    check = input('Is it ok (y=1|n=0)?  ');
  
    ctrial_thresh = bthresh; % save the single-trial threshold
    
    close all
    ok = check;
end

% trial_thresh(kk) = bthresh; % save the single-trial threshold

% find onsets based on derivative % --------------------------------------
Thresh = bthresh; % why not just keep the threshold chosen above?
all_onsets = find(deriv > Thresh); % find all onsets (this will be wayyyyyy too many)
tones = []; % = taps
distance = diff(all_onsets); % what is the distance between each "onset"?
for ii = 1:length(all_onsets)-1 % loop through all of them
    if ii == 1 % always keep the first one
        tones = [tones all_onsets(ii)];
    else
        if distance(ii-1) > .1*sf % this is another manual threshold, here I don't take onsets that are closer than 100 ms to each other (people can't tap that fast)
            tones = [tones all_onsets(ii)];
        end
    end
end


% make 1s & 0s out of the red circles (Td(tones))
Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
Tdlgth(tones) = 1;

ctaptimes_raw = Td(tones);
cints = diff(ctaptimes_raw);
cmean = mean(cints);
cstd = std(cints);

%taptimes_zeroed = taptimes_raw - rhythm_onsets(1); % now have to assume rhythm play is perfect until new info

end
%stimulusordersession1 = sprintf('/NRKstimorders/NRKstimorder_