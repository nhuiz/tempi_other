clear; close all; clc
 
 % % % hard-coding tempvec
 temp_lims = [100 300]; 

 % directories where to get stuff.
 % datadir = '/Users/molly/Downloads/laura pilot data/';
 % datadir = '/Users/molly/Downloads/TEMPITAPS_april/ece pilot/';
 % stimdir = '/Users/molly/Downloads/TEMPITAPS_april/stimorders/';

%  datadir = '/Users/molly.henry/Downloads/Molly/';
%  stimdir = '/Users/molly.henry/Downloads/';
 datadir = '/Users/nicoleahuizinga/Desktop/TEMPITAPS_april/data_old/Ece/';
 stimdir = '/Users/nicoleahuizinga/Desktop/TEMPITAPS_april/';

 subnr = 2; % populate this

 trial{1} = [1 2 3];
 trial{2} = [1 2 3];%1:40; % can use this structure to throw out some trials in case they are trash

 % for 
     s = subnr;
     [trial_thresh,num_outliers,time_outliers] = deal([]);

     stimorder = [zeros(length(trial{s}),1), 2*ones(length(trial{s}),1)]; % TEMPORARY
     load(strcat(stimdir,'BKKmirrorss.mat'))
 %     load(strcat(stimdir,'Rstimorder_',num2str(s),'.mat'))
 %     load(strcat(stimdir,'Rstimorder_3.mat'))

     for kk = 3;%%trial{s}
         file = strcat(datadir,'main_pt',num2str(s),'_trial',num2str(kk),'.wav');
 %         file = strcat(datadir,'main_pt42_trial',num2str(kk),'.wav');
         [data,sf] = audioread(file);
         t  = 1/sf:1/sf:length(data)/sf;         

         % get the stimulus ------------------------------------------------
         y = BKKmirrors{kk};
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
 % %             figure, plot(Td,deriv), hold on % plot the derivative
 % %             plot(Td(tones),0,'or') % with the onsets on top

         % make 1s & 0s out of the red circles (Td(tones))
         Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
         Tdlgth(tones) = 1;
 %         figure, stem(Td,Tdlgth)

         rhythm_onsets = Td(tones);
 %         aud_delay = rhythm_onsets(1);
 %         rhythm_onsets_zeroed = rhythm_onsets - aud_delay;
         rhythm_onsets_zeroes = rhythm_onsets - rhythm_onsets(1);

         nints = length(diff(rhythm_onsets));
 %         if stimorder(kk,2) == 1
 %             tempvec = temp_lims(1) * (temp_lims(2)/temp_lims(1)).^((0:nints-1)/(nints-1)); % can make only for individual rhythm, need to know stimulus identity.
 %         elseif stimorder(kk,2) == 2
 %             tempvec = temp_lims(2) * (temp_lims(1)/temp_lims(2)).^((0:nints-1)/(nints-1)); % can make only for individual rhythm, need to know stimulus identity.
 %         end  


         % get the tapping data (this plot if for inspection / thresholding)
         x = data(:,1);
 %         figure,  plot(t,x)

         % this is a new little interactive module for you to set the threshold
         % on a trial-by-trial basis. Not fully automated, but neither is my
         % brain right now.
         ok = 0;
         while ok == 0
             x = data(:,1);
             figure,  plot(t,x)

             % % get rid of babies (noise floor) -----------------------------------------
             bthresh = input('Enter threshold to try:  '); % this threshold is manual for now
             x(x < bthresh) = 0; % set everything below threshold to 0
 %                 figure, plot(t,x) % plot the nice clean signal

             % first derivative --------------------------------------------------------
             deriv = diff(x); % get the first derivative
             Td = t(1:end-1); % get a new time vector
             figure, plot(Td,deriv) % plot the derivative

             check = input('Is it ok (y=1|n=0)?  ');
             close all
             ok = check;
         end
         trial_thresh(kk) = bthresh; % save the single-trial threshold 

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
         %     figure, plot(Td,deriv), hold on % plot the derivative
         %     plot(Td(tones),0,'or') % with the onsets on top

         % make 1s & 0s out of the red circles (Td(tones))
         Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
         Tdlgth(tones) = 1;
         %     figure, stem(Td,Tdlgth) % plot the TAPS ONLY

         %% now how to analyze this?
         taptimes_raw = Td(tones);
 %         taptimes_raw = taptimes_raw - aud_delay; % we no longer know this
 %         taptimes_zeroed = taptimes_raw - taptimes_raw(1);
         taptimes_zeroed = taptimes_raw - rhythm_onsets(1); % now have to assume rhythm play is perfect until new info

         % need the specific rhythm info here (will need stimorder specific to
         % subject to decode)
         % get stimulus times raw and temp vec per subject

         % get the intervals between taps
         ints = diff(taptimes_zeroed);
         idxs = 1:length(ints);
 %         figure, subplot(1,2,1), plot(ints,'ok-'), hold on
 %         figure, plot(taptimes_zeroed(1:end-1),ints,'ok-'), hold on

 %         scatter(1:length(ints),ints,'k')%,lsline

         % plot the tempvec with the produced tap series.
         if stimorder(kk,2) == 1
             oversampled_tempvec = temp_lims(1) * (temp_lims(2)/temp_lims(1)).^((0:10000-1)/(10000-1)); % can make only for individual rhythm, need to know stimulus identity.
         elseif stimorder(kk,2) == 2
             oversampled_tempvec = temp_lims(2) * (temp_lims(1)/temp_lims(2)).^((0:10000-1)/(10000-1)); % can make only for individual rhythm, need to know stimulus identity.
         end
         oversampled_tempvec = oversampled_tempvec / 1000;
         oversampled_timevec = linspace(0,rhythm_onsets(end),10000);
         quart_temp = oversampled_tempvec;
         half_temp  = 2*oversampled_tempvec;
         whole_temp = 4*oversampled_tempvec;

 %         plot(oversampled_timevec,quart_temp,'m')
 %         plot(oversampled_timevec,half_temp,'c')
 %         plot(oversampled_timevec,whole_temp,'g')


         tap_temp = [];
         for tt = 1:length(taptimes_raw)
             %tap_temp = [tap_temp oversampled_tempvec(nearest_mjh(oversampled_timevec,taptimes_raw(tt)))];
             tap_temp = [tap_temp oversampled_tempvec(closeby(oversampled_timevec,taptimes_raw(tt)))];
         end

         % get the temp vecs for each metrical level
         quart_tap = tap_temp(1:end-1);
         half_tap  = 2*tap_temp(1:end-1);
         whole_tap = 4*tap_temp(1:end-1);

         % plot all the metrical levels on the produced interval time series
 %         plot(1:length(ints),quart_temp,'m')
 %         plot(1:length(ints),half_temp,'c')
 %         plot(1:length(ints),whole_temp,'g')
         %         figure, scatter(tap_temp(1:end-1),ints)


         % this is for removing outliers that are crazy, this threshold can
         % be adjusted. 
         out_perc = .5; % percentage of intended interval to be within, if outside of this, outlier
         % this is oversampled for plotting
         quart_range = [quart_temp + (quart_temp * out_perc); quart_temp - (quart_temp * out_perc)];
         half_range  = [half_temp + (half_temp * out_perc); half_temp - (half_temp * out_perc)];
         whole_range = [whole_temp + (whole_temp * out_perc); whole_temp - (whole_temp * out_perc)];

         % and this is for individual taps for calculations. 
         quart_tap_range = [quart_tap + (quart_tap * out_perc); quart_tap - (quart_tap * out_perc)];
         half_tap_range  = [half_tap + (half_tap * out_perc); half_tap - (half_tap * out_perc)];
         whole_tap_range = [whole_tap + (whole_tap * out_perc); whole_tap - (whole_tap * out_perc)];

 %         plot(oversampled_timevec,quart_range(1,:),'m--',oversampled_timevec,quart_range(2,:),'m--')
 %         plot(oversampled_timevec,half_range(1,:),'c--',oversampled_timevec,half_range(2,:),'c--')
 %         plot(oversampled_timevec,whole_range(1,:),'g--',oversampled_timevec,whole_range(2,:),'g--')

         % now remove outliers ---------------------------------------------
 %         clean_trialdata = ints;        
 %         clidxs = idxs;

         % this is to remove some at the beginning in case they aren't
         % stable yet
         ok = 0;
         while ok ~= 1
             figure, plot(taptimes_zeroed(1:end-1),ints,'ok-'), hold on
             plot(taptimes_zeroed(1:end-1),quart_tap,'m')
             plot(taptimes_zeroed(1:end-1),half_tap,'c')
             plot(taptimes_zeroed(1:end-1),whole_tap,'g')
             plot(taptimes_zeroed(1:end-1),quart_tap_range(1,:),'m--',taptimes_zeroed(1:end-1),quart_tap_range(2,:),'m--')
             plot(taptimes_zeroed(1:end-1),half_tap_range(1,:),'c--',taptimes_zeroed(1:end-1),half_tap_range(2,:),'c--')
             plot(taptimes_zeroed(1:end-1),whole_tap_range(1,:),'g--',taptimes_zeroed(1:end-1),whole_tap_range(2,:),'g--')
             how_many_unstable = input('How many unstable taps to remove?  ');
             plot(taptimes_zeroed(1:how_many_unstable),ints(1:how_many_unstable),'rx','MarkerSize',10)
             check = input('Is this ok?  ');
             close all
             ok = check;
         end
 %         clean_trialdata = ints(how_many_unstable+1:end);
 %         clidxs = idxs(how_many_unstable+1:end);

         % first need to correct the metrical level assignment in case there
         % is an accidental one (then it will be compared to the wrong
         % outlier criteria                
         % first choose which level
         met_lvl = [];
         met_lvl(1:how_many_unstable) = nan;
         for mla = how_many_unstable+1:length(ints)
 %             d = abs([ints(mla)-quart_tap(mla) ints(mla)-half_tap(mla) ints(mla)-whole_tap(mla)]);

             % revised to make metrical-lvl assignment based on relative error instead of abs
             % there are 2 ways to do this. 
 %             d = abs([(ints(mla)-quart_tap(mla))./quart_tap(mla) (ints(mla)-half_tap(mla))./half_tap(mla) (ints(mla)-whole_tap(mla))./whole_tap(mla)]); % relativize re: true tempo
             d = abs([(ints(mla)-quart_tap(mla))./ints(mla) (ints(mla)-half_tap(mla))./ints(mla) (ints(mla)-whole_tap(mla))./ints(mla)]); % relativize re: produced interval
             % I never know the right thing to relativize wrt, esp. when tempo and metrical level are potentially changing. Will have to test both ways. 
             met_lvl = [met_lvl find(d == min(d))];
 %             if met_lvl(end) == 1
 %                 plot(taptimes_zeroed(mla),ints(mla),'om')
 %             elseif met_lvl(end) == 2
 %                 plot(taptimes_zeroed(mla),ints(mla),'oc')
 %             elseif met_lvl(end) == 3
 %                 plot(taptimes_zeroed(mla),ints(mla),'og')
 %             end
         end


 %         [outliers, clean_trialdata] = deal([]);
         [outliers, clean_trialdata, clidxs] = deal([]);
         for ay = how_many_unstable+1:length(ints)
             if met_lvl(ay) == 1
                 if ints(ay) <= quart_tap_range(1,ay) && ints(ay) >= quart_tap_range(2,ay)
                     clean_trialdata = [clean_trialdata, ints(ay)];
                     clidxs = [clidxs ay];
                 else
 %                     clean_trialdata(ay) = nan;
                     outliers = [outliers; [ay ints(ay)]];
                 end
             elseif met_lvl(ay) == 2
                 if ints(ay) <= half_tap_range(1,ay) && ints(ay) >= half_tap_range(2,ay)
                     clean_trialdata = [clean_trialdata, ints(ay)];
                     clidxs = [clidxs ay];
                 else
 %                     clean_trialdata(ay) = nan;
                     outliers = [outliers; [ay ints(ay)]];
                 end
             elseif met_lvl(ay) == 3
                 if ints(ay) <= whole_tap_range(1,ay) && ints(ay) >= whole_tap_range(2,ay)
                     clean_trialdata = [clean_trialdata, ints(ay)];
                     clidxs = [clidxs ay];
                 else
 %                     clean_trialdata(ay) = nan;
                     outliers = [outliers; [ay ints(ay)]];
                 end
             end
         end

         % start from either side
 %         met_lvl(outliers(:,1)') = nan;
 %         first = met_lvl(how_many_unstable+1);
 %         final = met_lvl(end);
         met_lvl = met_lvl(clidxs);
         first = met_lvl(1);
         final = met_lvl(end);

 %         [preswitch,postswitch] = deal(zeros(size(met_lvl)));
         preswitch = met_lvl == first;
         postswitch = met_lvl == final;
 %         preswitch(outliers(:,1)') = nan;
 %         postswitch(outliers(:,1)') = nan;


         % so far I'm totally not able to think of a non-thresholded way to
         % handle the problem of identifying a switch point or region. 
         nmetsno = 3; % this may need to be adjusted. 
         [tmpout,tmpswitch] = deal([]);

         % coming from the "left"
         for nm = 1:length(preswitch)-nmetsno
             if preswitch(nm)==0
                 lvec = preswitch(nm+1:nm+nmetsno);
                 if sum(lvec==0)==nmetsno
                     tmpswitch = [tmpswitch, nm];
                     break
                 else
                     tmpout = [tmpout, nm];           

                 end
             end
         end

         % coming from the "right"
         for nm = length(postswitch):-1:nmetsno+1
             if postswitch(nm)==0
                 lvec = postswitch(nm-1:-1:nm-nmetsno);
                 if sum(lvec==0)==nmetsno
                     tmpswitch = [tmpswitch nm];
                     break
                 else
                     tmpout = [tmpout nm];
                 end
             end
         end

         % get the outliers out of the clean data. 
         allidx = ones(1,length(clean_trialdata));
         allidx(tmpout) = 0;
         clridxs = clidxs(logical(allidx));
         cleaner_trialdata = clean_trialdata(logical(allidx));        

         % combine all outlier info 
         % there might not be any
         if ~isempty(outliers)            
             outs_outs_baby = outliers(:,1)';
         else
             outs_outs_baby = [];
         end
         outsidx = sort([clidxs(tmpout) outs_outs_baby]);

         newtime = taptimes_zeroed(1:end-1);
         newints = ints;
         newtime(outsidx) = nan;
         newints(outsidx) = nan;
         newints(1:how_many_unstable) = nan;

         % save switch time
         % there might not be one. 
         if ~isempty(tmpswitch)
             switch_times_sym = taptimes_zeroed(clidxs(tmpswitch));
             switch_time = max(switch_times_sym); % SAVE THIS ONE.
             switch_temp = newints(min(clidxs(tmpswitch))); % the last interval size (tempo) produced before switch
         else
             switch_time = [];
             switch_temp = [];
         end
         total_range = max(newints) - min(newints);


         % combine switch time with other info (tempo dir, rhythm ID, etc)

         % analyze stability time course?


         % plot the full analysis summary
         figure, plot(taptimes_zeroed(1:end-1),ints,'ok-'), hold on

         plot(oversampled_timevec,quart_temp,'m')
         plot(oversampled_timevec,half_temp,'c')
         plot(oversampled_timevec,whole_temp,'g')

         plot(oversampled_timevec,quart_range(1,:),'m--',oversampled_timevec,quart_range(2,:),'m--')
         plot(oversampled_timevec,half_range(1,:),'c--',oversampled_timevec,half_range(2,:),'c--')
         plot(oversampled_timevec,whole_range(1,:),'g--',oversampled_timevec,whole_range(2,:),'g--')

         plot(taptimes_zeroed(1:how_many_unstable),ints(1:how_many_unstable),'rx','MarkerSize',10)

         for pi = 1:length(met_lvl)
             if met_lvl(pi) == 1
                 plot(taptimes_zeroed(clidxs(pi)),ints(clidxs(pi)),'om')
             elseif met_lvl(pi) == 2
                 plot(taptimes_zeroed(clidxs(pi)),ints(clidxs(pi)),'oc')
             elseif met_lvl(pi) == 3
                 plot(taptimes_zeroed(clidxs(pi)),ints(clidxs(pi)),'og')
             end
         end        

         plot(taptimes_zeroed(outsidx),ints(outsidx),'rx','MarkerSize',10)
         plot(taptimes_zeroed(clidxs(tmpswitch)),ints(clidxs(tmpswitch)),'*k','MarkerSize',20)

         xlabel('Time (s)'), ylabel('Produced interval (s)')
         set(gca,'XLim',[0 taptimes_zeroed(end)])


     end

 % end
