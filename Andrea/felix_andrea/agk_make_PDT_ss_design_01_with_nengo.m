% agk_make_PDT_ss_design
% for one subject
% makes the ss design of the PDT task based on the Tom et al. (2007) paper
% pictures here are put in as categorical vectors
% to allow a cue reactivity contrast (between categories modeling)

% INPUT
% cur_sub       : name of current subject
% which_folders : name of folder with niftis of task (pciks first) ???
% aggr          : aggregation level gambling matrix (we usually take 3)
% run_it        : should the job run? (will be only save if not)
% acc_rec       : should accept be included as a factor in model (def: no)
% expl_mask     : usually use a gray matter mask here
% cur_tmpl      : the ss job template to be used here, cell, first field is
%                 the one without accept reject as factor and second with
% ow            : overwrite exisiting results files? (0,1)

% OUTPUT
% error_message: success message or explanation what hasn't worked

% IMPLICIT OUTPUT
% writes the ss model into the results_ssdir in the subject folder

% THIS MODEL ('PDT_ss_design_01')
% 4 task on: PIC, PIC+GAM, PIC+GAM+OPT, FEEDBACK
% 3 param modulators in PIC: category dummy coded (gam, neg, neu)
% 6 param mod in all others: gain, abs(loss), ed; cat dummy coded

% DETAILS
% preproc_job is needed to get the slice timing parameters
% (to create time bins), will be looked for in the nifti folder provided
% will look for "preprocessing..." pattern
% abs values of losses are used
% nifti data need to be SPM12 preprocessed
% takes in an aggregation factor: e.g. 1 means no change; 2 means, the 12
% by 12 matrix will be aggregated to a 6-by-6 matrix, and so on
% the picture category is modeled by three dummy variables (i.e. all 0 is
% neutral)

% AUTHORSHIP
% author: Alexander Genauck
% date  : 08.08.2016
% email : alexander.genauck@charite.de

% CAVE: at input you should give as variable the text file with the
% amygdala time series (aus nengo)

function error_message = agk_make_PDT_ss_design_01_with_nengo(cur_sub,cur_tmpl,aggr, ...
    run_it,acc_rec,expl_mask,ow,path_to_nengo_output)
try
    %% PREPARATIONS
    % change into subs directory
    root_dir = pwd;
    cd(cur_sub)
    
    % name of this analysis determines the results directory
    results_ssdir = 'PDT_ss_design_01';
    
    % messaging
    disp(['Estimating subject ' cur_sub ' using PDT_ss_design_01'])
    
    % gain and loss ranges
    gain_min = 14;
    gain_max = 36;
    loss_max = 18;
    loss_min = 7;
    
    % for euclidean distance
    vec = [2;1;0];   % slope vector [gain; loss; 0]
    sp  = [26;13;0]; % sp support point; point on the diagonal [gain; loss; 0]
    
    % calculate the aggregated possible values
    cur_gains = (gain_min:2:gain_max);
    cur_losss = loss_min:loss_max    ;
    [vg,osg,nsg]  = agk_downsample_steps(cur_gains,aggr);
    gain_min  = min(vg);
    gain_max  = max(vg);
    [vl,osl,nsl]  = agk_downsample_steps(cur_losss,aggr);
    loss_min  = min(vl);
    loss_max  = max(vl);
    
    % read the nengo out put (model amy time series)
    nengo_amy = dlmread(path_to_nengo_output);
    
    %% FILLING BATCH
    % preparing the design information variables
    if acc_rec == 0
        names       = {'Pic.on','Pic.gam.on','Pic.gam.opt.on','feedback'};
        names_2     = names; % second session?
        onsets      = cell(1,4);
        onsets_2    = onsets;
        durations   = cell(1,4);
        durations_2 = durations;
    elseif acc_rec == 1 % needs to be adjusted later !
        names = {'accepted' 'rejected'};
        onsets = cell(1,2);
        durations = cell(1,2);
    end
    
    pmod = struct('name',{},'poly',{},'param',{});
    pmod(1).name{1} = 'gam';
    pmod(1).poly{1} = 1;
    pmod(1).name{2} = 'neg';
    pmod(1).poly{2} = 1;
    pmod(1).name{3} = 'pos';
    pmod(1).poly{3} = 1;
    
    % NENGO PARAM; here you need 100 ONSETS, FELIX
    % nga... nengo amygdala
    for kk = 2:51
        pmod(kk).name{1} = ['nga' num2str(kk)];
        pmod(kk).poly{1} = 1;
        pmod(kk).name{2} = 'gam';
        pmod(kk).poly{2} = 1;
        pmod(kk).name{3} = 'neg';
        pmod(kk).poly{3} = 1;
        pmod(kk).name{4} = 'pos';
        pmod(kk).poly{4} = 1;
    end

    
    for ii = 3:numel(names)
        pmod(ii).name{1} = 'gain';
        pmod(ii).poly{1} = 1;
        pmod(ii).name{2} = 'loss';
        pmod(ii).poly{2} = 1;
        pmod(ii).name{3} = 'ed';
        pmod(ii).poly{3} = 1;
        pmod(ii).name{4} = 'gam';
        pmod(ii).poly{4} = 1;
        pmod(ii).name{5} = 'neg';
        pmod(ii).poly{5} = 1;
        pmod(ii).name{6} = 'pos';
        pmod(ii).poly{6} = 1;
    end
    
    pmod_2=pmod;
    
    % get behav data
    behav_data = 0;
    try
        cd('Behav\PDT')
        load(ls('P*'))
        behav_data = 1;
    catch
        behav_data = 0;
        error_message = [results_ssdir ' ' cur_sub ' No P struct found.'];
        cd(root_dir)
        return
    end
    
    % prep the param mod regressors for each task_on regressor
    % for pic on
    p11 = []; % gam
    p12 = []; % neg
    p13 = []; % pos
    % for pic and gamble on
    p21 = []; % nengo amygd
    p22 = []; % gam
    p23 = []; % neg
    p24 = []; % pos
    % for pic and gamble and options on
    p31 = []; % gain
    p32 = []; % loss
    p33 = []; % ed
    p34 = []; % gam
    p35 = []; % neg
    p36 = []; % pos
    % for feedback (missings will be discarded anyways)
    p41 = []; % gain
    p42 = []; % loss
    p43 = []; % ed
    p44 = []; % gam
    p45 = []; % neg
    p46 = []; % pos
    
    % prep the param mod regressors for each task_on regressor (2nd run)
    % for pic on
    p11_2 = []; % gam
    p12_2 = []; % neg
    p13_2 = []; % pos
    % for pic and gamble on
    p21_2 = []; % nengo amygd
    p22_2 = []; % gam
    p23_2 = []; % neg
    p24_2 = []; % pos
    % for pic and gamble and options on
    p31_2 = []; % gain
    p32_2 = []; % loss
    p33_2 = []; % ed
    p34_2 = []; % gam
    p35_2 = []; % neg
    p36_2 = []; % pos
    % for feedback (missings will be discarded anyways)
    p41_2 = []; % gain
    p42_2 = []; % loss
    p43_2 = []; % ed
    p44_2 = []; % gam
    p45_2 = []; % neg
    p46_2 = []; % pos
    
    if acc_rec == 1 % needs to be reviewed still
        p21 = [];
        p22 = [];
        p23 = [];
    end
    
    % getting the param modulators split by acc_rec or not (NO NOT INCLUDED)
    
    % getting the gain mean and the loss mean
    for hh = 1:length(P.cur.choice)
        all_gains(hh) = str2double(cell2mat(P.gain.strings(P.cur.gamble{hh}(1))));
        all_losss(hh) = str2double(cell2mat(P.loss.strings(P.cur.gamble{hh}(2))));
    end
    mean_gain = mean(all_gains);
    mean_loss = mean(all_losss);
    tworuns   = isfield(P.t,'triggerscannertpostpause');
    
    if tworuns
        time_corr_second_run = (P.t.triggerscannert0 - P.t.triggerscannertpostpause); % new starting point in time for second run
    end
    try
        if behav_data == 1
            for ii = 1 : length(P.cur.choice)
                if tworuns
                    if ii <= 101 %% CHANGE!!!
                        if P.cur.choice(ii) > 0 && P.cur.choice(ii) < 5
                            % PIC ON
                            onsets(1)    = {[cell2mat(onsets(1)) P.t.cur_trial.stim_on(ii)]};
                            durations(1) = {[cell2mat(durations(1)) (P.t.cur_trial.stimLA_watch(ii) - P.t.cur_trial.stim_on(ii))]};
                            switch P.cur.cat(ii)
                                case 1
                                    cur_gam = 1;
                                    cur_neg = 0;
                                    cur_pos = 0;
                                case 2
                                    cur_gam = 0;
                                    cur_neg = 1;
                                    cur_pos = 0;
                                case 3
                                    cur_gam = 0;
                                    cur_neg = 0;
                                    cur_pos = 1;
                                case 6
                                    cur_gam = 0;
                                    cur_neg = 0;
                                    cur_pos = 0;
                            end
                            % pic dummy code
                            p11          = [p11 cur_gam];
                            p12          = [p12 cur_neg];
                            p13          = [p13 cur_pos];
                            
                            % PIC PLUS GAMBLE ON
                            % %% FELIX: this phase has to be a fixed number
                            % of onsets (100)
                            onsets(2)    = {[cell2mat(onsets(2)) P.t.cur_trial.stimLA_watch(ii)]};
                            durations(2) = {[cell2mat(durations(2)) (P.t.cur_trial.stimLA_on(ii) - P.t.cur_trial.stimLA_watch(ii))]};
                            % gain
                            cur_gain     = agk_recode(str2double(cell2mat(P.gain.strings(P.cur.gamble{ii}(1)))),osg,nsg) - mean_gain;
                            p21          = [p21 cur_gain];
                            % loss (here changed loss to abs. loss)
                            cur_loss     = agk_recode(abs(str2double(cell2mat(P.loss.strings(P.cur.gamble{ii}(2))))),osl,nsl) - abs(mean_loss);
                            p22          = [p22 cur_loss];
                            % ed
                            cur_point    = [cur_gain;cur_loss;0];
                            ed           = agk_get_ed(cur_point,sp,vec);
                            p23          = [p23 ed];
                            % pic dummy code
                            p24          = [p24 cur_gam];
                            p25          = [p25 cur_neg];
                            p26          = [p26 cur_pos];
                            
                            % PIC PLUS GAMBLE ON PLUS OPTIONS ON
                            onsets(3)    = {[cell2mat(onsets(3)) P.t.cur_trial.stimLA_on(ii)]};
                            durations(3) = {[cell2mat(durations(3)) (P.t.cur_trial.stimLA_off(ii) - P.t.cur_trial.stimLA_on(ii))]}; % or use: P.cur.rt(ii)
                            % gain, abs. loss, ed
                            p31          = [p31 cur_gain];
                            p32          = [p32 cur_loss];
                            p33          = [p33 ed];
                            % pic dummy code
                            p34          = [p34 cur_gam];
                            p35          = [p35 cur_neg];
                            p36          = [p36 cur_pos];
                            
                            % FEEDBACK
                            onsets(4)    = {[cell2mat(onsets(4)) P.t.cur_trial.stimLA_off(ii)]};
                            durations(4) = {[cell2mat(durations(4)) P.instr.slow.time]};
                            % gain, abs. loss, ed
                            p41          = [p41 cur_gain];
                            p42          = [p42 cur_loss];
                            p43          = [p43 ed];
                            % pic dummy code
                            p44          = [p44 cur_gam];
                            p45          = [p45 cur_neg];
                            p46          = [p46 cur_pos];
                        end
                    else % second run
                        if P.cur.choice(ii) > 0 && P.cur.choice(ii) < 5
                            % PIC ON
                            onsets_2(1)    = {[cell2mat(onsets_2(1)) P.t.cur_trial.stim_on(ii) + time_corr_second_run]};
                            durations_2(1) = {[cell2mat(durations_2(1)) (P.t.cur_trial.stimLA_watch(ii) - P.t.cur_trial.stim_on(ii))]};
                            switch P.cur.cat(ii)
                                case 1
                                    cur_gam = 1;
                                    cur_neg = 0;
                                    cur_pos = 0;
                                case 2
                                    cur_gam = 0;
                                    cur_neg = 1;
                                    cur_pos = 0;
                                case 3
                                    cur_gam = 0;
                                    cur_neg = 0;
                                    cur_pos = 1;
                                case 6
                                    cur_gam = 0;
                                    cur_neg = 0;
                                    cur_pos = 0;
                            end
                            % pic dummy code
                            p11_2          = [p11_2 cur_gam];
                            p12_2          = [p12_2 cur_neg];
                            p13_2          = [p13_2 cur_pos];
                            
                            % PIC PLUS GAMBLE ON
                            onsets_2(2)    = {[cell2mat(onsets_2(2)) P.t.cur_trial.stimLA_watch(ii) + time_corr_second_run]};
                            durations_2(2) = {[cell2mat(durations_2(2)) (P.t.cur_trial.stimLA_on(ii) - P.t.cur_trial.stimLA_watch(ii))]};
                            % gain
                            cur_gain       = agk_recode(str2double(cell2mat(P.gain.strings(P.cur.gamble{ii}(1)))),osg,nsg) - mean_gain;
                            p21_2          = [p21_2 cur_gain];
                            % loss (here changed loss to abs. loss)
                            cur_loss       = agk_recode(abs(str2double(cell2mat(P.loss.strings(P.cur.gamble{ii}(2))))),osl,nsl) - abs(mean_loss);
                            p22_2          = [p22_2 cur_loss];
                            % ed
                            cur_point      = [cur_gain;cur_loss;0];
                            ed             = agk_get_ed(cur_point,sp,vec);
                            p23_2          = [p23_2 ed];
                            % pic dummy code
                            p24_2          = [p24_2 cur_gam];
                            p25_2          = [p25_2 cur_neg];
                            p26_2          = [p26_2 cur_pos];
                            
                            % PIC PLUS GAMBLE ON PLUS OPTIONS ON
                            onsets_2(3)    = {[cell2mat(onsets_2(3)) P.t.cur_trial.stimLA_on(ii) + time_corr_second_run]};
                            durations_2(3) = {[cell2mat(durations_2(3)) (P.t.cur_trial.stimLA_off(ii) - P.t.cur_trial.stimLA_on(ii))]}; % or use: P.cur.rt(ii)
                            % gain, abs. loss, ed
                            p31_2          = [p31_2 cur_gain];
                            p32_2          = [p32_2 cur_loss];
                            p33_2          = [p33_2 ed];
                            % pic dummy code
                            p34_2          = [p34_2 cur_gam];
                            p35_2          = [p35_2 cur_neg];
                            p36_2          = [p36_2 cur_pos];
                            
                            % FEEDBACK
                            onsets_2(4)    = {[cell2mat(onsets_2(4)) P.t.cur_trial.stimLA_off(ii) + time_corr_second_run]};
                            durations_2(4) = {[cell2mat(durations_2(4)) P.instr.slow.time]}; % P.instr.slow.time is actually wrong; it's always a bit more than 500ms; use measrued times instead
                            % gain, abs. loss, ed
                            p41_2          = [p41_2 cur_gain];
                            p42_2          = [p42_2 cur_loss];
                            p43_2          = [p43_2 ed];
                            % pic dummy code
                            p44_2          = [p44_2 cur_gam];
                            p45_2          = [p45_2 cur_neg];
                            p46_2          = [p46_2 cur_pos];
                        end
                    end
                else % the normal case (not 2 PDT runs)
                    if P.cur.choice(ii) > 0 && P.cur.choice(ii) < 5
                        % PIC ON
                        onsets(1)    = {[cell2mat(onsets(1)) P.t.cur_trial.stim_on(ii)]};
                        durations(1) = {[cell2mat(durations(1)) (P.t.cur_trial.stimLA_watch(ii) - P.t.cur_trial.stim_on(ii))]};
                        switch P.cur.cat(ii)
                            case 1
                                cur_gam = 1;
                                cur_neg = 0;
                                cur_pos = 0;
                            case 2
                                cur_gam = 0;
                                cur_neg = 1;
                                cur_pos = 0;
                            case 3
                                cur_gam = 0;
                                cur_neg = 0;
                                cur_pos = 1;
                            case 6
                                cur_gam = 0;
                                cur_neg = 0;
                                cur_pos = 0;
                        end
                        % pic dummy code
                        p11          = [p11 cur_gam];
                        p12          = [p12 cur_neg];
                        p13          = [p13 cur_pos];
                        
                        % PIC PLUS GAMBLE ON
                        onsets(2)    = {[cell2mat(onsets(2)) P.t.cur_trial.stimLA_watch(ii)]};
                        durations(2) = {[cell2mat(durations(2)) (P.t.cur_trial.stimLA_on(ii) - P.t.cur_trial.stimLA_watch(ii))]};
                        % gain
                        cur_gain     = agk_recode(str2double(cell2mat(P.gain.strings(P.cur.gamble{ii}(1)))),osg,nsg) - mean_gain;
                        p21          = [p21 cur_gain];
                        % loss (here changed loss to abs. loss)
                        cur_loss     = agk_recode(abs(str2double(cell2mat(P.loss.strings(P.cur.gamble{ii}(2))))),osl,nsl) - abs(mean_loss);
                        p22          = [p22 cur_loss];
                        % ed
                        cur_point    = [cur_gain;cur_loss;0];
                        ed           = agk_get_ed(cur_point,sp,vec);
                        p23          = [p23 ed];
                        % pic dummy code
                        p24          = [p24 cur_gam];
                        p25          = [p25 cur_neg];
                        p26          = [p26 cur_pos];
                        
                        % PIC PLUS GAMBLE ON PLUS OPTIONS ON
                        onsets(3)    = {[cell2mat(onsets(3)) P.t.cur_trial.stimLA_on(ii)]};
                        durations(3) = {[cell2mat(durations(3)) (P.t.cur_trial.stimLA_off(ii) - P.t.cur_trial.stimLA_on(ii))]}; % or use: P.cur.rt(ii)
                        % gain, abs. loss, ed
                        p31          = [p31 cur_gain];
                        p32          = [p32 cur_loss];
                        p33          = [p33 ed];
                        % pic dummy code
                        p34          = [p34 cur_gam];
                        p35          = [p35 cur_neg];
                        p36          = [p36 cur_pos];
                        
                        % FEEDBACK
                        onsets(4)    = {[cell2mat(onsets(4)) P.t.cur_trial.stimLA_off(ii)]};
                        durations(4) = {[cell2mat(durations(4)) P.instr.slow.time]};
                        % gain, abs. loss, ed
                        p41          = [p41 cur_gain];
                        p42          = [p42 cur_loss];
                        p43          = [p43 ed];
                        % pic dummy code
                        p44          = [p44 cur_gam];
                        p45          = [p45 cur_neg];
                        p46          = [p46 cur_pos];
                    end
                end
            end
            
            pmod(1).param{1} = p11;
            pmod(1).param{2} = p12;
            pmod(1).param{3} = p13;
            
            pmod(2).param{1} = p21;
            pmod(2).param{2} = p22;
            pmod(2).param{3} = p23;
            pmod(2).param{4} = p24;
            pmod(2).param{5} = p25;
            pmod(2).param{6} = p26;
            
            pmod(3).param{1} = p31;
            pmod(3).param{2} = p32;
            pmod(3).param{3} = p33;
            pmod(3).param{4} = p34;
            pmod(3).param{5} = p35;
            pmod(3).param{6} = p36;
            
            pmod(4).param{1} = p41;
            pmod(4).param{2} = p42;
            pmod(4).param{3} = p43;
            pmod(4).param{4} = p44;
            pmod(4).param{5} = p45;
            pmod(4).param{6} = p46;
            
            % second run
            pmod_2(1).param{1} = p11_2;
            pmod_2(1).param{2} = p12_2;
            pmod_2(1).param{3} = p13_2;
            
            pmod_2(2).param{1} = p21_2;
            pmod_2(2).param{2} = p22_2;
            pmod_2(2).param{3} = p23_2;
            pmod_2(2).param{4} = p24_2;
            pmod_2(2).param{5} = p25_2;
            pmod_2(2).param{6} = p26_2;
            
            pmod_2(3).param{1} = p31_2;
            pmod_2(3).param{2} = p32_2;
            pmod_2(3).param{3} = p33_2;
            pmod_2(3).param{4} = p34_2;
            pmod_2(3).param{5} = p35_2;
            pmod_2(3).param{6} = p36_2;
            
            pmod_2(4).param{1} = p41_2;
            pmod_2(4).param{2} = p42_2;
            pmod_2(4).param{3} = p43_2;
            pmod_2(4).param{4} = p44_2;
            pmod_2(4).param{5} = p45_2;
            pmod_2(4).param{6} = p46_2;
        end
    catch MExc
        disp([results_ssdir ' ' cur_sub ' Problems with specifying param mod (behav data not ok?).']);
        error_message = {MExc,[results_ssdir ' ' cur_sub ' Problems with specifying param mod (behav data not ok?).']};
        cd(root_dir)
        return
    end
    
    % get the EPIs of this subject
    cd(root_dir)
    cd(cur_sub)
    try
        cd('MRT\NIFTI');
    catch
        disp([results_ssdir ' ' cur_sub ' no NIFTI dir. Skipping.'])
        error_message =  [results_ssdir ' ' cur_sub '  no NIFTI dir. Skipped.']; % ??? !!!
        return
    end
    found_epi = 0;
    cur_epi_dirs  = cellstr(ls('*_epi*'));
    cur_MoCo_dirs = cellstr(ls('*_MoCo*'));
    
    if ~tworuns
        if length(cur_epi_dirs{1}) > 0 % check if we're working with epis
            cd(cur_epi_dirs{1})        % PDT task is first EPI run
            found_epi = 1;
        elseif length(cur_MoCo_dirs{1}) > 0 % check if we're working with MoCo
            cd(cur_MoCo_dirs{1})            % PDT task is first EPI run
            found_epi = 1;
        end
        
        if found_epi == 1
            load(ls('Preprocessing*')); % for microtime resolution
            preproc_job = matlabbatch;
            cur_nslices   = preproc_job{1}.spm.temporal.st.nslices;
            cur_refslice  = preproc_job{1}.spm.temporal.st.refslice;
            % load template for ss model
            if acc_rec == 1
                load(cur_tmpl{2})
            elseif acc_rec == 0
                load(cur_tmpl{1})
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('FPList',pwd,'swuaepi'));
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = [];
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(spm_select('FPList',pwd,'^rp'));
        end
        
        if found_epi == 0 % incase no scans were found: throw away
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {''};
        end
        
    else % for two PDT runs
        
        % FIRST RUN
        found_epi = 0;
        if length(cur_epi_dirs{1}) > 0 % check if we're working with epis
            cd(cur_epi_dirs{1})     % PDT task is first EPI run
            found_epi = 1;
        elseif length(cur_MoCo_dirs{1}) > 0 % check if we're working with MoCo
            cd(cur_MoCo_dirs{1})    % PDT task is first EPI run
            found_epi = 1;
        end
        
        if found_epi == 1
            load(ls('Preprocessing*')); % for microtime resolution
            preproc_job = matlabbatch;
            cur_nslices   = preproc_job{1}.spm.temporal.st.nslices;
            cur_refslice  = preproc_job{1}.spm.temporal.st.refslice;
            % load template for ss model
            if acc_rec == 1
                % needs to be revised!
                load(cur_tmpl{2})
            elseif acc_rec == 0
                load(cur_tmpl{3})
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans  = cellstr(spm_select('FPList',pwd,'swuaepi'));
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = [];
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(spm_select('FPList',pwd,'^rp'));
        end
        
        if found_epi == 0 % incase no scans were found: throw away
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {''};
        end
        
        % SECOND RUN
        cd ..
        found_epi = 0;
        if length(cur_epi_dirs{1}) > 0 % check if we're working with epis
            cd(cur_epi_dirs{2})     % PDT task is first EPI run
            found_epi = 1;
        elseif length(cur_MoCo_dirs{1}) > 0 % check if we're working with MoCo
            cd(cur_MoCo_dirs{2})    % PDT task is first EPI run
            found_epi = 1;
        end
        
        if found_epi == 1
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans  = cellstr(spm_select('FPList',pwd,'swuaepi'));
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = [];
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = cellstr(spm_select('FPList',pwd,'^rp'));
        end
        
        if found_epi == 0 % in case no scans were found: throw away
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {''};
        end
    end
    
    
    cd .. % go back to NIFTI folder
    
    % create the results dir
    mkdir(pwd,'results')
    cd('results')
    
    if exist([pwd filesep results_ssdir])
        if ow == 1
            cmd_rmdir([pwd filesep results_ssdir])
        elseif ow == 0
            disp('Results dir already present. Overwrite not allowed, I will skip this subject.')
            error_message =  [results_ssdir ' ' cur_sub ' results dir already present. Skipped.'];
            return
        end
    end
    
    agk_mkdir_ex(pwd,results_ssdir)
    cd(results_ssdir)
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(pwd);
    
    if ~tworuns
        if behav_data == 1
            save('mult_cond.mat','names','onsets','durations','pmod');
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = cellstr(spm_select('FPList',pwd,'mult_cond.mat'));
        end
    else % in case of two sessions
        if behav_data == 1
            save('mult_cond.mat','names','onsets','durations','pmod');
            names     = names_2;
            onsets    = onsets_2;
            durations = durations_2;
            pmod      = pmod_2;
            save('mult_cond_2.mat','names','onsets','durations','pmod');
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = cellstr(spm_select('FPList',pwd,'mult_cond.mat'));
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = cellstr(spm_select('FPList',pwd,'mult_cond_2.mat'));
        end
    end
    
    % fill in microtime resolution
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = cur_nslices;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = cur_refslice;
    % specifiy an explicit mask, which voxels shall only be analysed
    matlabbatch{1}.spm.stats.fmri_spec.mask = cellstr(expl_mask);
    % specify TR-time
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.0; %% CHANGE!
    
    % contrast manager
    con.contrastType = {'t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t','t'};
    con.contrastNames = {'pic.on','pic.on.gam','pic.on.neg','pic.on.pos', ...
        'picgam.on','picgam.on.gain', 'picgam.on.loss','picgam.on.ed','picgam.on.gam','picgam.on.neg','picgam.on.pos', ...
        'picgamopt.on','picgamopt.on.gain', 'picgamopt.on.loss','picgamopt.on.ed','picgamopt.on.gam','picgamopt.on.neg','picgamopt.on.pos', ...
        'fb.on','fb.on.gain', 'fb.on.loss','fb.on.ed','fb.on.gam','fb.on.neg','fb.on.pos'};
    con.contrastWeights = agk_basic_t_cons(length(con.contrastType));
    if tworuns
        con.contrastWeights = agk_basic_t_cons_2sess(length(con.contrastType),length(con.contrastType));
    end
    con.contrastRep = {'none', 'none','none','none','none', ...
        'none', 'none','none','none','none', ...
        'none', 'none','none','none','none', ...
        'none', 'none','none','none','none', ...
        'none', 'none','none','none','none'
        };
    
    for jj = 1:length(con.contrastType)
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.name    = con.contrastNames{jj};
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.convec  = con.contrastWeights{jj};
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.sessrep = con.contrastRep{jj};
    end
    
    if behav_data == 1
        save('design.mat','matlabbatch');
    end
    
    try
        if run_it == 1
            spm_jobman('run',matlabbatch);
        end
        error_message = [results_ssdir ' ' cur_sub ' Estimation successfull.'];
    catch MExc
        error_message = {MExc,[results_ssdir ' ' cur_sub ' Estimation not successfull.']};
        cd(root_dir)
        return
    end
    cd(root_dir);
catch MExc
    disp('Something went wrong');
    error_message = MExc;
end
end