import scipy.io as sio

# load P structure and extract the relevant values:
# 'squeezing' eliminates singular (redundant) levels in the structure
# omitting the transformation to 'record', preserves the MATLAB-like structure/accessing of elements

VP_mat = sio.loadmat('P_0901', squeeze_me = True, struct_as_record = False)
VP_P = VP_mat['P']

gamble_list = VP_P.cur.gamble                                        # STILL NEED TO BE TRANSLATED PROPERLY (cheap workaround below)
gamble_length = VP_P.gamble_watch.durations							 # refers to the duration of the presentation of the gamble (without decision scale?)
iti_length = VP_P.cur_iti_dur_corr 
gamble_start_times = VP_P.t.cur_trial.stimLA_watch

fix_time = 0.2														 # a priori assumption for fixation times
ignore_time = 0.01

# TIMEPOINTS IN TRIAL:
# (0) ONSET ITI
# (1) OFFSET ITI / STIM ONSET
# (2) OFFSET STIM / LA_WATCH ONSET
# ---- TIME WINDOW: GAMBLE_WATCH.DURATIONS
# (3) ONSET LA PRESS (BU PRESS POSSIBLE)
# (4) OFFSET LA (== BUTTON PRESS) / ONSETT FEEDBACK
# (5) OFFSET FEEDBACK / ONSET (NEXT) ITI

# Let's work through this: 
# (0):
iti_on = VP_P.t.cur_trial.iti_on
# (1):
stim_on = VP_P.t.cur_trial.stim_on
# (2):
LA_watch_on = VP_P.t.cur_trial.stimLA_watch # ( == gamble_start_times)
# (3):
LA_on = VP_P.t.cur_trial.stimLA_on
## Alternative: LA_on = LA_watch_on + VP_P.gamble_watch.durations	# This is not 100% equvalent as there are smallish deviations (M = 19.8ms, SD = 4.3ms)) probably due to PTB delays
# (4):
LA_off = VP_P.t.cur_trial.stimLA_off
# (5):																			# This is reversely defined anyways, so not really needed
## iti_on(+1)																	# Problems, however, might arise with the last trials before breaks/end; CHECK FOR THIS!!


## the following is really ugly, but works to transfer the gamble values
for p in range(0,2):
	for q in range(len(gamble_list)):
		gamble_list[q] = gamble_list[q].astype(float)
		if p == 0:
			gamble_list[q][p] = ((gamble_list[q][p]*2) + 12) / 20
		if p == 1:
			gamble_list[q][p] = ((gamble_list[q][p]+6) * -1) / 20
		


def new_stim_time(gamble_list, gamble_length, iti_length, gamble_start_times, fix_time):     # len_act_rew
	"""Takes in a list of 2-element-arrays (1.: Gain, 2.: Loss).
	Wants to know how long the stimuli tuples are shown.
	And length of ITI.
	And the time of stimulus-onset of the stimuli-tuples.
	And how long each fixation time is.
	Value is always zero, unless during assumed fixations on
	either Gain or Loss stimulus.
	Gives out dictionnary combining times with values of 
	stimulus representation."""
	
    # initialize the dictionary
	cur_dic = {0:0}
		
	# for every trial create the stimulation
	n_trials = len(gamble_list)
	cur_time = 0
	for ii in range(n_trials):
		cur_gamble = gamble_list[ii]
		cur_gamble_length = gamble_length[ii]		                 
		# cur_iti_length = iti_length[ii]             				 # At the moment, this doesn't do anything useful. 
		cur_dic[cur_time] = 0
		# cur_time = cur_time + cur_iti_length 						 # Is omitted, as in this version, the next line let's cur_time jump to stimulus onset time directly.
		cur_time = gamble_start_times[ii]
		cur_base_time = 0											 # Within-Trial time <- is reset on everey iteration/trial.
		while cur_base_time < cur_gamble_length:		             # saccades only while gamble-watch period  
			for kk in range(len(cur_gamble)):                        # loop through gain and loss in this trial <--- still needs proper specification
				if cur_base_time < cur_gamble_length:
					cur_dic[cur_time] = cur_gamble[kk]
					cur_time = cur_time + fix_time 		             
					cur_base_time = cur_base_time + fix_time	
					cur_dic[cur_time] = 0                            # always set to zero after feature fix
					cur_time = cur_time + fix_time
					cur_base_time = cur_base_time + fix_time 		 # can also be omitted if the while statement is dropped due to fixed fix_times
		# measure how much time left to end of trial
		time_left = cur_gamble_length - cur_base_time				 # these two lines can probably be omitted?! 
		cur_time  = cur_time + time_left # + len_act_rew 			  
	return cur_dic

cur_dic = new_stim_time(gamble_list,gamble_length,iti_length,gamble_start_times,fix_time)





def control_time(iti_on, stim_on, LA_watch_on, LA_on, LA_off, contr_iti = 1, 
				contr_stim = 1, contr_LA_watch = 0, contr_LA = 0, 
				contr_feedback = 1, post_integ = 0.05, max_post_integ = 1):
	"""Makes a dictionary to set control time for integrator.
	1: resetting (to 0; i.e. integrator OFF)
	0: integrating. (integrator ON)
	Starts NON-INTEGRATING.
	Needs:
	(a) Five vectors with the starting times of the 5 trial subsections
	-> (*_on vectors)
	(b) Logical value for each trial subsection, declaring whether 
	during this section the integrator shall be on (0) or off (1) 
	-> contr_* vectors
		default: 	subsection 			integration
					ITI					OFF
					stimulus			OFF
					LA watch			ON
					LA 					ON
					feedback phase		OFF
	(c) post_integ: a scalar that defines the delay with which 
	integrating starts and ends (in seconds; default: 0.05)
	(d) max_post_integ: the maximum of post integration in seconds
	(e.g. after the last trial; default: 1)
	"""
    # initialize the dictionary - starting with NO inegration
	contr_dic  = {0:1}
	n_trials = len(stim_on)														# Taking the length of one of the timing vectors
		
	# let's run over all trials
	for ii in range(n_trials):
		contr_dic[iti_on[ii]+post_integ] = contr_iti
		contr_dic[stim_on[ii]+post_integ] = contr_stim
		contr_dic[LA_watch_on[ii]+post_integ] = contr_LA_watch
		contr_dic[LA_on[ii]+post_integ] = contr_LA
		contr_dic[LA_off[ii]+post_integ] = contr_feedback
		
		# In the following, I want to control for trials which were 
		# followed by a irregularly long interval 
		# (e.g. last trial, before breaks) 
		
		if ii < n_trials-1:														# you can easily add the trials before the breaks here, if needed; now this only controls for the very last trial.
			if LA_off[ii] + max_post_integ < iti_on[ii+1]:
				contr_dic[LA_off[ii] + max_post_integ] = 1						# if the 'break' is longer than max_post_integ, the integrater is stopped automatically
		else:
			contr_dic[LA_off[ii] + post_integ] = 1								# after the last trial, the integrater is stopped immediately 
		
	return contr_dic
	
	
def new_control_PE(gamble_length, gamble_start_times, fix_time, ignore_time):
	"""Makes a dictionary to set control time for ignoring/
	taking into account the OFC signal.
	ignore-time: 0
	eval time:   1
	To avoid taking into account PEs when just setting to baseline.
	Needs:
	- Vector with the lengths of the presentations of stimuli tuples
	- The respective start times of these presentations
	- The length of the single fixation
	- The length of the time window while the prediction error 
	  calculation shall be suspended
	"""
    # initialize the dictionary - starting with evaluation
	ctrl_PE_dic = {0:1}
		
	# for every trial create the stimulation
	n_trials  = len(gamble_length)
	cur_time  = 0
	for ii in range(n_trials):
		cur_gamble_length = gamble_length[ii]
		ctrl_PE_dic[cur_time] = 1									 
		# cur_time = cur_time + cur_iti_length 						 # Is omitted, as in this version, the next line let's cur_time jump to stimulus onset time directly.
		cur_time = gamble_start_times[ii]							 # Jump to stimulus onset time directly
		cur_base_time = 0											 # Within-Trial time <- is reset on everey iteration/trial.
		while cur_base_time < cur_gamble_length:		             # saccades only while gamble-watch period  
			cur_time = cur_time + fix_time 
			cur_base_time = cur_base_time + fix_time					
			ctrl_PE_dic[cur_time] = 0	   					         # start ignoring at end of first fixation 
			cur_time = cur_time + ignore_time						 
			cur_base_time = cur_base_time + ignore_time
			ctrl_PE_dic[cur_time] = 1                                # end ignoring after the duration of ignore_time
			cur_time = cur_time - ignore_time + fix_time 
			cur_base_time = cur_base_time - ignore_time + fix_time   # wait for the fixation on the middle to end before start new
	    ## measure how much time left to end of trial
		#time_left = cur_gamble_length - cur_base_time				 # the next two lines could be omitted, but might be useful if a break should be included above ( in the while loop) 
		#cur_time  = cur_time + time_left
	return ctrl_PE_dic
	
	
	