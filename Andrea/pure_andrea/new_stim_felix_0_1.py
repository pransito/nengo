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

fix_time = 200														 # a priori assumption for fixation times


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
		while cur_base_time < cur_gamble_length:		             # saccades only while gamble-watch period  <--- can be omitted if fixed fix_times are assumed
			for kk in range(len(cur_gamble)):                        # loop through gain and loss in this trial <--- still needs proper specification
				if cur_base_time < cur_gamble_length:
					cur_dic[cur_time] = cur_gamble[kk]
					cur_time = cur_time + fix_time 		             
					cur_base_time = cur_base_time + fix_time	
					cur_dic[cur_time] = 0                            # always set to zero after feature fix
					cur_time = cur_time + fix_time
					cur_base_time = cur_base_time + fix_time 		 # can also be omitted if the while statement is dropped due to fixed fix_times
		# measure how much time left to end of trial
		time_left = cur_gamble_length - cur_base_time
		cur_time  = cur_time + time_left # + len_act_rew
	return cur_dic

cur_dic = new_stim_time(gamble_list,gamble_length,iti_length,gamble_start_times,fix_time)