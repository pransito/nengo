## FUNCTIONS FOR ANDREA BG MODEL ##
import 	nengo
from	nengo.dists import Uniform
import 	numpy as np
from 	nengo.utils.functions import piecewise
import  matplotlib.pyplot as plt
import  random as rn
import  collections
import  andrea_func as abf
from    nengo.processes import WhiteSignal

## BG Model functions

# normalize input to BG
def normalize(x):
    OldRange = x.max() - x.min()
    if x.size == 1:
        return 0.8
    if OldRange == 0:
        x_new = np.ones(x.size)*0.8
    else:
        NewRange = (0.8 - 0.5)
        x_new = (((x - x.min())*0.3)/OldRange)+0.5
    return x_new

class EnvironmentReward():

        '''
        Here we define how the environment behaves
        '''
    
	def __init__(self):
		self.reward_gamble0_win = 1
		self.reward_gamble0_lose = -1
		self.reward_gamble0_win_prob = 0.3
		self.reward_gamble0_lose_prob = 1 - self.reward_gamble0_win_prob
		self.reward_gamble1_win = 0.2
		self.reward_gamble1_lose = -0.2
		self.reward_gamble1_win_prob = 0.7
		self.reward_gamble1_lose_prob = 1 - self.reward_gamble1_win_prob

		self.current_reward  =  [0, 0]
		self.choice_time     = 0
		self.reward_duration = 0.01
		self.len_act_rew     = 0.5
		self.len_stim        = 0.8
		self.len_pause       = 0.1
		self.more_than_one   = 0
		self.wait_time       = self.len_act_rew + self.len_pause + self.len_stim
		self.cur_choice      = [0,0]

	def gamble(self, t, x):
		self.cur_choice = x[0:2]
		
		# change prob of gambling
		#if t%15 == 0:
		#	self.reward_gamble0_win_prob = 0.7
		#	self.reward_gamble0_lose_prob = 1 - self.reward_gamble0_win_prob
		#	self.reward_gamble1_win_prob = 0.3
		#	self.reward_gamble1_lose_prob = 1 - self.reward_gamble1_win_prob
			
		if self.more_than_one is 0 and t > (self.choice_time + self.wait_time - self.len_act_rew):
			#print 'First Trial!'
			## determine reward
			max_index  = np.argmax(self.cur_choice)
			if self.cur_choice[max_index] > 0.5:
				#print 'a choice was made'
				self.more_than_one = 1
				self.choice_time = t
				if max_index == 1:   # gamble 2
					if np.random.random() < self.reward_gamble1_win_prob:
						self.current_reward[1] = x[2]
					else:
						self.current_reward[1] = x[3]
				elif max_index == 0: # gamble 1					
					if np.random.random() < self.reward_gamble0_win_prob:
						self.current_reward[0] = self.reward_gamble0_win
					else:
						self.current_reward[0] = self.reward_gamble0_lose
				else:
					self.current_reward =  [0, 0]
			else:
				#print 'no choice was made'
				if self.cur_choice[0] < 0.2 and self.cur_choice[1] < 0.2:
					self.current_reward = [0, 0]
			return self.current_reward
			##
		elif (self.more_than_one is 1) and t > (self.choice_time + self.wait_time):
			#print '2nd or higher trial!'
			## determine reward
			max_index  = np.argmax(self.cur_choice)
			if self.cur_choice[max_index] > 0.5:
				#print 'a choice was made'
				self.choice_time = t
				if max_index == 1:   # gamble 2
					if np.random.random() < self.reward_gamble1_win_prob:
						self.current_reward[1] = x[2]
					else:
						self.current_reward[1] = x[3]
				elif max_index == 0: # gamble 1					
					if np.random.random() < self.reward_gamble0_win_prob:
						self.current_reward[0] = self.reward_gamble0_win
					else:
						self.current_reward[0] = self.reward_gamble0_lose
				else:
					self.current_reward =  [0, 0]
			else:
				#print 'no choice was made'
				if self.cur_choice[0] < 0.2 and self.cur_choice[1] < 0.2:
					self.current_reward = [0, 0]
			return self.current_reward
			##
		elif (self.more_than_one is 1) and t > (self.choice_time + self.len_act_rew):
			self.current_reward = [0,0]
			return self.current_reward
		else:
			return self.current_reward
		
	def gamble_old(self, t, x):
	
		#x[0] = gamble0
		#x[1] = gamble1
		print t
		print x
		print self.current_reward
		if cmp(self.current_reward, [0,0]) is 0:
			print 'current reward is [0,0]!'
			cur_choice = x[0:2]
			max_index  = np.argmax(cur_choice)
			if cur_choice[max_index] > 0.5:
				print 'a choice was made'
				if max_index == 1:   # gamble 2
					if np.random.random() < self.reward_gamble1_win_prob:
						self.current_reward[1] = x[2]
					else:
						self.current_reward[1] = x[3]
				elif max_index == 0: # gamble 1
					
					if np.random.random() < self.reward_gamble0_win_prob:
						self.current_reward[0] = self.reward_gamble0_win
					else:
						self.current_reward[0] = self.reward_gamble0_lose
				else:
					self.current_reward =  [0, 0]
			else:
				print 'no choice was made'
				if self.cur_choice[0] < 0.2 and self.cur_choice[1] < 0.2:
					self.current_reward = [0, 0]
		if cmp(self.current_reward,[0,0]) is not 0:
			self.current_reward = [0,0]
			return self.current_reward
		else:
			return [0,0]
# ANDREA model functions
def get_maxrates(n,x,y): # do not understand this; why everywhere (?)s
    return ((np.random.rand(n)*x+y).astype(int))

def input_stim(cur_list,n_trials):
	"""Takes in a list of n same-length lists.
	Every of those lists has as elements values.
	E.g. put in a vector of possible gains and 
	possible losses.
	Samples from those vectors with replacement,
	to produce tupels which are of same dim as cur_list."""
	output_list = []
	
	for ii in range(n_trials):
		cur_trial = []
		for jj in range(len(cur_list)):
			cur_vec       = cur_list[jj,:]
			cur_trial.append(cur_vec[rn.randrange(1,len(cur_vec), 1)])
		output_list.append(cur_trial)
	return output_list

def stim_time(cur_list,stim_length,pause_length,len_act_rew, fix_time):
	"""Takes in a list of stimuli vectors (sim to tuples).
	Wants to know how long the stimuli tuples are shown.
	And length of pause.
	And length of action__reward_period.
	And how long each fixation time is.
	After every fixations will go to zero,
	to avoid PEs dependent of difference between,
	features. Give times so that 0 period is possible,
	towards end of trial."""
    # initialize the dictionary
	cur_dic = {0:0}
		
	# for every trial create the stimulation
	n_trials  = len(cur_list)
	cur_time  = 0
	for ii in range(n_trials):
		cur_stim          = cur_list[ii]
		cur_stim_length   = stim_length             # add jitter here later
		cur_pause_length  = pause_length            # add jitter here later
		cur_dic[cur_time] = 0
		cur_time          = cur_time + pause_length # first add the pause
		cur_base_time     = 0
		while cur_base_time < cur_stim_length:		# saccades only as long as overall stim time not over
			for kk in range(len(cur_stim)):
				if cur_base_time < cur_stim_length:
					cur_dic[cur_time] = cur_stim[kk]
					cur_time          = cur_time + fix_time 		# overall time
					cur_base_time     = cur_base_time + fix_time	# within-trial time
					cur_dic[cur_time] = 0                           # always set to zero after feature fix
					cur_time          = cur_time + fix_time
					cur_base_time     = cur_base_time + fix_time
		# measure how much time left to end of trial
		time_left = cur_stim_length - cur_base_time
		cur_time  = cur_time + time_left + len_act_rew
	return cur_dic

# Felix stim function 
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
		while cur_base_time < 3.2:		         	 # saccades only while gamble-watch period  <--- can be omitted if fixed fix_times are assumed
			for kk in range(len(cur_gamble)):                        # loop through gain and loss in this trial <--- still needs proper specification
				if cur_base_time < 3.5:
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
	
def control_time(n_trials,stim_length,pause_length,len_act_rew):
	"""Makes a dictionary to set control time for integrator.
	Pause time and lenght of action and reward: 1: resetting (to 0; i.e. integrator OFF)
	Value time: 								0: integrating. (integrator ON)
	Starts with pause."""
    # initialize the dictionary
	cur_dic  = {0:1}
	cur_time = 0
		
	# for every trial create the stimulation
	for ii in range(n_trials):
		cur_stim_length   = stim_length                 # add jitter here later
		cur_pause_length  = pause_length                # add jitter here later
		cur_dic[cur_time+0.05] = 1
		cur_time          = cur_time + cur_pause_length # first add the pause
		cur_dic[cur_time+0.05] = 0
		cur_time          = cur_time + cur_stim_length # first add the pause
		cur_time          = cur_time + len_act_rew
	return cur_dic
	
# FELIX CONTROL_TIME FUNCTION
def new_control_time(iti_on, stim_on, LA_watch_on, LA_on, LA_off, contr_iti = 1, 
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
    # initialize the dictionary - starting with NO integration
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
	
def control_PE(cur_list,stim_length,pause_length, len_act_rew, fix_time):
	"""Makes a dictionary to set control time for ignoring/
	taking into account the OFC signal.
	ignore-time: 0
	eval time:   1
	To avoid taking into account PEs when just setting to baseline."""
    # initialize the dictionary
	cur_dic = {0:1}
		
	# for every trial create the stimulation
	n_trials  = len(cur_list)
	cur_time  = 0
	for ii in range(n_trials):
		cur_stim          = cur_list[ii]
		cur_stim_length   = stim_length             # add jitter here later
		cur_pause_length  = pause_length            # add jitter here later
		cur_dic[cur_time] = 0
		cur_time          = cur_time + pause_length # first add the pause
		cur_base_time     = 0
		while cur_base_time < cur_stim_length:		# saccades only as long as overall stim time not over
			for kk in range(len(cur_stim)):
				if cur_base_time < cur_stim_length:
					cur_dic[cur_time] = 0
					cur_time          = cur_time + fix_time 		# overall time
					cur_base_time     = cur_base_time + fix_time	# within-trial time
					cur_dic[cur_time] = 1                           # always set to zero after feature fix
					cur_time          = cur_time + fix_time
					cur_base_time     = cur_base_time + fix_time
		# measure how much time left to end of trial
		time_left = cur_stim_length - cur_base_time
		cur_time  = cur_time + time_left + len_act_rew
	return cur_dic
	
def new_control_PE(gamble_length, gamble_start_times, fix_time, ignore_time,early_start):
	"""Makes a dictionary to set control time for ignoring/
	taking into account the OFC signal.
	ignore-time: 1
	eval time:   0
	To avoid taking into account PEs when just setting to baseline.
	Needs:
	- Vector with the lengths of the presentations of stimuli tuples
	- The respective start times of these presentations
	- The length of the single fixation
	- The length of the time window while the prediction error 
	  calculation shall be suspended
	- early_start: how many seconds before jump to baseline should ignore PE time start?
	"""
    # initialize the dictionary - starting with evaluation
	ctrl_PE_dic = {0:0}
		
	# for every trial create the stimulation
	n_trials  = len(gamble_length)
	cur_time  = 0
	for ii in range(n_trials):
		cur_gamble_length = gamble_length[ii]
		ctrl_PE_dic[cur_time] = 0									 
		# cur_time = cur_time + cur_iti_length 						 # Is omitted, as in this version, the next line let's cur_time jump to stimulus onset time directly.
		cur_time = gamble_start_times[ii]							 # Jump to stimulus onset time directly
		cur_base_time = 0											 # Within-Trial time <- is reset on everey iteration/trial.
		while cur_base_time < cur_gamble_length:		             # saccades only while gamble-watch period  
			cur_time = cur_time + fix_time 
			cur_base_time = cur_base_time + fix_time					
			ctrl_PE_dic[cur_time-early_start] = 1	   					         # start ignoring at end of first fixation 
			cur_time = cur_time + ignore_time						 
			cur_base_time = cur_base_time + ignore_time
			ctrl_PE_dic[cur_time] = 0                                # end ignoring after the duration of ignore_time
			cur_time = cur_time - ignore_time + fix_time 
			cur_base_time = cur_base_time - ignore_time + fix_time   # wait for the fixation on the middle to end before start new
	    ## measure how much time left to end of trial
		#time_left = cur_gamble_length - cur_base_time				 # the next two lines could be omitted, but might be useful if a break should be included above ( in the while loop) 
		#cur_time  = cur_time + time_left
	return ctrl_PE_dic

# FUNCTIONS FOR NEURONS #
def simple_product(x):
    return x[0]*x[1]
def simple_product_ctrl(x):
    return x[0]*x[1]
def simple_product_tau(x):
	return x[0]*x[1]*tau
def triple_product(x):
    return x[0]*x[1]*x[2]
def neg_triple_product(x):
    return (x[0]*x[1]*x[2])*(-1)
def neg_product(x):
    return -x[0]*x[1]
def neg_product_ctrl(x):
    return -x[0]*x[1]