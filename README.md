# muzziolab_sleep_tetrodes_objecttask
 Code to perform the analysis of the tetrode data from the object task for the sleep paper in 2021. Requires Mulana package.

# Requirement of the data
There must be a text file called "experiment_description" for each mouse dataset. Here is an example:
"""{
	"animal":"ANIMAL001",
	"imaging_region":"CA1",
	"experiment":"object_task_consecutive_trials",
	"apparatus_type": "neuralynx_tetrodes",
	"has_digs": 0,
	"arena":{
		"shape":"square",
		"length_cm":35.0
	},
	"num_contexts":1,
	"session_folders":["s1","s2","s3","s4","s5"],
	"mclust_tfile_bits":32,
	"nvt_file_trial_separation_threshold_s": 10.0,
	"nvt_filename": "VT1.nvt"
}"""

