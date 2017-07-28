import sys
sys.path.append("src/protocol")
sys.path.append("../GEMSTAT_scripts/python/src/")

import gemstat_pb2

import scipy

#TARGET EXPRESSION
def make_labeled_matrix(row_names,in_data,colnames = None):
	N_i = len(row_names)
	N_j = in_data.shape[1]
	assert N_i == in_data.shape[0]
	target_condnames=map(str,range(1,N_j+1))

	ret_labeled_data = gemstat_pb2.GSLabeledMatrixMessage()

	for n in row_names:
		ret_labeled_data.row_names.append(n)
	for cn in target_condnames:
		ret_labeled_data.col_names.append(cn)

	stovo = ret_labeled_data.storage
	stovo.matrix_storage_type = stovo.DENSE

	stovo.i = N_i
	stovo.j = N_j

	for val in in_data.ravel():
		stovo.dense_storage.append(val)
	return ret_labeled_data


def make_sequence(name,seq):
	seq_message = gemstat_pb2.sequence__pb2.SequenceMessage()
	seq_message.name = name
	seq_message.seq = seq
	return seq_message

def make_random_sequence(name,length,alph=["A","C","G","T"]):
	data = scipy.random.randint(0,len(alph),length)
	seq = "".join([alph[i] for i in data])
	return make_sequence(name,seq)

#MOTIF
def make_a_motif(in_name,in_pseudo,in_counts):
	a = gemstat_pb2.motif__pb2.MotifMessage()
	a.name = in_name
	a.pseudocount = in_pseudo

	for i in range(in_counts.shape[0]):
		pos = a.counts.add()
		pos.A = in_counts[i,0]
		pos.C = in_counts[i,1]
		pos.G = in_counts[i,2]
		pos.T = in_counts[i,3]

	return a

def make_sparse_bool_mat(in_bool_mat,ret_mat=None):
	if ret_mat is None:
		ret_mat = gemstat_pb2.math__types__pb2.SparseBoolMatrix()

	for i in range(in_bool_mat.shape[0]):
		for j in range(in_bool_mat.shape[1]):
			#an_entry = gemstat_pb2.math__types__pb2.SparseBoolMatrix.SparseBoolStorage()

			if in_bool_mat[i,j]:
				an_entry = ret_mat.storage.add()
				an_entry.i = i
				an_entry.j = j
				an_entry.val = in_bool_mat[i,j]
			#ret_mat.storage.add(an_entry)
	return ret_mat

from itertools import izip_longest, count
def make_param_msg(in_pardict):
	par = gemstat_pb2.ExprParMessage()
	for name, data in in_pardict["tfs"]:
		#tf_entry = gemstat_pb2.ExprParMessage.TFDataMessage()
		tf_entry = par.tf_data.add()
		tf_entry.k = data[0]
		tf_entry.alpha_txp = data[1]
		tf_entry.alpha_rep = data[2]
		#par.tf_data.append(tf_entry)

	for qbtm, pi, beta in izip_longest(in_pardict["qbtms"],in_pardict["pis"],in_pardict["betas"]):
		if None is qbtm:
			qbtm = in_pardict["qbtms"][0]
		if None is pi:
			pi = in_pardict["pis"][0]
		if None is beta:
			beta = in_pardict["betas"][0]

		#promoter_message = gemstat_pb2.ExprParMessage.PromoterDataMessage()
		promoter_message = par.promoter_data.add()
		promoter_message.basal_transcription = qbtm
		promoter_message.pi = pi
		promoter_message.beta = beta

		#par.promoter_data.append(promoter_message)

	for i in in_pardict["thresholds"]:
		par.tf_annotation_cutoffs.append(i)

	tf_index_dict = dict(zip([n for n,_ in in_pardict["tfs"]],count()))
	par.interaction_weights.i=len(in_pardict["tfs"])
	par.interaction_weights.j=len(in_pardict["tfs"])
	par.interaction_weights.matrix_storage_type = par.interaction_weights.SPARSE
	for (name_a,name_b),interact_val in in_pardict["coops"]:
		sp_entry = par.interaction_weights.sparse_storage.add()
		sp_entry.i = tf_index_dict[name_a]
		sp_entry.j = tf_index_dict[name_b]
		sp_entry.val = interact_val

	return par

def make_a_model(model_option_str,act,rep,coop,repress,coop_dist_thr,rep_dist_thr,coop_type="BINARY",shared_scaling=False,oq=False,max_contact=1):
	m = gemstat_pb2.ExprModelMessage()
	m.model_option = gemstat_pb2.ExprModelMessage.ModelOptionMessage.Value(model_option_str)
	m.one_qbtm_per_crm = oq
	m.shared_scaling = shared_scaling
	make_sparse_bool_mat(coop, m.coop_mat)
	make_sparse_bool_mat(repress, m.repression_mat)
	m.max_contact = max_contact
	m.repression_dist_thr = rep_dist_thr

	m.interaction_function.interaction_type = gemstat_pb2.FactorIntFuncMessage.FactorIntType.Value(coop_type)
	m.interaction_function.coop_dist_thr = coop_dist_thr

	for i in act:
		m.act_indicators.append(i)
	for i in rep:
		m.rep_indicators.append(i)

	return m



#### MAIN #####
import sys
import os
import Bio.SeqIO as BSIO

sys.path.append

COOP_DIST = 50
REP_DIST = 250


N_PREDICTIONS = 1

#DIRNAME="."
DIRNAME="data_streamtest"
foo = scipy.loadtxt(os.path.join(DIRNAME,"factor_info.txt"),dtype=scipy.str_)
factor_names = foo[:,0]
factor_act_indicators = scipy.array(foo[:,1],dtype=scipy.int_)
factor_rep_indicators = scipy.array(foo[:,2],dtype=scipy.int_)

enhancer_seqrecords = [i for i in BSIO.parse(os.path.join(DIRNAME,"seqs.fa"),"fasta")]

for i in enhancer_seqrecords:
	#print "'{}'".format(i.name)
	#print "'{}'".format(i.seq.tostring())
	pass

#sys.exit(1)

factor_pseudo = [0.25,1.0,1.0,1.0,1.0,1.0]
factors = list()
for one_factor_n,one_factor_p in zip(factor_names,factor_pseudo):
	factors.append(
		(one_factor_n,one_factor_p,
				scipy.loadtxt(
					os.path.join(DIRNAME,"factors","{}.txt".format(one_factor_n))
				)
		)
	)

def load_tabfile(in_filename):
	foo = scipy.loadtxt(in_filename,dtype=scipy.str_)
	foo_data = scipy.array(foo[1:,1:],dtype=scipy.float_)
	foo_rownames = [i.strip() for i in foo[1:,0]]
	foo_colnames = [i.strip() for i in foo[0,1:]]
	return foo_rownames, foo_colnames, foo_data


labeled_exp_data_foo = load_tabfile(os.path.join(DIRNAME,"expr.tab"))
labeled_input_data_foo = load_tabfile(os.path.join(DIRNAME,"factor_expr.tab"))

coop_pairs = scipy.loadtxt(os.path.join(DIRNAME,"coop.txt"),dtype=scipy.str_)


import gemstat.model as GMOD
example_parfile = GMOD.ParFile.read_1_6a(os.path.join(DIRNAME,"example.par"))
#print example_parfile

test_parameters = make_param_msg(example_parfile)
#test_str = test_parameters.SerializeToString()
#print test_str
#exit(1)

#N_CONDITIONS = 30
#ENHANCER_NAMES=["en1","en2"]
#MOTIF_NAMES=["zld","baz"]
#EN_LENGTH=500



#test_model = make_a_model("DIRECT",[True,False],[False,True],
#	scipy.array([[True,False],[False,True]],dtype=scipy.bool_),scipy.zeros((2,2),dtype=scipy.bool_),
#	100,250)


some_enhancers = [make_sequence(i.name,i.seq.tostring().strip()) for i in enhancer_seqrecords]
some_motifs = [make_a_motif(*i) for i in factors]
labeled_exp_data = make_labeled_matrix(labeled_exp_data_foo[0],labeled_exp_data_foo[2])
labeled_input_data = make_labeled_matrix(labeled_input_data_foo[0],labeled_input_data_foo[2])

coop_matrix = scipy.zeros((len(some_motifs),len(some_motifs)),dtype=scipy.bool_)
from itertools import count
factor_name_to_index = dict(zip(factor_names,count()))
for a,b in coop_pairs:
	coop_matrix[factor_name_to_index[a],factor_name_to_index[b]] = True
	coop_matrix[factor_name_to_index[b],factor_name_to_index[a]] = True


rep_matrix = scipy.zeros((len(some_motifs),len(some_motifs)),dtype=scipy.bool_)

test_model = make_a_model("DIRECT",factor_act_indicators,factor_rep_indicators,coop_matrix,rep_matrix,
				COOP_DIST,REP_DIST)



def send_message(stream,name,msg=None):
	stream.write("{}\n".format(name))
	if msg is not None:
		serialized_str = msg.SerializeToString()
		stream.write("{}\n".format(1+len(serialized_str)))
		stream.write("{}\n".format(serialized_str))

#with open("./test_serialized.out","w") as outfile:
with sys.stdout as outfile:

	send_message(outfile,"TARGET_EXPRESSION",labeled_exp_data)

	send_message(outfile,"FACTOR_CONCENTRATIONS",labeled_input_data)

	send_message(outfile,"CLEAR_SEQUENCES")
	for one_enhancer in some_enhancers:
		send_message(outfile,"SEQUENCE",one_enhancer)

	send_message(outfile,"CLEAR_MOTIFS")
	for one_motif in some_motifs:
		send_message(outfile,"MOTIF",one_motif)


	send_message(outfile,"MODEL",test_model)

	send_message(outfile,"PARAMETER",test_parameters)

	send_message(outfile,"ANNOTATE")

	send_message(outfile,"CHECKPOINT")
	for i in range(N_PREDICTIONS):
		send_message(outfile,"PREDICT")
