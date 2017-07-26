import sys
sys.path.append("src/protocol")

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


N_CONDITIONS = 30
ENHANCER_NAMES=["en1","en2"]
MOTIF_NAMES=["zld","baz"]
EN_LENGTH=500

some_enhancers = [make_random_sequence(one_name,EN_LENGTH) for one_name in ENHANCER_NAMES]
some_motifs = [make_a_motif(n,1.0,scipy.absolute(10.0*scipy.randn(8,4))) for n in MOTIF_NAMES]

labeled_exp_data = make_labeled_matrix(ENHANCER_NAMES,scipy.absolute(scipy.randn(len(ENHANCER_NAMES),N_CONDITIONS)))
labeled_input_data = make_labeled_matrix(MOTIF_NAMES,scipy.absolute(scipy.randn(len(MOTIF_NAMES),N_CONDITIONS)))

test_model = make_a_model("DIRECT",[True,False],[False,True],
	scipy.array([[True,False],[False,True]],dtype=scipy.bool_),scipy.zeros((2,2),dtype=scipy.bool_),
	100,250)

def send_message(stream,name,msg=None):
	stream.write("{}\n".format(name))
	if msg is not None:
		serialized_str = msg.SerializeToString()
		stream.write("{}\n".format(1+len(serialized_str)))
		stream.write("{}\n".format(serialized_str))

with open("./test_serialized.out","w") as outfile:

	send_message(outfile,"TARGET_EXPRESSION",labeled_exp_data)

	send_message(outfile,"FACTOR_CONCENTRATIONS",labeled_input_data)

	send_message(outfile,"CLEAR_SEQUENCES")
	for one_enhancer in some_enhancers:
		send_message(outfile,"SEQUENCE",one_enhancer)

	send_message(outfile,"CLEAR_MOTIFS")
	for one_motif in some_motifs:
		send_message(outfile,"MOTIF",one_motif)
	send_message(outfile,"ANNOTATE")

	send_message(outfile,"MODEL",test_model)

	send_message(outfile,"CHECKPOINT")
	send_message(outfile,"PREDICT")
