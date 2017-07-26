AUTOMAKE_OPTIONS = foreign

#optional_nlopt
if DISABLE_NLOPT
NLOPT_CPPFLAGS = -DDISABLE_NLOPT
else
NLOPT_CPPFLAGS =
endif

#beta optimization stuff
if BETAOPTTOGETHER
BETAOPTFLAGS = -DBETAOPTTOGETHER
endif

if BETAOPTSEPARATE
BETAOPTFLAGS = -DBETAOPTSEPARATE
endif

if BETAOPTBROKEN
BETAOPTFLAGS = -DBETAOPTBROKEN
endif

AM_CPPFLAGS = $(GSL_CFLAGS) $(BETAOPTFLAGS) $(NLOPT_CPPFLAGS) $(NLOPT_CFLAGS)
bin_PROGRAMS = seq2expr seqannot seq2expr_streaming
seq2expr_SOURCES = seq2expr.cpp FactorIntFunc.cpp FactorIntFunc.h ExprModel.cpp ExprPar.cpp ExprPar.h ExprModel.h PredictorTrainer.h PredictorTrainer.cpp ExprPredictor.cpp ExprPredictor.h SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h IO.cpp IO.h ObjFunc.cpp ObjFunc.h ExprFunc.cpp ExprFunc.h DataSet.cpp DataSet.h
seq2expr_LDADD = $(GSL_LIBS) $(NLOPT_LIBS)
seqannot_SOURCES = seqannot.cpp SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h
seqannot_LDADD = $(GSL_LIBS) $(NLOPT_LIBS)
seq2expr_streaming_SOURCES = seq2expr_streaming.cpp FactorIntFunc.cpp FactorIntFunc.h ExprModel.cpp ExprPar.cpp ExprPar.h ExprModel.h PredictorTrainer.h PredictorTrainer.cpp ExprPredictor.cpp ExprPredictor.h SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h IO.cpp IO.h ObjFunc.cpp ObjFunc.h ExprFunc.cpp ExprFunc.h DataSet.cpp DataSet.h protocol/gemstat.pb.cc protocol/gemstat.pb.h protocol/motif.pb.cc protocol/motif.pb.h protocol/math_types.pb.cc protocol/math_types.pb.h protocol/sequence.pb.cc protocol/sequence.pb.h
seq2expr_streaming_LDADD = $(GSL_LIBS) $(NLOPT_LIBS) -L/usr/local/lib -lprotobuf
