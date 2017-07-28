/*****************************************************
 * Train and test the expression model
 * Input: sequence, expression, motif, factor expr, cooperativity rules,
 *   activator list, repression rules
 * Output: trained model and expression of training data
 * File formats:
 * (1) Sequence: Fasta format
 * (2) Expression: one line per sequence
 *       <seq_name expr_1 ... expr_m>
 * (3) Motif: Stubb format
 * (4) Factor expression: one line per factor
 *       <factor expr_1 ... expr_m>
 * (5) Cooperativity: list of cooperative factor pairs
 * (6) Factor roles: the role field is either 1 (yes) or 0 (no)
 *       <factor activator_role repressor_role>
 * (7) Repression: list of <R A> where R represses A
 * (8) Parameters: the format is:
 *       <factor binding activation repression>
 *       <factor1 factor2 coop>
 *       <basal_transcription = x>
 *     where repression is optional, and the coop. lines are optional.
 * Note that (5), (6), (7) and (8) may be empty
 ******************************************************/
#include "Utils.h"
#include "IO.h"

#include "ExprModel.h"
#include "ExprPredictor.h"

#include "ObjFunc.h"

#include "protocol/gemstat.pb.h"
#include "protocol/sequence.pb.h"

#include <iostream>

#include "SeqAnnotator.h"

static string read_to_string(istream &stream, uint32_t count)
{
    std::string result(count, ' ');
    stream.read(&result[0], count);

    return result;
}

static string message_str_from_stream(istream &stream){
    std::string limits;
    std::getline(cin,limits);
    cout << "READ LIMITS " << limits << endl;
    int the_limit = std::stoi(limits);

    std::string result(the_limit, ' ');
    stream.read(&result[0], the_limit);
    return result;
}

Matrix message_to_labeled_matrix(const gemstat::GSLabeledMatrixMessage& in_labeled_matrix_msg, vector<string>& rownames, vector<string>& colnames);
Motif message_to_motif(const gemstat::MotifMessage& in_motif_msg, vector< double > in_background);
ExprModel* message_to_new_exprmodel(const gemstat::ExprModelMessage& in_model_msg, vector< Motif >& default_motifs);

void delete_model(ExprModel*& model_to_delete){
  delete &(model_to_delete->coopMat);
  delete &(model_to_delete->repressionMat);
  delete &(model_to_delete->actIndicators);
  delete &(model_to_delete->repressionMat);
  delete model_to_delete;
  model_to_delete = NULL;
}



void free_fix_from_par(const ExprPar& param_ff, vector<bool>& indicator_bool);

int main( int argc, char* argv[] )
{
    // command line processing
    string seqFile, r_seqFile, annFile, exprFile, motifFile, factorExprFile, coopFile, factorInfoFile, repressionFile, parFile, axis_wtFile;
    string outFile;                               // output file
    string dnase_file;
    string factor_thr_file;
    string par_out_file; // the learned parameters will get stored here
    ofstream par_out_stream; // Uninitialized at first.

    ModelType cmdline_modelOption = LOGISTIC;
    double coopDistThr = 50;
    double factorIntSigma = 50.0;                 // sigma parameter for the Gaussian interaction function
    double repressionDistThr = 250;
    int maxContact = 1;
    bool read_factor_thresh = false;
    double eTF = 0.60;
    unsigned long initialSeed = time(0);

    double l1 = 0.0;
    double l2 = 0.0;

    bool cmdline_one_qbtm_per_crm = false;
    bool cmdline_one_beta = false;

    bool cmdline_write_gt = true;

    string lower_bound_file; ExprPar lower_bound_par; bool lower_bound_par_read = false;
    string upper_bound_file; ExprPar upper_bound_par; bool upper_bound_par_read = false;
    string free_fix_indicator_filename;
    ExprPar::one_qbtm_per_crm = false;
    ExprFunc::one_qbtm_per_crm = false;

    ExprPredictor::nAlternations = 5;

    GOOGLE_PROTOBUF_VERIFY_VERSION;

    for ( int i = 1; i < argc; i++ )
    {

        if ( !strcmp( "-oo", argv[ i ] ) )
            ExprPredictor::objOption = getObjOption( argv[++i] );
        else if ( !strcmp( "-wt", argv[ i ]) )
            axis_wtFile = argv[ ++ i ];
        else if ( !strcmp( "-na", argv[i] ) )
            ExprPredictor::nAlternations = atoi( argv[++i] );
        else if ( !strcmp( "-ff", argv[i] ) )
            free_fix_indicator_filename = argv[++i];
        else if ( !strcmp( "-df", argv[ i ]))
            dnase_file = argv[ ++i ];
        else if ( !strcmp( "-ft", argv[ i ]))
            factor_thr_file = argv[ ++i ];
	else if ( !strcmp( "--seed", argv[ i ]))
	    initialSeed = atol( argv[++i] );
  else if ( !strcmp("-l1", argv[ i ]))
      l1 = atof(argv[ ++i ]);
  else if ( !strcmp("-l2", argv[ i ]))
      l2 = atof(argv[ ++i ]);
	else if ( !strcmp("-lower_bound", argv[ i ]))
	    lower_bound_file = argv[ ++i ];
  else if ( !strcmp("-upper_bound", argv[ i ]))
	    upper_bound_file = argv[ ++i ];
        else if( !strcmp("-no_gt_out", argv[ i ]))
            cmdline_write_gt = false;
    }

    //     bool readSites = false;     // whether read sites (if true) or read sequences

    // additional control parameters
    double gcContent = 0.5;
    FactorIntType intOption = BINARY;             // type of interaction function
    ExprPar::searchOption = CONSTRAINED;          // search option: unconstrained; constrained.

    ExprPredictor::nRandStarts = 0;
    ExprPredictor::min_delta_f_SSE = 1.0E-10;
    ExprPredictor::min_delta_f_Corr = 1.0E-10;
    ExprPredictor::min_delta_f_CrossCorr = 1.0E-10;
    ExprPredictor::nSimplexIters = 400;
    ExprPredictor::nGradientIters = 50;

    int rval;
    vector< vector< double > > data;              // buffer for reading matrix data
    vector< string > labels;                      // buffer for reading the labels of matrix data
    string factor1, factor2;                      // buffer for reading factor pairs

    // read the sequences
    vector< Sequence > seqs;
    vector< string > seqNames;
    vector< SiteVec > seqSites( seqs.size() );
    vector< int > seqLengths( seqs.size() );

    //read the random sequences
    vector< Sequence > r_seqs;
    vector< string > r_seqNames;
    int r_nSeqs = r_seqs.size();




    // read the motifs
    vector< Motif > motifs;
    vector< string > motifNames;
    vector< double > background = createNtDistr( gcContent );

    // read the factor expression data
    Matrix factorExprData;

    // read the expression data
    vector< string > condNames;
    Matrix exprData;
    int nConds = exprData.nCols();

    vector < string > expr_condNames = condNames;

    //TODO
    //Initialize the dataset that is actually provided
    DataSet *training_dataset = NULL;

    vector < int > axis_start;
    vector < int > axis_end;
    vector < double > axis_wts;


    //TODO: R_SEQ Either remove this feature or un-comment it.
    //site representation of the random sequences
    vector< SiteVec > r_seqSites( r_nSeqs );
    vector< int > r_seqLengths( r_nSeqs );

    ExprModel *the_model = NULL;
    ExprPredictor *predictor = NULL;
    ParFactory *param_factory = NULL;
    ExprPar the_parameters;
    vector <bool> indicator_bool(0,false);


    bool dataset_dirty = true;
    bool param_factory_dirty = true;
    bool predictor_dirty = true;
    bool annotations_dirty = false;
    bool user_annotations = false;
    /************************************
    * Begin parsing the interactive input
    *************************************/
    std::string command_line;
    vector< string > tokens;
    while(std::getline(cin,command_line)){
      tokens.clear();
      std::istringstream line_ss(command_line);
      copy(istream_iterator<string>(line_ss),
            istream_iterator<string>(),
            back_inserter(tokens));

    cout << "COMMAND " << tokens << endl;

    if(0 == tokens[0].compare("EXIT")){
      cout << "Goodbye." << endl;
      break;
    }
    else if(0 == tokens[0].compare("MODEL")){
        std::string one_message = message_str_from_stream(cin);
        gemstat::ExprModelMessage input_model_msg = gemstat::ExprModelMessage();
        input_model_msg.ParseFromString(one_message);

        if( NULL != the_model ) {
          delete_model(the_model);
          the_model = NULL;
        }

        the_model = message_to_new_exprmodel(input_model_msg, motifs);
        //In case motifs were updated
        motifNames.clear();
        for(int k=0;k<motifs.size();k++){
          motifNames.push_back(motifs[k].getName());
        }


        param_factory_dirty = true;
    }
    else if(0 == tokens[0].compare("TARGET_EXPRESSION")){
    //TODO read_expression
      std::string one_message = message_str_from_stream(cin);
      gemstat::GSLabeledMatrixMessage read_a_labeled_matrix = gemstat::GSLabeledMatrixMessage();
      read_a_labeled_matrix.ParseFromString(one_message);

      vector<string> target_rownames;
      vector<string> target_colnames;

      if(NULL != training_dataset){
        delete training_dataset;
        training_dataset = NULL;
      }

      exprData = message_to_labeled_matrix(read_a_labeled_matrix,target_rownames,target_colnames);
      dataset_dirty = true;
    }
    else if(0 == tokens[0].compare("FACTOR_CONCENTRATIONS")){
      //read factor concentrations
      std::string one_message = message_str_from_stream(cin);
      gemstat::GSLabeledMatrixMessage read_a_labeled_matrix = gemstat::GSLabeledMatrixMessage();
      read_a_labeled_matrix.ParseFromString(one_message);

      vector<string> factor_rownames;
      vector<string> factor_colnames;

      if(NULL != training_dataset){
        delete training_dataset;
        training_dataset = NULL;
      }

      factorExprData = message_to_labeled_matrix(read_a_labeled_matrix,factor_rownames,factor_colnames);
      dataset_dirty = true;
    }
    else if(0 == tokens[0].compare("CLEAR_MOTIFS")){
      motifs.clear();
      motifNames.clear();
      for(int i;i<seqSites.size();i++){//also invalidates annotations
        seqSites[i].clear();
      }
      annotations_dirty=true;
    }
    else if(0 == tokens[0].compare("MOTIF")){
      //read one motif
      std::string one_message = message_str_from_stream(cin);
      gemstat::MotifMessage cmdl_motif_message = gemstat::MotifMessage();
      cmdl_motif_message.ParseFromString(one_message);

      Motif mot_from_msg = message_to_motif(cmdl_motif_message,background);

      motifs.push_back( mot_from_msg );
      motifNames.push_back( mot_from_msg.getName() );
      annotations_dirty = true;
    }
    else if(0 == tokens[0].compare("CLEAR_SEQUENCES")){
      //vector< Sequence > seqs;
      //vector< string > seqNames;
      //int nSeqs = seqs.size();
      seqs.clear();
      seqNames.clear();
      seqLengths.clear();
      seqSites.clear();
      //nSeqs = 0;
      annotations_dirty = true;
    }
    else if(0 == tokens[0].compare("SEQUENCE")){
      //read_and_annotate_sequences()
      std::string one_message = message_str_from_stream(cin);
      gemstat::SequenceMessage cmdl_seq_message = gemstat::SequenceMessage();
      cmdl_seq_message.ParseFromString(one_message);
      string tmp_seqname = string(cmdl_seq_message.name());

      Sequence tmp_sequence = Sequence(cmdl_seq_message.seq());

      seqs.push_back(tmp_sequence);
      seqNames.push_back(tmp_seqname);
      seqLengths.push_back(tmp_sequence.size());
      seqSites.push_back(SiteVec());
      //nSeqs += 1;
      annotations_dirty=true;
    }
    else if(0 == tokens[0].compare("ANNOTATIONS")){
      //Take straight annotations
      cout << "Sorry, this feature is not yet implemented." << endl; continue;
      user_annotations = true;
      annotations_dirty = false;
    }else if(0 == tokens[0].compare("ANNOTATE")){
      //force annotation
      annotations_dirty = true;
    }else if(0 == tokens[0].compare("WEIGHTS")){
      //input_ax_weights();
      cerr << "WEIGHTS command is not implemented" << endl;
    }
    else if(0 == tokens[0].compare("PARAMETER")){
      std::string one_message = message_str_from_stream(cin);
      gemstat::ExprParMessage par_msg = gemstat::ExprParMessage();
      par_msg.ParseFromString(one_message);

      //TODO: check vs old to see if we need to invalidate the annotations
      //the_parameters
      ExprPar new_parameters = param_factory->create_expr_par(); //Currently, code further down expects par_init to be in PROB_SPACE.
      new_parameters = param_factory->changeSpace(new_parameters, PROB_SPACE);

      cerr << "TF DATA SIZE " << par_msg.tf_data_size() << endl;
      //WTF WTF WTF? tf_data_size is 1 too big, but it seems right on the python end.
      //damnit
      for(int i = 0;i<par_msg.tf_data_size();i++){

        gemstat::ExprParMessage::TFDataMessage one_tf_entry = par_msg.tf_data(i);
        if(i>=new_parameters.maxBindingWts.size()){
          cerr << "DEBUG Warning received extra tf parameters when getting par vector." << endl;
          cerr << "i=" << i << " : " << one_tf_entry.k() << " : " << one_tf_entry.alpha_txp() << " : " << one_tf_entry.alpha_rep() << endl;
          continue;
        }

        //cerr << "GyA : " << one_tf_entry.k() << " : " << one_tf_entry.alpha_txp() << " : " << one_tf_entry.alpha_rep() << endl;
        new_parameters.maxBindingWts[i] = one_tf_entry.k();
        new_parameters.txpEffects[i] = one_tf_entry.alpha_txp();
        new_parameters.repEffects[i] = one_tf_entry.alpha_rep();
      }

      cerr << "PROMOTER DATA SIZE " << par_msg.promoter_data_size() << endl;
      for(int i = 0;i<par_msg.promoter_data_size();i++){
        gemstat::ExprParMessage::PromoterDataMessage one_prom_data = par_msg.promoter_data(i);
        if(i < new_parameters.basalTxps.size()){ new_parameters.basalTxps[i] = one_prom_data.basal_transcription(); }
        if(i < new_parameters.pis.size()){ new_parameters.pis[i] = one_prom_data.pi(); }
        if(i < new_parameters.betas.size()){ new_parameters.betas[i] = one_prom_data.beta(); }
      }

      //cerr << "DEBUG " << par_msg.has_interaction_weights() << " : " << par_msg.interaction_weights().sparse_storage_size() << endl;
      //cerr << "DEBUG " << the_parameters.factorIntMat; cerr << endl;
      //cerr << "DEBUG new_interact" << new_parameters.factorIntMat.nRows() << endl << new_parameters.factorIntMat; cerr << endl;
      gemstat::GSMatrixMessage sparse_interaction_matrix = par_msg.interaction_weights();
      for(int i = 0;i<sparse_interaction_matrix.sparse_storage_size();i++){
        gemstat::GSMatrixMessage::SparseEntry one_sparse_entry = sparse_interaction_matrix.sparse_storage(i);
        new_parameters.factorIntMat.setElement(one_sparse_entry.i(),one_sparse_entry.j(),one_sparse_entry.val());
        new_parameters.factorIntMat.setElement(one_sparse_entry.j(),one_sparse_entry.i(),one_sparse_entry.val());
      }


      for(int i = 0;i<par_msg.tf_annotation_cutoffs_size();i++){
        new_parameters.energyThrFactors[i] = par_msg.tf_annotation_cutoffs(i);
      }
      //Check if the thresholds differ, which would mean that the annotations need to be updated.
      for(int i = 0;i<new_parameters.energyThrFactors.size();i++){
        annotations_dirty = annotations_dirty || (the_parameters.energyThrFactors[i] != new_parameters.energyThrFactors[i]);
      }

      cerr << the_parameters.energyThrFactors << endl;
      cerr << new_parameters.energyThrFactors << endl;

      the_parameters = ExprPar(new_parameters);

    }
    else if(0 == tokens[0].compare("CHECKPOINT")){

        // CHECK POINT
        cout << "CHECKPOINT" << endl;
        cout << "Sequences:" << endl;
        for ( int i = 0; i < seqs.size(); i++ ) cout << seqNames[i] << endl << seqs[i] << endl;
        cout << "Site representation of sequences:" << endl;
        for ( int i = 0; i < seqs.size(); i++ ) {
          cout << ">" << seqNames[i] << endl;
          for ( int j = 0; j < seqSites[i].size(); j++ ){cout << seqSites[i][j] << endl;}
        }


        cout << "Expression: " << endl << exprData << endl;
        cout << "Factor expression:" << endl << factorExprData << endl;

        if(NULL == the_model){
          cout << "Factor motifs:" << endl;
          for ( int i = 0; i < motifs.size(); i++ ) cout << motifNames[i] << endl << motifs[i] << endl;
        }

        if(NULL != the_model){
          cout << "MODEL" << endl;
          cout << "Model Type: " << endl << getModelOptionStr(the_model->modelOption) << endl;
          cout << "Cooperativity matrix:" << endl << the_model->coopMat << endl;
          cout << "Activators:" << endl << the_model->actIndicators << endl;
          cout << "Repressors:" << endl << the_model->repIndicators << endl;
          cout << "Repression matrix:" << endl << the_model->repressionMat << endl;

          cout << "Coop distance threshold = " << the_model->intFunc->getMaxDist() << endl;

          if ( the_model->modelOption == QUENCHING || the_model->modelOption == CHRMOD_LIMITED )
          {
              cout << "Maximum_Contact = " << the_model->maxContact << endl;
          }
          if ( the_model->modelOption == QUENCHING || the_model->modelOption == CHRMOD_LIMITED || the_model->modelOption == CHRMOD_UNLIMITED )
          {
              cout << "Repression_Distance_Threshold = " << the_model->repressionDistThr << endl;
          }
          cout << "Factor motifs (in model):" << endl;
          for ( int i = 0; i < motifs.size(); i++ ) cout << motifNames[i] << endl << motifs[i] << endl;

        }

        cout << "PARAMETERS " << endl;
        cout << the_parameters.maxBindingWts << endl;
        cout << the_parameters.txpEffects << endl;
        cout << the_parameters.repEffects << endl;


      }else if(0 == tokens[0].compare("PREDICT")){
        // create the expression predictor

        if(annotations_dirty){
          cerr << "ERROR : somehow, annotations are still dirty!" << endl;
          continue;
        }
        if(NULL == predictor){
          cout << "ERROR no predictor!" << endl;
          continue;
        }

      //DO PREDICTION
      //TODO:
      //int predict( const SiteVec& targetSites, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const;
      vector< vector<double > > prediction_output;
      predictor->predict_all(the_parameters, prediction_output);



      //TODO: Actually output the prediction
      for(int i = 0;i<seqs.size();i++){
          cout << seqNames[i] << endl;
          for(int j=0;j<prediction_output[i].size();j++){ prediction_output[i][j] *= the_parameters.betas[i]; }
          cout << prediction_output[i] << endl;
      }

    }//END OF if tree , now handle things that were dirtied

    if(dataset_dirty){
      cerr << "CREATING DATASET " << endl;
      predictor_dirty = true;
      if(NULL != training_dataset){
        delete training_dataset;
        training_dataset = NULL;
      }

      if(factorExprData.nCols() == exprData.nCols()){
        training_dataset = new DataSet(factorExprData,exprData);
        dataset_dirty = false;
      }
    }

    if(param_factory_dirty && NULL != the_model){
      if( NULL != param_factory){
        delete param_factory;
        param_factory = NULL;
      }
      cerr << "UPDATING PARAM FACTORY" << endl;
      param_factory = new ParFactory(*the_model, seqs.size());
      the_parameters = param_factory->create_expr_par(); //Currently, code further down expects par_init to be in PROB_SPACE.
      the_parameters = param_factory->changeSpace(the_parameters, PROB_SPACE);

      param_factory_dirty = false;
      predictor_dirty = true;
    }

    //Handle DIRTY state.
    //ANNOTATION
    if((!user_annotations) && annotations_dirty && NULL != the_model ){
      //Annotate the sequences we have/MO
      //TODO: need to get this from the current model.

      if(motifs.size() != the_parameters.energyThrFactors.size()){
        cerr << "ERROR trying to annotate, motifs and pearameters don't agree." << endl;
        cerr << the_parameters.energyThrFactors << endl;
      }else{
        cerr << "DOING ANNOTATION" << endl;

        vector < double > tmp_energyThrFactors(6, 0.6);
        //SeqAnnotator ann( motifs, the_parameters.energyThrFactors );
        SeqAnnotator ann(motifs, tmp_energyThrFactors);
        for ( int i = 0; i < seqs.size(); i++ )
        {
            seqSites[i].clear();
            ann.annot( seqs[ i ], seqSites[ i ] );
        }
        annotations_dirty = false;
        predictor_dirty = true;
      }
    }

    if(predictor_dirty){
      if(NULL != predictor){
        delete predictor;
        predictor = NULL;
      }

      if( predictor_dirty && NULL != the_model && !dataset_dirty && !param_factory_dirty && (user_annotations || !annotations_dirty)){
        predictor = new ExprPredictor( seqs, seqSites, r_seqSites, seqLengths, r_seqLengths, *training_dataset, motifs, *the_model, indicator_bool, motifNames, axis_start, axis_end, axis_wts );
        predictor_dirty = false;
      }
    }
    //Do these later
    /*
    load_init_params();

    load_free_fix();


    load_upper_lower_bounds();

    setup_regularization();
    setup_bounds();

    train_the_predictor();
    */

    }

    google::protobuf::ShutdownProtobufLibrary();
    return 0;
}
//END MAIN

Matrix message_to_labeled_matrix(const gemstat::GSLabeledMatrixMessage& in_labeled_matrix_msg, vector<string>& rownames, vector<string>& colnames){
    gemstat::GSMatrixMessage its_storage = in_labeled_matrix_msg.storage();

    int max_i = its_storage.i();
    int max_j = its_storage.j();

    Matrix will_return(its_storage.i(),its_storage.j());

    rownames.clear();
    colnames.clear();

    for(int k=0;k<in_labeled_matrix_msg.row_names_size();k++){
      rownames.push_back(in_labeled_matrix_msg.row_names(k));
    }

    for(int k=0;k<in_labeled_matrix_msg.col_names_size();k++){
      colnames.push_back(in_labeled_matrix_msg.col_names(k));
    }

    if(its_storage.matrix_storage_type() == gemstat::GSMatrixMessage::DENSE){ //load from dense storage
      int storage_total_size = its_storage.dense_storage_size();

      for(int k=0;k<storage_total_size;k++){
        will_return.setElement(k/max_j,k%max_j,its_storage.dense_storage(k));
      }
    }else if(its_storage.matrix_storage_type() == gemstat::GSMatrixMessage::SPARSE){//load from sparse
      int sparse_total_size = its_storage.sparse_storage_size();
      for(int k=0;k<sparse_total_size;k++){
        gemstat::GSMatrixMessage::SparseEntry one_sparse_entry = its_storage.sparse_storage(k);
        will_return.setElement(one_sparse_entry.i(),one_sparse_entry.j(),one_sparse_entry.val());
      }
    }else{
      assert(false);
    }
    return will_return;
}

Motif message_to_motif(const gemstat::MotifMessage& in_motif_msg,vector< double > in_background){

  Matrix countMat( in_motif_msg.counts_size(), 4 );
  for(int i = 0;i< in_motif_msg.counts_size();i++){
    const gemstat::MotifMessage::MotifPosition& one_pos = in_motif_msg.counts(i);
    countMat(i,0) = one_pos.a();
    countMat(i,1) = one_pos.c();
    countMat(i,2) = one_pos.g();
    countMat(i,3) = one_pos.t();
  }

  Motif tmp_motif(countMat,in_motif_msg.pseudocount(),in_background);
  tmp_motif.setName(in_motif_msg.name());
  return tmp_motif;
}


ExprModel* message_to_new_exprmodel(const gemstat::ExprModelMessage& in_model_msg, vector< Motif >& default_motifs){

    int nActInd = in_model_msg.act_indicators_size();
    int nRepInd = in_model_msg.rep_indicators_size();
    assert(nActInd == nRepInd);
    int nFactors = nActInd;


    IntMatrix *coopMat = new IntMatrix( nFactors, nFactors, false );

    // read the roles of factors
    vector< bool > *actIndicators = new vector<bool>( nFactors, true );
    vector< bool > *repIndicators = new vector<bool>( nFactors, false );

    //get the factorinfo

    // read the repression matrix
    IntMatrix *repressionMat = new IntMatrix( nFactors, nFactors, false );
    FactorIntFunc* intFunc;

    gemstat::SparseBoolMatrix msg_coop_mat = in_model_msg.coop_mat();
    for(int k=0;k<msg_coop_mat.storage_size();k++){
      gemstat::SparseBoolMatrix::SparseBoolStorage one_storage = msg_coop_mat.storage(k);
      coopMat->setElement(one_storage.i(),one_storage.j(),one_storage.val());
    }

    for(int k=0;k<in_model_msg.act_indicators_size();k++){
      (*actIndicators)[k] = in_model_msg.act_indicators(k);
    }

    for(int k=0;k<in_model_msg.rep_indicators_size();k++){
      (*repIndicators)[k] = in_model_msg.rep_indicators(k);
    }

    gemstat::SparseBoolMatrix msg_rep_mat = in_model_msg.repression_mat();
    for(int k=0;k<msg_rep_mat.storage_size();k++){
      gemstat::SparseBoolMatrix::SparseBoolStorage one_storage = msg_rep_mat.storage(k);
      repressionMat->setElement(one_storage.i(),one_storage.j(),one_storage.val());
    }

    gemstat::FactorIntFuncMessage the_int_option = in_model_msg.interaction_function();

    switch(the_int_option.interaction_type()){
      case BINARY:
        intFunc = new FactorIntFuncBinary( the_int_option.coop_dist_thr() );
        break;
      case GAUSSIAN:
        assert(the_int_option.has_factor_int_sigma());
        intFunc = new FactorIntFuncGaussian( the_int_option.coop_dist_thr(), the_int_option.factor_int_sigma() );
        break;
      case HELICAL:
        intFunc = new FactorIntFuncHelical( the_int_option.coop_dist_thr() );
        break;
      default:
        assert(false);
        break;
    }

    ExprModel *ret_ptr = NULL;
    //Create a new ExprModel with all of the selected options.
    gemstat::ExprModelMessage_ModelOptionMessage model_option_msg = in_model_msg.model_option();
    string model_asked_for_in_msg_str = gemstat::ExprModelMessage_ModelOptionMessage_Name(in_model_msg.model_option());
    //DEBUG
    cout << "MODEL TYPE " << model_asked_for_in_msg_str << endl;
    ModelType what_kind_of_model = getModelOption( model_asked_for_in_msg_str );

    //TODO: Continue here
    ret_ptr = new ExprModel( what_kind_of_model, in_model_msg.one_qbtm_per_crm(),
      default_motifs, intFunc, in_model_msg.max_contact(),
      *coopMat, *actIndicators, *repIndicators, *repressionMat,
      in_model_msg.repression_dist_thr());
    ret_ptr->shared_scaling = in_model_msg.shared_scaling();
    return ret_ptr;
}


/**

  //Create a new ExprModel with all of the selected options.
  //TODO: Continue here

}




void setup_regularization(){
  //Setup regularization
  if(0.0 != l1 || 0.0 != l2){
    cerr << "INFO: Regularization was turned on and will be used. l1 = " << l1 << " l2 = " << l2 << " ."<< endl;

    ExprPar tmp_centers = predictor->param_factory->create_expr_par();
    ExprPar tmp_l1 = predictor->param_factory->create_expr_par();
    ExprPar tmp_l2 = predictor->param_factory->create_expr_par();

    //TODO: add an option to read l1 and l2 values from a file.
    vector< double > tmp_l12_vector;
    tmp_l1.getRawPars(tmp_l12_vector);
    std::fill(tmp_l12_vector.begin(),tmp_l12_vector.end(),l1);
    tmp_l1 = predictor->param_factory->create_expr_par(tmp_l12_vector, ENERGY_SPACE);

    tmp_l2.getRawPars(tmp_l12_vector);
    std::fill(tmp_l12_vector.begin(),tmp_l12_vector.end(),l2);
    tmp_l2 = predictor->param_factory->create_expr_par(tmp_l12_vector, ENERGY_SPACE);


    RegularizedObjFunc *tmp_reg_obj_func = new RegularizedObjFunc(predictor->trainingObjective,
                                            tmp_centers,
                                            tmp_l1,
                                            tmp_l2
                                          );
    predictor->trainingObjective = tmp_reg_obj_func;
  }
}

void load_upper_lower_bounds(){
  if ( !upper_bound_file.empty() ){
try{
  upper_bound_par = param_factory->load( upper_bound_file );
  upper_bound_par = param_factory->changeSpace(upper_bound_par, ENERGY_SPACE);
  upper_bound_par_read = true;
}catch (int& e){
  cerr << "Cannot read upper bounds from " << upper_bound_file << endl;
  exit( 1 );
}
  }

  if ( !lower_bound_file.empty() ){
try{
  lower_bound_par = param_factory->load( lower_bound_file );
  lower_bound_par = param_factory->changeSpace(lower_bound_par, ENERGY_SPACE);
  lower_bound_par_read = true;
}catch (int& e){
  cerr << "Cannot read lower bounds from " << lower_bound_file << endl;
  exit( 1 );
}
  }
}

void setup_bounds(){

      if(upper_bound_par_read){
      	predictor->param_factory->setMaximums(upper_bound_par);
      }
      if(lower_bound_par_read){
      	predictor->param_factory->setMinimums(lower_bound_par);
      }
}

void load_free_fix(){
  //Load free_fix from the same format as parameter vectors!
  vector< double > tmp_vector;
  par_init.getRawPars(tmp_vector);
  int num_indicators = tmp_vector.size();
  vector <bool> indicator_bool(num_indicators, true);
  if( !free_fix_indicator_filename.empty() )
  {
      ExprPar param_ff;
      try{
        param_ff = param_factory->load( free_fix_indicator_filename );
      }catch (int& e){
        cerr << "Could not parse/read the free_fix file " << free_fix_indicator_filename << endl;
        exit(1);
      }
      vector < double > tmp_ff;
      param_ff.getRawPars(tmp_ff);
      indicator_bool.clear();
      for(vector<double>::iterator iter = tmp_ff.begin();iter != tmp_ff.end();++iter){
        double one_val = *iter;
        if( -0.00000001 < one_val && 0.0000001 > one_val){ indicator_bool.push_back(false); }
        else if (0.9999999 < one_val && 1.0000001 > one_val){ indicator_bool.push_back(true);}
        else{ ASSERT_MESSAGE(false,"Illegal value in indicator_bool file");}
      }
  }
}

void train_the_predictor(){
  // random number generator
  gsl_rng* rng;
  gsl_rng_env_setup();
  const gsl_rng_type * T = gsl_rng_default;     // create rng type
  rng = gsl_rng_alloc( T );
  gsl_rng_set( rng, initialSeed );                // set the seed equal to simulTime(0)

  // model fitting

  predictor->train( par_init, rng );

  gsl_rng_free( rng );
  // print the training results
  ExprPar par = predictor->getPar();
  if( par_out_stream){
      par.print( par_out_stream, motifNames, coopMat );
      par_out_stream.close();
  }
  cout << "Estimated values of parameters:" << endl;
  par.print( cout, motifNames, coopMat );
  cout << "Performance = " << setprecision( 5 ) << ( ( ExprPredictor::objOption == SSE || ExprPredictor::objOption == PGP ) ? predictor->getObj() : -predictor->getObj() ) << endl;

}

*/


void free_fix_from_par(const ExprPar& param_ff, vector<bool>& indicator_bool){

      vector < double > tmp_ff;
      param_ff.getRawPars(tmp_ff);
      indicator_bool.clear();
      for(vector<double>::iterator iter = tmp_ff.begin();iter != tmp_ff.end();++iter){
        double one_val = *iter;
        if( -0.00000001 < one_val && 0.0000001 > one_val){ indicator_bool.push_back(false); }
        else if (0.9999999 < one_val && 1.0000001 > one_val){ indicator_bool.push_back(true);}
        else{ ASSERT_MESSAGE(false,"Illegal value in indicator_bool file");}
      }
}
