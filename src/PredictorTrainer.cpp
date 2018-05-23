/*
 * PredictorTrainer.cpp
 *
 *  Created on: Jul 31, 2015
 *      Author: lunt
 */

#include "PredictorTrainer.h"

ObjType getObjOption( const string& objOptionStr )
{
    if ( toupperStr( objOptionStr ) == "SSE" ) return SSE;
    if ( toupperStr( objOptionStr ) == "CORR" ) return CORR;
    if ( toupperStr( objOptionStr ) == "CROSS_CORR" ) return CROSS_CORR;
    if ( toupperStr( objOptionStr ) == "PGP" ) return PGP;
    if ( toupperStr( objOptionStr ) == "LOGISTIC_REGRESSION") return LOGISTIC_REGRESSION;
    if ( toupperStr( objOptionStr ) == "PEAK_WEIGHTED") return PEAK_WEIGHTED;
    if ( toupperStr( objOptionStr ) == "WEIGHTED_SSE") return WEIGHTED_SSE;


    cerr << "objOptionStr is not a valid option of objective function" << endl;
    exit(1);
}


string getObjOptionStr( ObjType objOption )
{
    if ( objOption == SSE ) return "SSE";
    if ( objOption == CORR ) return "Corr";
    if ( objOption == CROSS_CORR ) return "Cross_Corr";
    if ( objOption == PGP ) return "PGP";
    if ( objOption == LOGISTIC_REGRESSION ) return "LOGISTIC_REGRESSION";
    if ( objOption == PEAK_WEIGHTED ) return "PEAK_WEIGHTED";
    if ( objOption == WEIGHTED_SSE ) return "WEIGHTED_SSE";

    return "Invalid";
}


string getSearchOptionStr( SearchType searchOption )
{
    if ( searchOption == UNCONSTRAINED ) return "Unconstrained";
    if ( searchOption == CONSTRAINED ) return "Constrained";

    return "Invalid";
}

//bool TrainingCheckpointer::peek_checkpoint();
bool TrainingCheckpointer::checkpoint(){
	bool final_value = false;



	if(!this->cp_active){
		return false;
	}

	//std::cerr << "TEST CHECKPOINT " << this->cp_check_n_evals_elapsed << std::endl;

	this->cp_check_n_evals_elapsed += 1;
	if(this->cp_check_n_evals_elapsed >= this->cp_check_n_evals){
		//std::cerr << "NEVALS PASSED" << std::endl;
		std::time_t nowtime = std::time(NULL);
		if(this->cp_last_checkpoint + this->cp_interval <= nowtime){

			final_value = this->checkpoint_impl();

			this->cp_last_checkpoint = std::time(NULL);//Checkpointing may be time consuming.
			this->cp_check_n_evals_elapsed = 0;
			if( final_value ){
				this->times_checkpointed += 1;
			}
		}
	}

	return final_value;
}

bool TrainingCheckpointer::load_checkpoint(){
	//#ifdef DEBUG
	std::cerr << "Uninplemented load_checkpoint in TrainingCheckpointer subclass" << endl;
	//#endif

	return false;
}
bool TrainingCheckpointer::checkpoint_impl(){
	//#ifdef DEBUG
	std::cerr << "Uninplemented checkpoint_helper in TrainingCheckpointer subclass" << endl;
	//#endif

	return false;
}

/*
PredictorTrainer::PredictorTrainer() {
	// TODO Auto-generated constructor stub

}

PredictorTrainer::~PredictorTrainer() {
	// TODO Auto-generated destructor stub
}

*/
