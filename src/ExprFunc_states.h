#ifndef EXPR_FUNC_STATES_H
#define EXPR_FUNC_STATES_H

#include "ExprFunc.h"

typedef struct {gemstat_dp_t O, A, AB; } states_opt_return_t;

class States_ExprFunc : public ExprFunc {
  public:
      // constructors
      States_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
      double predictExpr( const vector< double >& factorConcs );
  protected:

		states_opt_return_t compAllParts() const;
};


#endif //EXPR_FUNC_STATES_H
