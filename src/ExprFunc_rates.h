#ifndef EXPR_FUNC_RATES_H
#define EXPR_FUNC_RATES_H

#include "ExprFunc.h"


class Rates_ExprFunc : public ExprFunc {
  public:
      // constructors
      Rates_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
      double predictExpr( const vector< double >& factorConcs );
  protected:
    // compute the partition function when the BTM is bound
    double compPartFuncAB() const;
    double compPartFuncA() const;
    double compPartFuncB() const;
    double compPartFuncO() const;

};


#endif //EXPR_FUNC_RATES_H
