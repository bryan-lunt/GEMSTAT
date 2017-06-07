#ifndef RATES_EXPR_FUNC_H
#define RATES_EXPR_FUNC_H

#include "../ExprFunc.h"


class Rates_ExprFunc : public ExprFunc {
  public:
      // constructors
      Rates_ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ExprFunc( _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par){} ;
      double predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num );
  protected:
    double compPartFunc_RatesO() const;
    double compPartFunc_RatesA() const;
    double compPartFunc_RatesB() const;
    double compPartFunc_RatesAB() const;
};


#endif
