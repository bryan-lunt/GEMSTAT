#ifndef PIONEERS_H
#define PIONEERS_H

#include "ExprModel.h"
#include "ExprFunc.h"

class PioneerExprModel : public ExprModel {
public: //TODO: Implement good accessors / mutators instead.
	PioneerExprModel( ModelType _modelOption, bool _one_qbtm_per_crm, vector< Motif>& _motifs, FactorIntFunc* _intFunc, int _maxContact, IntMatrix& _coopMat, vector< bool >& _actIndicators, vector< bool>& _repIndicators, vector< bool >& _pioneer_indicators, IntMatrix& _repressionMat, double _repressionDistThr );

    ExprFunc* createNewExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const;
    vector< bool >& pioneer_indicators;
};

class Pioneer_ExprFunc : public ExprFunc {
  public:
      Pioneer_ExprFunc( ExprFunc *to_wrap, const PioneerExprModel& in_model, const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par );// : ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par);
      ~Pioneer_ExprFunc();
      virtual double predictExpr( const Condition& in_condition );
      int pioneer_dist_thresh;
  protected:
      ExprFunc *wrapped_expr_func;
      const PioneerExprModel& expr_model;

      vector< vector< int > > site_pioneers;

  private:
      void estimate_pioneer_activity(const Condition& in_condition, vector< double > &destination) const;
      void alter_priors_by_pioneer_activity(vector< double > &in_pioneer_activity);
      double calculate_pioneer_effect_on_i_from_j(const int i, const int j, const double its_activity) const;
};

#endif
