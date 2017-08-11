#include "ExprPredictor.h"
#include "ExprModel.h"

#include "Pioneers.h"

PioneerExprModel::PioneerExprModel( ModelType _modelOption, bool _one_qbtm_per_crm, vector< Motif>& _motifs, FactorIntFunc* _intFunc, int _maxContact, IntMatrix& _coopMat, vector< bool >& _actIndicators, vector< bool>& _repIndicators, vector< bool>& _pioneer_indicators, IntMatrix& _repressionMat, double _repressionDistThr ) : ExprModel( _modelOption, _one_qbtm_per_crm, _motifs, _intFunc, _maxContact, _coopMat, _actIndicators, _repIndicators, _repressionMat, _repressionDistThr ), pioneer_indicators( _pioneer_indicators)
{
}

Pioneer_ExprFunc::Pioneer_ExprFunc( ExprFunc *to_wrap, const PioneerExprModel& in_model, const SiteVec& sites_, const int seq_length, const int seq_num, const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : ExprFunc( sites_, seq_length, seq_num, _motifs,  _intFunc, _actIndicators, _maxContact, _repIndicators, _repressionMat, _repressionDistThr, _par), expr_model(in_model)
{
    wrapped_expr_func = to_wrap;

    pioneer_dist_thresh = 100;

    //cerr << "gamma: " << par.extra_params[0] << endl;

    //cerr << "blah" << endl;//TODO: Segfaults somewhere in here.
    //Setup a mapping that allows a site to quickly find who its pioneers are.
    // site_pioneers[i] will be a vector listing all the pioneers who exert an influence on site[i];
    site_pioneers.clear();
    site_pioneers.reserve(sites.size());    //allocate all the necessary memory in advance.
    site_pioneers.push_back(vector<int>(0));    //no pioneers affect the beginning pseudo-site
    for(int i = 1;i<=n_sites;i++){
        site_pioneers.push_back(vector<int>(0));    //the vector that stores those pioneers that affect site[i]
        //pioneers are not affected by pioneers

        if(expr_model.pioneer_indicators[ sites[i].factorIdx ]){
            //If this site is a pioneer, consider it immune to pioneer activity. //TODO: Possibly stupid.
            continue; }

        for(int j = i-1;j>0;j--){
            //search backwards over previous sites to find interacting pioneers up to a distance of pioneer_dist_thresh
            //Ignoring the beginning pseudo-site
            if( sites[i].start - sites[j].start > pioneer_dist_thresh ){  break; }
            if( expr_model.pioneer_indicators[ sites[j].factorIdx ] ){
                site_pioneers[i].push_back(j);
                //cerr << "site " << j << " shall affect site " << i << endl;
            }
        }
        for(int j = i+1;j<=n_sites;j++){
            //Ignoring the terminal pseudo-site
            //search forwards over subsequent sites to find interacting pioneers up to a distance of pioneer_dist_thresh
            if( sites[j].start - sites[i].start > pioneer_dist_thresh ){  break; }
            if( expr_model.pioneer_indicators[ sites[j].factorIdx ] ){
                site_pioneers[i].push_back(j);
                //cerr << "site " << j << " shall affect site " << i << endl;
            }
        }
    }
    site_pioneers.push_back(vector<int>(0));    //no pioneers affect the terminal pseudo-site
    assert(site_pioneers.size() == sites.size());   //There should be one (potentially empty) list of interacting pioneers per site.
    //cerr << "blah" << endl;
}

Pioneer_ExprFunc::~Pioneer_ExprFunc(){
    delete wrapped_expr_func;
}

double Pioneer_ExprFunc::predictExpr( const Condition& in_condition )
{
    //TODO: calculate the priors based on pioneer effects and concentrations and modify the priors of the wrapped function.

    //Estimate the strengths of pioneer factor binding/effects.
    vector< double > pioneer_binding_and_activity(sites.size(),0.0);
    vector< double > saved_priors(sites.size(),0.0);
    //Save out the priors for the sites, because we are going to alter that.
    for(int i = 0;i<sites.size();i++){
        saved_priors[i] = sites[i].prior_probability;
    }

    estimate_pioneer_activity(in_condition, pioneer_binding_and_activity);

    //Use that to alter the priors of the other binding sites.
    alter_priors_by_pioneer_activity(pioneer_binding_and_activity);


    double to_return = wrapped_expr_func->predictExpr( in_condition.concs );
    //Restor the saved prior probabilities for the sites.
    for(int i = 0;i<sites.size();i++){
        sites[i].prior_probability = saved_priors[i];
    }
    return to_return;
}

void Pioneer_ExprFunc::estimate_pioneer_activity(const Condition& in_condition, vector< double > &destination) const {
    //Just do it like the logistic function for now.
    assert(destination.size() == sites.size());
    destination[0] = 1.0;
    destination[destination.size()-1] = 1.0;
    for(int i=1;i<n_sites;i++){
        if(!expr_model.pioneer_indicators[ sites[i].factorIdx  ]){
            destination[i] = 0.0;
            continue;
        }
        double this_site_binding_weight = par.maxBindingWts[ sites[i].factorIdx ] * in_condition.concs[sites[i].factorIdx] * sites[i].prior_probability * sites[i].wtRatio;
        double this_site_binding_prob = this_site_binding_weight / ( 1.0 + this_site_binding_weight );
        destination[i] = this_site_binding_prob;
        //TODO: we got the binding probabilities, how about the strength of the pioneer activity?
    }

    //cerr << "FOOBAR " << destination << endl;
}

void Pioneer_ExprFunc::alter_priors_by_pioneer_activity(vector< double > &in_pioneer_activity) {
    double gamma = par.extra_params[0];
    //TODO: Dammit, we will use a naive approach here.
    for(int i=1;i<=n_sites;i++){
        //once again, pioneers are not affected by pioneer activity (at this stage)
        if(expr_model.pioneer_indicators[ sites[i].factorIdx ]){ continue; }

        double total_pioneer_activity_for_this_site = 0.0;
        //So, do want to take the maximum of pioneers in range?
        //The sum? Soft^sum?
        //For no, I guess the max, but that isn't very smooth.
        for(int j_index = 0;j_index<site_pioneers[i].size();j_index++){
            //TODO: FInish this.
            int j = site_pioneers[i][j_index];
            double one_pioneer_effect = calculate_pioneer_effect_on_i_from_j(i,j,in_pioneer_activity[j]);
            //cerr << "one effect : " << one_pioneer_effect << endl;//DEBUG
            total_pioneer_activity_for_this_site = max(one_pioneer_effect,total_pioneer_activity_for_this_site);

        }

        double new_prior = (1.0 - gamma) + gamma*total_pioneer_activity_for_this_site;

        /*
        if(site_pioneers[i].size() > 0){//DEBUG
            cerr << "total pioneer effect : " << total_pioneer_activity_for_this_site << endl;
        }
        */
        //cerr << "modified : " << new_prior - sites[i].prior_probability << endl;
        sites[i].prior_probability = new_prior;
    }
}

double Pioneer_ExprFunc::calculate_pioneer_effect_on_i_from_j(const int i, const int j, const double its_activity) const
{
    if(!expr_model.pioneer_indicators[ sites[j].factorIdx ]){
        cerr << "indicators : " << expr_model.pioneer_indicators << endl;
        cerr << "non pioneer site: " << sites[j] << endl;
    }
    assert(expr_model.pioneer_indicators[ sites[j].factorIdx ]);
    //cerr << "site " << i << " will feel an effect from site " << j << " of : " << its_activity << endl;
    return its_activity;
}

ExprFunc* PioneerExprModel::createNewExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const
{

    //consider filtering the sites to remove pioneer-only sites.
    ExprFunc *func_to_wrap = ExprModel::createNewExprFunc(par,sites_,seq_length,seq_num);

    ExprPar parToPass = par.my_factory->changeSpace(par, PROB_SPACE );
    Pioneer_ExprFunc *func_to_return = new Pioneer_ExprFunc(func_to_wrap, *this,
                    sites_,seq_length,seq_num,
                      this->motifs,
                      this->intFunc,
                      this->actIndicators,
                      this->maxContact,
                      this->repIndicators,
                      this->repressionMat,
                      this->repressionDistThr,
                      parToPass );

    //TODO: Wrap this and return a new PioneerExprModel.
    return func_to_return;
}
