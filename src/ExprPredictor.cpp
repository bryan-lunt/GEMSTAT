#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include <nlopt.hpp>

#include "ExprPredictor.h"
#include "ExprPar.h"

double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data);

ExprFunc::ExprFunc( const vector< Motif >& _motifs, const FactorIntFunc* _intFunc, const vector< bool >& _actIndicators, int _maxContact, const vector< bool >& _repIndicators, const IntMatrix& _repressionMat, double _repressionDistThr, const ExprPar& _par ) : motifs( _motifs ), intFunc( _intFunc ), actIndicators( _actIndicators ), maxContact( _maxContact ), repIndicators( _repIndicators ), repressionMat( _repressionMat ), repressionDistThr( _repressionDistThr )
{
    par = _par;

    int nFactors = par.nFactors();
    assert( motifs.size() == nFactors );
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );

}


double ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num )
{
    bindingWts.clear(); boundaries.clear();

    if( !one_qbtm_per_crm )
        seq_num = 0;
    // store the sequence
    int n = _sites.size();
    sites = SiteVec(_sites);
    sites.insert( sites.begin(), Site() );        // start with a pseudo-site at position 0
    boundaries.push_back( 0 );
    double range = max( intFunc->getMaxDist(), repressionDistThr );
    for ( int i = 1; i <= n; i++ )
    {
        int j;
        for ( j = i - 1; j >= 1; j-- )
        {
            if ( ( sites[i].start - sites[j].start ) > range ) break;
        }
        int boundary = j;
        boundaries.push_back( boundary );
    }

    // compute the Boltzman weights of binding for all sites
    bindingWts.push_back( 1.0 );
    for ( int i = 1; i <= n; i++ )
    {
        bindingWts.push_back( par.maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].prior_probability * sites[i].wtRatio );
    }


    double Z_O = compPartFunc_RatesO();
    double Z_A = compPartFunc_RatesA();
    double Z_B = compPartFunc_RatesB();
    double Z_AB = compPartFunc_RatesAB();

    double Z_total = Z_O + Z_A + Z_B + Z_AB;

    double prob_A = (Z_A + Z_AB) / Z_total;
    double prob_B = (Z_B + Z_AB) / Z_total;



    //TODO: Handle K_max idea

    return (prob_A*prob_B*par.pis[0])/(prob_A + prob_B*par.pis[0]);


}

ModelType ExprFunc::modelOption = QUENCHING;

double ExprFunc::compPartFuncOff() const
{
    if ( modelOption == CHRMOD_UNLIMITED || modelOption == CHRMOD_LIMITED ) return compPartFuncOffChrMod();

    int n = sites.size() - 1;
    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];

        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }

        Z[i] = bindingWts[ i ] * sum;

        Zt[i] = Z[i] + Zt[i - 1];
    }

    // the partition function
    // 	double Z_bind = 1;
    // 	for ( int i = 0; i < sites.size(); i++ ) {
    // 		Z_bind += Z[ i ];
    // 	}
    return Zt[n];
}


double ExprFunc::compPartFuncOffChrMod() const
{
    int n = sites.size()- 1;

    // initialization
    vector< double > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< double > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        double sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            double dist = sites[i].start - sites[j].start;

            // sum for Z0
            sum0 += compFactorInt( sites[i], sites[j] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j];

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] )
            {
                sum1 += compFactorInt( sites[i], sites[j] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];
}


double ExprFunc::compPartFuncOn() const
{
    if ( modelOption == DIRECT ) return compPartFuncOnDirect();
    if ( modelOption == QUENCHING ) return compPartFuncOnQuenching();
    if ( modelOption == CHRMOD_UNLIMITED) return compPartFuncOnChrMod_Unlimited();
    if ( modelOption == CHRMOD_LIMITED ) return compPartFuncOnChrMod_Limited();
//TODO: A compiler warning is generated here. Shouldn't there be some defensive coding?
}


double ExprFunc::compPartFuncOnDirect() const
{
    int n = sites.size() - 1;

    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }
        //Z[i] = bindingWts[ i ] * par.txpEffects[ sites[i].factorIdx ] * sum;
        if( actIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = bindingWts[ i ] * par.txpEffects[ sites[ i ].factorIdx ] * sum;
            //cout << "1: " << par.txpEffects[ sites[ i ].factorIdx ] << endl;
        }
        if( repIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = bindingWts[ i ] * par.repEffects[ sites[ i ].factorIdx ] * sum;
            //cout << "2: " << par.repEffects[ sites[ i ].factorIdx ] << endl;
        }
        //cout << "DEBUG 0: " << sum << "\t" << Zt[ i - 1] << endl;
        Zt[i] = Z[i] + Zt[i - 1];
        /*if( actIndicators[ sites[ i ].factorIdx ] )
            cout << "DEBUG 1: " << Zt[i] << "\t" << bindingWts[i]*par.txpEffects[sites[i].factorIdx]*(Zt[ i - 1] + 1) << endl;
        if( repIndicators[ sites[ i ].factorIdx ] )
            cout << "DEBUG 2: " << Zt[i] << "\t" << bindingWts[i]*par.repEffects[sites[i].factorIdx]*(Zt[ i - 1] + 1) << endl;*/
    }

    return Zt[n];
}

double ExprFunc::compPartFunc_RatesO() const
{
    int n = sites.size() - 1;
    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];

        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }

        Z[i] = bindingWts[ i ] * sum;

        Zt[i] = Z[i] + Zt[i - 1];
    }

    // the partition function
    // 	double Z_bind = 1;
    // 	for ( int i = 0; i < sites.size(); i++ ) {
    // 		Z_bind += Z[ i ];
    // 	}
    return Zt[n];
}

double ExprFunc::compPartFunc_RatesA() const
{
    int n = sites.size() - 1;

    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * par.txpEffects[ sites[ i ].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    return Zt[n];
}

double ExprFunc::compPartFunc_RatesB() const
{
    int n = sites.size() - 1;

    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * par.repEffects[ sites[ i ].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    return Zt[n];
}

double ExprFunc::compPartFunc_RatesAB() const
{
    int n = sites.size() - 1;

    // initialization
    vector< double > Z( n + 1 );
    Z[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ i ], sites[ j ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * par.txpEffects[ sites[ i ].factorIdx ] * par.repEffects[ sites[ i ].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    return Zt[n];
}


double ExprFunc::compPartFuncOnQuenching() const
{
    int n = sites.size() - 1;
    int N0 = maxContact;
    Matrix Z1(n+1, N0+1);
    Matrix Z0(n+1, N0+1);

    // k = 0
    for ( int i = 0; i <= n; i++ )
    {
        double sum1 = 1, sum0 = 0;
        for ( int j = 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            bool R = testRepression( sites[j], sites[i] );
            double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1.getElement(j,0) + Z0.getElement(j,0) );
            sum1 += ( 1 - R )* term;
            sum0 += R * term;
        }
	Z1.setElement(i,0, bindingWts[i] * sum1);
	Z0.setElement(i,0, bindingWts[i] * sum0);
    }

    // k >= 1
    for ( int k = 1; k <= N0; k++ )
    {
        for ( int i = 0; i <= n; i++ )
        {
            if ( i < k )
            {
                Z1.setElement(i,k,0.0);
                Z0.setElement(i,k,0.0);
                continue;
            }
            double sum1 = 0, sum0 = 0;
            for ( int j = 1; j < i; j++ )
            {
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
                bool R = testRepression( sites[j], sites[i] );
                double effect = actIndicators[sites[j].factorIdx] * ( 1 - testRepression( sites[i], sites[j] ) ) * Z1.getElement(j,k-1) * par.txpEffects[sites[j].factorIdx];
                double term = compFactorInt( sites[ i ], sites[ j ] ) * ( Z1.getElement(j,k) + Z0.getElement(j,k) + effect );
                sum1 += ( 1 - R )* term;
                sum0 += R * term;
            }
            Z1.setElement(i,k,bindingWts[i] * sum1);
            Z0.setElement(i,k,bindingWts[i] * sum0);
        }
    }

    //     for ( int i = 1; i <= n; i++ ) {
    //         for ( int k = 0; k <= N0; k++ ) {
    //             cout << "Z1(" << i << ", " << k << ") = " << Z1[i][k] << "\t";
    //             cout << "Z0(" << i << ", " << k << ") = " << Z0[i][k] << endl;
    //         }
    //         cout << endl;
    //     }

    // the partition function
    double Z_on = 1;
    for ( int i = 1; i <= n; i++ )
    {
        for ( int k = 0; k <= N0; k++ )
        {
            double term = Z1.getElement(i,k) + Z0.getElement(i,k);
            Z_on += term;
        }
        for ( int k = 0; k <= N0 - 1; k++ )
        {
            Z_on += actIndicators[sites[i].factorIdx] * Z1.getElement(i,k) * par.txpEffects[sites[i].factorIdx];
        }
    }
    return Z_on;
}


double ExprFunc::compPartFuncOnChrMod_Unlimited() const
{
    int n = sites.size()- 1;

    // initialization
    vector< double > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< double > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< double > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        double sum = Zt[boundaries[i]];
        double sum0 = sum, sum1 = sum;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            double dist = sites[i].start - sites[j].start;
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;

            // sum for Z0
            sum0 += compFactorInt( sites[i], sites[j] ) * Z0[j];
            if ( dist > repressionDistThr ) sum0 += Z1[j];

            // sum for Z1
            if ( repIndicators[ sites[i].factorIdx ] )
            {
                sum1 += compFactorInt( sites[i], sites[j] ) * Z1[j];
                if ( dist > repressionDistThr ) sum1 += Z0[j];
            }
        }
        Z0[i] = bindingWts[i] * par.txpEffects[ sites[i].factorIdx ] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }

    // the partition function
    return Zt[n];
}


double ExprFunc::compPartFuncOnChrMod_Limited() const
{
    int n = sites.size()- 1;

    // initialization
    int N0 = maxContact;
    Matrix Z0( n + 1, N0 + 1 );
    Matrix Z1( n + 1, N0 + 1 );
    Matrix Zt( n + 1, N0 + 1 );
    Z0.setElement(0,0,0.0);
    Z1.setElement(0,0,0.0);
    Zt.setElement(0,0,1.0);
    for ( int k = 1; k <= N0; k++ )
    {
        Z0.setElement(0,k,0.0);
        Z1.setElement(0,k,0.0);
        Zt.setElement(0,k,0.0);
    }

    // recurrence
    for ( int k = 0; k <= N0; k++ )
    {
        for ( int i = 1; i <= n; i++ )
        {
            //             cout << "k = " << k << " i = " << i << endl;
            double sum0 = Zt.getElement(boundaries[i],k);
	    double sum0A = k > 0 ? Zt.getElement(boundaries[i],k-1) : 0.0;
	    double sum1 = sum0;

            for ( int j = boundaries[i] + 1; j < i; j++ )
            {
                double dist = sites[i].start - sites[j].start;
                if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;

                // sum for Z0
                sum0 += compFactorInt( sites[i], sites[j] ) * Z0.getElement(j,k);
                sum0A += k > 0 ? compFactorInt( sites[i], sites[j] ) * Z0.getElement(j,k-1) : 0;
                if ( dist > repressionDistThr )
                {
                    sum0 += Z1.getElement(j,k);
                    sum0A += k > 0 ? Z1.getElement(j,k-1) : 0;
                }

                // sum for Z1
                if ( repIndicators[ sites[i].factorIdx ] )
                {
                    sum1 += compFactorInt( sites[i], sites[j] ) * Z1.getElement(j,k);
                    if ( dist > repressionDistThr ) sum1 += Z0.getElement(j,k);
                }
            }
            Z0.setElement(i,k,bindingWts[i] * sum0);
            if ( actIndicators[sites[i].factorIdx] ) Z0(i,k) += k > 0 ? bindingWts[i] * par.txpEffects[sites[i].factorIdx] * sum0A : 0;
            if ( repIndicators[ sites[i].factorIdx ] ) Z1(i,k) = bindingWts[i] * par.repEffects[ sites[i].factorIdx ] * sum1;
            else Z1.setElement(i,k,0.0);
            Zt.setElement(i,k,Z0.getElement(i,k) + Z1.getElement(i,k) + Zt.getElement(i - 1,k));
            //             cout << "i = " << i << " k = " << k << " Z0 = " << Z0[i][k] << " Z1 = " << Z1[i][k] << " Zt = " << Zt[i][k] << endl;
        }
    }

    // the partition function
    //     cout << "Zt[n] = " << Zt[n] << endl;
    return sum( Zt.getRow(n) );//And we end up with a vector anyway. See about fixing this.
}


double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
    // 	assert( !siteOverlap( a, b, motifs ) );
    double maxInt = par.factorIntMat( a.factorIdx, b.factorIdx );
    double dist = abs( a.start - b.start );
    bool orientation = ( a.strand == b.strand );
    return intFunc->compFactorInt( maxInt, dist, orientation );
}


bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
    // 	assert( !siteOverlap( a, b, motifs ) );

    double dist = abs( a.start - b.start );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}


ExprPredictor::ExprPredictor( const vector <Sequence>& _seqs, const vector< SiteVec >& _seqSites, const vector < SiteVec >& _r_seqSites, const vector< int >& _seqLengths, const vector <int>& _r_seqLengths, const Matrix& _exprData, const vector< Motif >& _motifs, const Matrix& _factorExprData, const ExprModel& _expr_model,
		const vector < bool >& _indicator_bool, const vector <string>& _motifNames, const vector < int >& _axis_start, const vector < int >& _axis_end, const vector < double >& _axis_wts ) : seqs(_seqs), seqSites( _seqSites ), r_seqSites( _r_seqSites ), seqLengths( _seqLengths ), r_seqLengths( _r_seqLengths ), exprData( _exprData ), factorExprData( _factorExprData ),
	expr_model( _expr_model),
	indicator_bool ( _indicator_bool ), motifNames ( _motifNames ), axis_start ( _axis_start ), axis_end( _axis_end ), axis_wts( _axis_wts )
{
    //TODO: Move appropriate lines from this block to the ExprModel class.
    assert( exprData.nRows() == nSeqs() );
    assert( factorExprData.nRows() == nFactors() && factorExprData.nCols() == nConds() );
    assert( expr_model.coopMat.isSquare() && expr_model.coopMat.isSymmetric() && expr_model.coopMat.nRows() == nFactors() );
    assert( expr_model.actIndicators.size() == nFactors() );
    assert( expr_model.maxContact > 0 );
    assert( expr_model.repIndicators.size() == nFactors() );
    assert( expr_model.repressionMat.isSquare() && expr_model.repressionMat.nRows() == nFactors() );
    assert( expr_model.repressionDistThr >= 0 );

    //gene_crm_fout.open( "gene_crm_fout.txt" );

    // set the model option for ExprPar and ExprFunc
    ExprPar::modelOption = expr_model.modelOption;//TODO: Remove both of these.
    ExprFunc::modelOption = expr_model.modelOption;

    // set the values of the parameter range according to the model option
    if ( expr_model.modelOption != LOGISTIC && expr_model.modelOption != DIRECT )
    {
        //ExprPar::min_effect_Thermo = 0.99;
        //ExprPar::min_interaction = 0.99;
    }

    // set the option of parameter estimation
    ExprPar::estBindingOption = estBindingOption;

    //expr_model was already initialized. Setup the parameter factory.
    param_factory = new ParFactory(expr_model, nSeqs(), _indicator_bool);

    //TODO: Move this to the front-end or something?
    //Maybe make it have a default SSE score objective, but anything else gets specified in the front-end.
    switch(objOption){
      case CORR:
        trainingObjective = new AvgCorrObjFunc();
        break;
      case PGP:
        trainingObjective = new PGPObjFunc();
        break;
      case CROSS_CORR:
        trainingObjective = new AvgCrossCorrObjFunc(ExprPredictor::maxShift, ExprPredictor::shiftPenalty);
        break;
      case LOGISTIC_REGRESSION:
        trainingObjective = new LogisticRegressionObjFunc();
        break;
      case SSE:
      default:
        trainingObjective = new RMSEObjFunc();
        break;
    }

    /* DEBUG
    cout << setprecision(10);
    ExprPar foo = param_factory->createDefaultMinMax(true);
    cout << " MAXIMUMS " << endl;
    printPar(foo);
    foo = param_factory->createDefaultMinMax(false);
    cout << " MINIMUMS " << endl;
    printPar(foo);
    */
}

ExprPredictor::~ExprPredictor()
{
  delete param_factory;
  delete trainingObjective;
}

double ExprPredictor::objFunc( const ExprPar& par )
{
    double objective_value = evalObjective( par );

    return objective_value;
}


int ExprPredictor::train( const ExprPar& par_init )
{
    par_model = par_init;

    cout << "*** Diagnostic printing BEFORE adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << objFunc( par_model ) << endl;
    cout << "*******************************************" << endl << endl;

    if ( nAlternations > 0 && ExprPar::searchOption == CONSTRAINED ){
      ExprPar tmp_constrained = param_factory->changeSpace(par_model,CONSTRAINED_SPACE);
      par_model = param_factory->changeSpace(par_model,PROB_SPACE);
    }
    obj_model = objFunc( par_model );

    cout << "*** Diagnostic printing AFTER adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << objFunc( par_model ) << endl;
    cout << "*******************************************" << endl << endl;

    if ( nAlternations == 0 ) return 0;

    // alternate between two different methods
    ExprPar par_result;
    double obj_result;
    for ( int i = 0; i < nAlternations; i++ )
    {
        simplex_minimize( par_result, obj_result );
        par_model = par_result;

        gradient_minimize( par_result, obj_result );
        par_model = par_result;
    }

    #ifdef BETAOPTBROKEN
    optimize_beta( par_model, obj_result );
    #endif

    // commit the parameters and the value of the objective function
    //par_model = par_result;
    obj_model = obj_result;

    return 0;
}


int ExprPredictor::train( const ExprPar& par_init, const gsl_rng* rng )
{
    /*
        //for random starts:
        ExprPar par_rand_start = par_init;
        par_rand_start = param_factor->randSamplePar( rng );
        train( par_rand_start );*/
    // training using the initial values
    train( par_init );

    cout << "Initial training:\tParameters = "; printPar( par_model );
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;

    // training with random starts
    ExprPar par_best = par_model;
    double obj_best = obj_model;
    for ( int i = 0; i < nRandStarts; i++ )
    {
        ExprPar par_curr = par_init;
        par_curr = param_factory->randSamplePar( rng );
        train( par_curr );
        cout << "Random start " << i + 1 << ":\tParameters = "; printPar( par_model );
        cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;
        if ( obj_model < obj_best )
        {
            par_best = par_model;
            obj_best = obj_model;
        }
    }

    // training using the best parameters so far
    if ( nRandStarts ) train( par_best );
    cout << "Final training:\tParameters = "; printPar( par_model );
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;

    //gene_crm_fout.close();

    return 0;
}


int ExprPredictor::train()
{
    // random number generator
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;     // create rng type
    rng = gsl_rng_alloc( T );
    gsl_rng_set( rng, time( 0 ) );                // set the seed equal to simulTime(0)

    // training using the default initial values with random starts
    ExprPar par_default( nFactors(), nSeqs() );
    train( par_default, rng );

    return 0;
}


int ExprPredictor::predict( const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const
{
    targetExprs.clear();

    // create site representation of the target sequence
    SiteVec targetSites;
    #ifdef REANNOTATE_EACH_PREDICTION
    SeqAnnotator ann( expr_model.motifs, par_model.energyThrFactors );
    ann.annot( seqs[ seq_num ], targetSites );
    #else
    targetSites = targetSites_;
    #endif

    // predict the expression
    ExprFunc* func = createExprFunc( par_model );
    for ( int j = 0; j < nConds(); j++ )
    {
        vector< double > concs = factorExprData.getCol( j );
        double predicted = func->predictExpr( targetSites, targetSeqLength, concs, seq_num );
        targetExprs.push_back( predicted );
    }

    delete func;
    return 0;
}

int ExprPredictor::estBindingOption = 1;          // 1. estimate binding parameters; 0. not estimate binding parameters
ObjType ExprPredictor::objOption = SSE;

int ExprPredictor::maxShift = 5;
double ExprPredictor::shiftPenalty = 0.8;

int ExprPredictor::nAlternations = 4;
int ExprPredictor::nRandStarts = 5;
double ExprPredictor::min_delta_f_SSE = 1.0E-8;
double ExprPredictor::min_delta_f_Corr = 1.0E-8;
double ExprPredictor::min_delta_f_CrossCorr = 1.0E-8;
double ExprPredictor::min_delta_f_PGP = 1.0E-8;
int ExprPredictor::nSimplexIters = 200;
int ExprPredictor::nGradientIters = 50;




bool ExprPredictor::testPar( const ExprPar& par ) const
{
    // test binding weights
    if ( estBindingOption )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( par.maxBindingWts[i] < ExprPar::min_weight * ( 1.0 + ExprPar::delta ) || par.maxBindingWts[i] > ExprPar::max_weight * ( 1.0 - ExprPar::delta ) )
                return false;
        }
    }

    //cout << "dbg1" << endl;

    // test the interaction matrix
    if ( expr_model.modelOption != LOGISTIC )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            for ( int j = 0; j <= i; j++ )
            {
                if ( par.factorIntMat( i, j ) < ExprPar::min_interaction * ( 1.0 + ExprPar::delta ) || par.factorIntMat( i, j ) > ExprPar::max_interaction * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    // test the transcriptional effects
    for ( int i = 0; i < nFactors(); i++ )
    {
        if ( expr_model.modelOption == LOGISTIC )
        {
            if ( par.txpEffects[i] < ExprPar::min_effect_Logistic + ExprPar::delta || par.txpEffects[i] > ExprPar::max_effect_Logistic - ExprPar::delta )
                return false;
        }
        else if ( expr_model.modelOption == DIRECT )
        {
            if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
        }
        else
        {
            if ( expr_model.actIndicators[i] )
            {
                if ( par.txpEffects[i] < ExprPar::min_effect_Thermo * ( 1.0 + ExprPar::delta ) || par.txpEffects[i] > ExprPar::max_effect_Thermo * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    //cout << "dbg2" << endl;
    // test the repression effects
    if ( expr_model.modelOption == CHRMOD_UNLIMITED || expr_model.modelOption == CHRMOD_LIMITED || expr_model.modelOption == DIRECT )
    {
        for ( int i = 0; i < nFactors(); i++ )
        {
            if ( expr_model.repIndicators[i] )
            {
                if ( par.repEffects[i] < ExprPar::min_repression * ( 1.0 + ExprPar::delta ) || par.repEffects[i] > ExprPar::max_repression * ( 1.0 - ExprPar::delta ) )
                    return false;
            }
        }
    }

    //cout << "dbg3" << endl;
    // test the basal transcription
    for( int _i = 0; _i < par.basalTxps.size(); _i++ )
    {
        if ( expr_model.modelOption == LOGISTIC )
        {
            if ( par.basalTxps[ _i ] < ExprPar::min_basal_Logistic + ExprPar::delta || par.basalTxps[ _i ] > ExprPar::max_basal_Logistic - ExprPar::delta )
                return false;
        }
        else
        {
            if ( par.basalTxps[ _i ] < ExprPar::min_basal_Thermo * ( 1.0 + ExprPar::delta ) || par.basalTxps[ _i ] > ExprPar::max_basal_Thermo * ( 1.0 - ExprPar::delta ) )
                return false;
        }
    }
    //test the pi values

    for( int _i = 0; _i < par.pis.size(); _i++ )
    {
        if( par.pis[ _i ] < 0 || par.pis[ _i ] > ExprPar::max_pi * ( 1.0 - ExprPar::delta ))
        {
            return false;
        }
    }
    //cout << "dbg4" << endl;
    //test the beta values
    for( int _i = 0; _i < par.betas.size(); _i++ )
    {
        if ( par.betas[ _i ] < ExprPar::min_beta * ( 1.0 + ExprPar::delta ) || par.betas[ _i ] > ExprPar::max_beta * ( 1.0 - ExprPar::delta ) )
            return false;
    }
    //cout << "dbg5" << endl;
    //adjust the energyThrFactors
    for( int i = 0; i < nFactors(); i++ )
    {
        if( par.energyThrFactors[ i ] < ExprPar::min_energyThrFactors * ( 1.0 + ExprPar::delta ) ) return false;
        if( par.energyThrFactors[ i ] > ExprPar::max_energyThrFactors * ( 1.0 - ExprPar::delta ) ) return false;
    }
    return true;
}

void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
    cout.precision( 8 );
    //     cout.width( 8 );

    // print binding weights
    cout << "MAXBIND : " << par.maxBindingWts << endl;
    cout << "INTER : " ;
    // print the interaction matrix
    for ( int i = 0; i < nFactors(); i++ )
    {
        for ( int j = 0; j <= i; j++ )
        {
           cout << par.factorIntMat( i, j ) << "\t";
        }
    }
    cout << endl;

    // print the transcriptional effects
    cout << "TXP : " << par.txpEffects << endl;

    // print the repression effects
    cout << "REP : " << par.repEffects << endl;

    // print the basal transcriptions
    cout << "BASAL : " << par.basalTxps << endl;

    //print the pi values
    cout << "PIS : " << par.pis << endl;

    //print the beta values
    cout << "BETAS : " << par.betas << endl;
    //assert( par.betas.size() == nSeqs() );

    cout << "THRESH : " << par.energyThrFactors << endl;
    cout << flush;
}


ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par ) const
{
	//TODO: just make it take an expr_model as a parameter.
	//Also, it should use a factory to create the ExprFunc, but this is already close to a factory method.

    //TODO: This makes ExprPredictor act as the factory for ExprFunc objects.
    //Since a par factory is not needed to go between ENERGY_SPACE and PROB_SPACE, this could get moved into the ExprFunc constructor. That would be better for inheritance.
    ExprPar parToPass = param_factory->changeSpace(par, expr_model.modelOption == LOGISTIC ? ENERGY_SPACE : PROB_SPACE );

    return new ExprFunc( expr_model.motifs, expr_model.intFunc, expr_model.actIndicators, expr_model.maxContact, expr_model.repIndicators, expr_model.repressionMat, expr_model.repressionDistThr, parToPass );
}


int indices_of_crm_in_gene[] =
{
};

double ExprPredictor::evalObjective( const ExprPar& par )
{

    vector< int > seqLengths( seqs.size() );
    for( int i = 0; i < seqs.size(); i++ ){
      seqLengths[i] = seqs[i].size();
    }

    ExprFunc* func = createExprFunc( par );

    vector< SiteVec > seqSites( seqs.size() ); //
    #ifdef REANNOTATE_EACH_PREDICTION
    SeqAnnotator ann( expr_model.motifs, par.energyThrFactors );
    for ( int i = 0; i < seqs.size(); i++ ) {
       	ann.annot( seqs[ i ], seqSites[ i ] );
    }
    #else
    seqSites = this->seqSites;
    #endif

    vector<vector<double> > ground_truths;
    vector<vector<double> > predictions;

    //Create predictions for every sequence and condition
    for ( int i = 0; i < nSeqs(); i++ ) {
        ground_truths.push_back(exprData.getRow(i));

        vector< double > predictedExprs;
        for ( int j = 0; j < nConds(); j++ ) {
        		double predicted = -1;
            	vector< double > concs = factorExprData.getCol( j );
            	predicted = func->predictExpr( seqSites[ i ], seqLengths[i], concs, i );
            // predicted expression for the i-th sequence at the j-th condition
            predictedExprs.push_back( predicted );
        }
        predictions.push_back(predictedExprs);
    }

    //Evaluate the objective function on that.
    double ret_val = trainingObjective->eval(ground_truths, predictions, &par);
    delete func;
    return ret_val;

}

int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector < double > pars;

    //ExprPar tmp_par_model = param_factory->changeSpace(par_model, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar tmp_par_model = param_factory->changeSpace(par_model, ENERGY_SPACE);
    param_factory->separateParams(tmp_par_model, free_pars, fix_pars, indicator_bool );

    pars.clear();
    pars = free_pars;

    //SIMPLEX MINIMIZATION with NLOPT
    nlopt::opt optimizer(nlopt::LN_NELDERMEAD, pars.size());
    optimizer.set_min_objective(nlopt_obj_func, this);
    //optimizer.set_maxmaxeval(nSimplexIters);
    optimizer.set_initial_step(1.0);//TODO: enforce simplex staring size.
    optimizer.set_maxeval(nSimplexIters);//TODO: enforce nSimplexIters

    if(ExprPar::searchOption == CONSTRAINED){
      vector<double> free_mins;
      vector<double> fix_mins;

      param_factory->separateParams(param_factory->getMinimums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_lower_bounds(free_mins);

      param_factory->separateParams(param_factory->getMaximums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_upper_bounds(free_mins);
    }

    nlopt::result result = optimizer.optimize(free_pars, obj_result);
    obj_result = optimizer.last_optimum_value();
    //Done Minimizing

    param_factory->joinParams(free_pars, fix_pars, pars, indicator_bool);
    //tmp_par_model = param_factory->create_expr_par(pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    tmp_par_model = param_factory->create_expr_par(pars, ENERGY_SPACE);
    par_result = param_factory->changeSpace(tmp_par_model, PROB_SPACE);

    printPar( par_result );

    return 0;
}


int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector< double > pars;
    //cout << "DEBUG: in getFreePars()" << endl;
    //par_model.getFreePars( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators );
    //cout << "DEBUG: out getFreePars()" << endl;
    //ExprPar tmp_par_model = param_factory->changeSpace(par_model, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar tmp_par_model = param_factory->changeSpace(par_model, ENERGY_SPACE);

    param_factory->separateParams(tmp_par_model, free_pars, fix_pars, indicator_bool );

    pars.clear();
    pars = free_pars;

    //GRADIENT MINIMIZATION with NLOPT
    nlopt::opt optimizer(nlopt::LD_LBFGS, pars.size());
    optimizer.set_min_objective(nlopt_obj_func, this);

    //TODO: Move this to a nice lookup table or something.
    //Set the stopping criterion
    double ftol;
    switch(objOption){
      SSE:
        ftol = min_delta_f_SSE;
        break;
      CORR:
        ftol = min_delta_f_Corr;
        break;
      CROSS_CORR:
        ftol = min_delta_f_CrossCorr;
        break;
      PGP:
        ftol = min_delta_f_PGP;
        break;
      default:
        ftol = 1e-5;
        break;
    }
    optimizer.set_ftol_abs(ftol);

    if(ExprPar::searchOption == CONSTRAINED){
      vector<double> free_mins;
      vector<double> fix_mins;

      param_factory->separateParams(param_factory->getMinimums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_lower_bounds(free_mins);

      param_factory->separateParams(param_factory->getMaximums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_upper_bounds(free_mins);
    }


    //TODO: enforce nGradientIters
    optimizer.set_maxeval(nGradientIters);
    try{
      nlopt::result result = optimizer.optimize(free_pars, obj_result);
      obj_result = optimizer.last_optimum_value();
    }catch(std::runtime_error){
      cerr << "There was an exception in the gradient descent!" << endl;
    }

    //Done Minimizing
    //pars now contains the optimal parameters

    param_factory->joinParams(free_pars, fix_pars, pars, indicator_bool);
    //tmp_par_model = param_factory->create_expr_par(pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    tmp_par_model = param_factory->create_expr_par(pars, ENERGY_SPACE);

    par_result = param_factory->changeSpace(tmp_par_model, PROB_SPACE);
    cout << "DEBUG" << endl;
    return 0;
}


double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data){
        gsl_vector *xv = vector2gsl(x); //TODO: Ugly, remove (Make all objective functions use native STL vectors)
        double objective = gsl_obj_f(xv, f_data);

        if(!grad.empty()){
                gsl_vector *dxv = vector2gsl(grad);

                gsl_obj_df(xv,f_data,dxv);

                for(int i = 0;i< grad.size();i++){
                        grad[i] = dxv->data[i];
                }
                cerr << " obj called, derivative : " << grad << endl;
                gsl_vector_free(dxv);
        }

        cerr << " obj called " << objective << endl;

        gsl_vector_free(xv);
        return objective;
}

double gsl_obj_f( const gsl_vector* v, void* params )
{
    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;

    // parse the variables (parameters to be optimized)
    //     vector< double > expv;
    //     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
    vector <double> temp_free_pars = gsl2vector(v);
    vector < double > all_pars;

    predictor->param_factory->joinParams(temp_free_pars, predictor->fix_pars, all_pars, predictor->indicator_bool);
    //ExprPar par = predictor->param_factory->create_expr_par(all_pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar par = predictor->param_factory->create_expr_par(all_pars, ENERGY_SPACE);
    par = predictor->param_factory->changeSpace(par, PROB_SPACE); //TODO: WTF? This shouldn't be required because it's done in the createExprFunc method. Stack corruption or something?


    // call the ExprPredictor object to evaluate the objective function
    double obj = predictor->objFunc( par );
    return obj;
}

void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
    double step = 1.0E-3;
    numeric_deriv( grad, gsl_obj_f, v, params, step );
}


void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f( v, params );
    gsl_obj_df( v, params, grad );
}
