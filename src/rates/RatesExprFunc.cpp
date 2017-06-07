#include "RatesExprFunc.h"

double Rates_ExprFunc::predictExpr( const SiteVec& _sites, int length, const vector< double >& factorConcs, int seq_num )
{
    //initialize the sites and boundaries and whatnot.
    setupSitesAndBoundaries(_sites,length,seq_num);

    // compute the Boltzman weights of binding for all sites
    setupBindingWeights(factorConcs);


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

double Rates_ExprFunc::compPartFunc_RatesO() const
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

double Rates_ExprFunc::compPartFunc_RatesA() const
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

double Rates_ExprFunc::compPartFunc_RatesB() const
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

double Rates_ExprFunc::compPartFunc_RatesAB() const
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
