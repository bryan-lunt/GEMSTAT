



#include "ExprFunc_rates.h"

double Rates_ExprFunc::predictExpr( const vector< double >& factorConcs )
{

    // compute the Boltzman weights of binding for all sites
    setupBindingWeights(factorConcs);

    double Z_O  = compPartFuncO();
    double Z_A  = compPartFuncA();
    double Z_B  = compPartFuncB();
    double Z_AB = compPartFuncAB();

    gemstat_dp_t Z_total = Z_O + Z_A + Z_B + Z_AB;

    gemstat_dp_t prob_A = (Z_A + Z_AB) / Z_total;
    gemstat_dp_t prob_B = (Z_B + Z_AB) / Z_total;



    //TODO: Handle K_max idea

    /** I have no idea what the pis were for in this originally */
    //return (prob_A*prob_B*par.pis[seq_num])/(prob_A + prob_B*par.pis[seq_num]);
    
    return (prob_A*prob_B)/(prob_A + prob_B);
}

double Rates_ExprFunc::compPartFuncO() const
{
    #ifdef DEBUG
      //assert(modelOption != CHRMOD_UNLIMITED && modelOption != CHRMOD_LIMITED );
    #endif

    int n = n_sites;
    // initialization
    vector< gemstat_dp_t > Z( n + 1 );
    Z[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        if( sum != sum )
        {
            cout << "DEBUG: sum nan" << "\t" << Zt[ boundaries[i] ] <<  endl;
            exit(1);
        }
        //cout << "DEBUG: sum = " << n << endl;
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
            //cout << "compFactorInt: " << compFactorInt( sites[ j ], sites[ i ] ) << "\t";
            //cout << "Z[j]: " << Z[ j ] << endl;
            #ifdef DEBUG
            gemstat_dp_t old_sum = sum;
            #endif //DEBUG
            
            sum += compFactorInt( sites[ j ], sites[ i ] ) * Z[ j ];
            
            #ifdef DEBUG
            if( sum != sum || isinf( sum )) //DEBUG
            {
                cout << "Old sum:\t" << old_sum << endl;
                cout << "Factors:\t" << sites[ i ].factorIdx << "\t" << sites[ j ].factorIdx << endl;
                cout << "compFactorInt:\t" << compFactorInt( sites[ j ], sites[ i ] ) << endl;
                cout << "Z[j]:\t" << Z[ j ] << endl;
                cout << i << "\t" << j << "\t" << factorIntMat( (sites[i]).factorIdx, (sites[j]).factorIdx ) << endl;
                cout << "DEBUG: sum nan/inf\t"<< sum << endl;
                exit(1);
            }
            #endif //DEBUG
        }

        Z[i] = bindingWts[ i ] * sum;
        
        #ifdef DEBUG
        if( Z[i]!=Z[i] ) //DEBUG
        {
            cout << "DEBUG: Z bindingWts[i]: " << sites[i].factorIdx << "\t" << bindingWts[ sites[i].factorIdx ] <<"\t" << sum << endl;
            exit(1);
        }
        #endif //DEBUG
        Zt[i] = Z[i] + Zt[i - 1];
        //cout << "debug: Zt[i] = " << Zt[i] << endl;
    }

    // the partition function
    // 	gemstat_dp_t Z_bind = 1;
    // 	for ( int i = 0; i < sites.size(); i++ ) {
    // 		Z_bind += Z[ i ];
    // 	}
    return Zt[n];
}

double Rates_ExprFunc::compPartFuncA() const
{
    int n = sites.size() - 1;

    // initialization
    vector< gemstat_dp_t > Z( n + 1 );
    Z[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ j ], sites[ i ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * this->txpEffects[ sites[ i ].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    return Zt[n];
}

double Rates_ExprFunc::compPartFuncB() const
{
    int n = sites.size() - 1;

    // initialization
    vector< gemstat_dp_t > Z( n + 1 );
    Z[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ j ], sites[ i ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * this->repEffects[ sites[ i ].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    return Zt[n];
}

double Rates_ExprFunc::compPartFuncAB() const
{
    int n = sites.size() - 1;

    // initialization
    vector< gemstat_dp_t > Z( n + 1 );
    Z[0] = 1.0;
    vector< gemstat_dp_t > Zt( n + 1 );
    Zt[0] = 1.0;

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        gemstat_dp_t sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            sum += compFactorInt( sites[ j ], sites[ i ] ) * Z[ j ];
        }
        Z[ i ] = bindingWts[ i ] * this->txpEffects[ sites[ i ].factorIdx ] * this->repEffects[ sites[ i ].factorIdx ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    return Zt[n];
}
