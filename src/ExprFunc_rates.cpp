#ifdef BENCHMARK
    #include <chrono>
    #include <iostream>
#endif //BENCHMARK

#include "ExprFunc_rates.h"

double Rates_ExprFunc::predictExpr( const vector< double >& factorConcs )
{
    
    #ifdef BENCHMARK
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time, end_time;
    std::chrono::duration<double> run_time;
    
    start_time = std::chrono::high_resolution_clock::now();
    #endif //BENCHMARK
    
    // compute the Boltzman weights of binding for all sites
    setupBindingWeights(factorConcs);

    /*
    double Z_O  = compPartFuncO();
    double Z_A  = compPartFuncA();
    double Z_B  = compPartFuncB();
    double Z_AB = compPartFuncAB();
    */
    
    rates_opt_return_t all_parts = compAllParts();
    gemstat_dp_t Z_O = all_parts.O;
    gemstat_dp_t Z_A = all_parts.A;
    gemstat_dp_t Z_B = all_parts.B;
    gemstat_dp_t Z_AB = all_parts.AB;
    
    

    gemstat_dp_t Z_total = Z_O + Z_A + Z_B + Z_AB;

    gemstat_dp_t prob_A = (Z_A + Z_AB) / Z_total;
    gemstat_dp_t prob_B = (Z_B + Z_AB) / Z_total;


    //TODO: Handle K_max idea

    /** I have no idea what the pis were for in this originally */
    //return (prob_A*prob_B*par.pis[seq_num])/(prob_A + prob_B*par.pis[seq_num]);
    
    #ifdef BENCHMARK
    end_time = std::chrono::high_resolution_clock::now();
    run_time = end_time - start_time;
    std::cout << "Duration: " << run_time.count() << " sec" << std::endl;
    #endif //BENCHMARK
    
    return (prob_A*prob_B)/(prob_A + prob_B);
}

/*****************
*Slow version. Only kept for debugging / confirmation.
*The version that does all at once is at least 1/3 faster.
*********************/

/*


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
*/

/*********
*Faster version of the reccurance.
* IMHO, it is also easier to read and could be made yet faster by using
* vector intrinsics.
*
**********/

typedef struct {gemstat_dp_t O, A, B, AB;} all_parts_dp_t;

rates_opt_return_t Rates_ExprFunc::compAllParts() const
{
    int n = n_sites;

    // initialization
    vector< all_parts_dp_t > Z( n + 1 );
    Z[0] = {1.0,1.0,1.0,1.0};
    vector< all_parts_dp_t > Zt( n + 1 );
    Zt[0] = {1.0,1.0,1.0,1.0};

    // recurrence
    for ( int i = 1; i <= n; i++ )
    {
        auto site_i_alpha_a = this->txpEffects[ sites[ i ].factorIdx ];
        auto site_i_alpha_r = this->repEffects[ sites[ i ].factorIdx ];
        
        all_parts_dp_t sum = Zt[boundaries[i]];
        for ( int j = boundaries[i] + 1; j < i; j++ )
        {
            if ( siteOverlap( sites[ j ], sites[ i ], motifs ) ) continue;
            gemstat_dp_t factor_interaction = compFactorInt( sites[ j ], sites[ i ] );
            
            //Could be vector intrinsics...
            sum.O  += factor_interaction * Z[ j ].O;
            sum.A  += factor_interaction * Z[ j ].A;
            sum.B  += factor_interaction * Z[ j ].B;
            sum.AB += factor_interaction * Z[ j ].AB;
        }
        
        //O
        Z[ i ].O  = bindingWts[ i ] * sum.O;
        //A
        Z[ i ].A  = bindingWts[ i ] * sum.A  * site_i_alpha_a ;
        //B
        Z[ i ].B  = bindingWts[ i ] * sum.B  *                  site_i_alpha_r ;
        //AB
        Z[ i ].AB = bindingWts[ i ] * sum.AB * site_i_alpha_a * site_i_alpha_r ;//could be precomputed
        
        
        Zt[i].O = Z[i].O + Zt[i - 1].O;
        Zt[i].A = Z[i].A + Zt[i - 1].A;
        Zt[i].B = Z[i].B + Zt[i - 1].B;
        Zt[i].AB = Z[i].AB + Zt[i - 1].AB;
    }

    return {Zt[n].O,Zt[n].A,Zt[n].B,Zt[n].AB};
}
