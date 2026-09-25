#ifdef BENCHMARK
    #include <chrono>
    #include <iostream>
#endif //BENCHMARK

#include "ExprFunc_states.h"

double States_ExprFunc::predictExpr( const vector< double >& factorConcs )
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

		states_opt_return_t all_parts = compAllParts();
    gemstat_dp_t Z_O = all_parts.O;
    gemstat_dp_t Z_A = all_parts.A;
    gemstat_dp_t Z_AB = all_parts.AB;

    //Apply qbtm, there is one for both arcs. recycle pi as the one for the other arc.
    GEMSTAT_PROMOTER_DATA_T my_promoter = par.getPromoterData( this->seq_number );//should happen at object construction
    Z_A *= my_promoter.basal_trans;
    Z_AB *= my_promoter.basal_trans*my_promoter.pi;

    gemstat_dp_t Z_total = Z_O + Z_A + Z_AB;

		gemstat_dp_t prob_AB = Z_AB / Z_total;


    #ifdef BENCHMARK
    end_time = std::chrono::high_resolution_clock::now();
    run_time = end_time - start_time;
    std::cout << "Duration: " << run_time.count() << " sec" << std::endl;
    #endif //BENCHMARK

    //TODO: Handle K_max idea

    /** I have no idea what the pis were for in this originally */
    //return (prob_A*prob_B*par.pis[seq_num])/(prob_A + prob_B*par.pis[seq_num]);

		return prob_AB;


}


typedef struct {gemstat_dp_t O, A, B, AB;} all_parts_dp_t;//keeping the B for memory alignment or something.

states_opt_return_t States_ExprFunc::compAllParts() const
{
    int n = n_sites;

    // initialization
    vector< all_parts_dp_t > Z( n + 1 );
    Z[0] = {1.0,1.0,1.0,1.0};
    vector< all_parts_dp_t > Zt( n + 1 );
		Zt[0] = {1.0,1.0,0.0,1.0};

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
						//sum.B  += factor_interaction * Z[ j ].B;
            sum.AB += factor_interaction * Z[ j ].AB;
        }

        //O
        Z[ i ].O  = bindingWts[ i ] * sum.O;
        //A
        Z[ i ].A  = bindingWts[ i ] * sum.A  * site_i_alpha_a ;
        //B
				//Z[ i ].B  = bindingWts[ i ] * sum.B  *                  site_i_alpha_r ;
        //AB
        Z[ i ].AB = bindingWts[ i ] * sum.AB * site_i_alpha_a * site_i_alpha_r ;//could be precomputed


        Zt[i].O = Z[i].O + Zt[i - 1].O;
        Zt[i].A = Z[i].A + Zt[i - 1].A;
				//Zt[i].B = Z[i].B + Zt[i - 1].B;
        Zt[i].AB = Z[i].AB + Zt[i - 1].AB;
    }

		return {Zt[n].O,Zt[n].A,Zt[n].AB};
}
