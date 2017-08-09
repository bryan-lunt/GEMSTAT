#ifndef GEMSTAT_IO_H
#define GEMSTAT_IO_H

#include <string>
#include <map>

#include "ExprPredictor.h"

#include "Tools.h"

int readEdgelistGraph( const string& filename, const map<string, int>& factorIdxMap, IntMatrix& destination, bool directed);

int readFactorThresholdFile( const string& filename, vector< double >& destination, int nFactors);

int readFactorRoleFile(const string& filename, const map<string, int>& factorIdxMap,  vector< bool>& actIndicators, vector<bool>& repIndicators);
int readFactorRoleFile(const string& filename, const map<string, int>& factorIdxMap,  vector< bool>& actIndicators, vector<bool>& repIndicators, vector<bool>& pioneer_indicators);

//TODO: Shouldn't this function check that it is reading the correct number of sites?
int readAxisWeights(const string& filename, vector< int >& axis_start, vector< int >& axis_end, vector< double >& axis_wts);

int writePredictions(const string& filename, ExprPredictor& predictor, const Matrix& exprData, vector< string >& expr_condNames, bool write_gt, bool fix_beta = false);

#endif
