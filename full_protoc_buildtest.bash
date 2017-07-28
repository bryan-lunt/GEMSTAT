#bash protoc.bash
#make && python make_test_input.py && bash example_run.bash < test_serialized.out

rm baseline.log streaming.log baseline.out
make
#comparison run
./src/seq2expr -s data_streamtest/seqs.fa -e data_streamtest/expr.tab -m data/factors.wtmx -f data/factor_expr.tab -c data/coop.txt -i data/factor_info.txt -o Direct -oo SSE -na 0 -fo baseline.out -p data_streamtest/example.par -onebeta > baseline.log

python make_test_input2.py | ./src/seq2expr_streaming > streaming.log
