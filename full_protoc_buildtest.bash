#bash protoc.bash
make && python make_test_input.py && bash example_run.bash < test_serialized.out
