#!/bin/bash

protoc -I=src/streaming/ --cpp_out=src/streaming/ --python_out=src/streaming/ src/streaming/gemstat.proto
