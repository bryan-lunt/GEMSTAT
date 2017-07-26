#!/bin/bash

protoc -I=src/protocol/ --cpp_out=src/protocol/ --python_out=src/protocol/ src/protocol/*.proto
