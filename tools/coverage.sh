#!/usr/bin/env bash

OBJ="${PWD}/CMakeFiles/abcd.dir/src/"

# update files
rsync -ra --delete --include='*.gc*' --exclude='*.o'  ${OBJ}/* cov/gcv

GCOVFILES="${PWD}/cov/gcv"

if [ "$1" == "clean" ]
then
    make clean
    find $GCOVFILES -iname '*.gcno' -exec rm -f {} \;
    find $GCOVFILES -iname '*.gcda' -exec rm -f {} \;
    if [ "$2" != "run" ]
    then
        exit 0
    fi

    cmake .. -DTEST=ON -DDEBUG=ON
    make -j9
fi

export GCOVOUTPUT="${PWD}/cov/"
d=$(date +"%m-%d-%y_%H-%M-%S")

if [ "$1" == "main" ]
then
    d="main"
fi

#sequential tests
./Debug/runTests

#parallel tests
mpirun -np 1 ./Debug/runParallelTests
mpirun -np 2 ./Debug/runParallelTests
mpirun -np 7 ./Debug/runParallelTests
mpirun -np 16 ./Debug/runParallelTests

echo "Extracting data..."
lcov --directory $GCOVFILES --capture --output-file $GCOVOUTPUT/app.info  &> /dev/null
lcov --remove $GCOVOUTPUT/app.info "*/mmio*" -o $GCOVOUTPUT/abcd_test.info &> /dev/null

# remove this because it's meant to be in production
lcov --remove $GCOVOUTPUT/abcd_test.info "*/s_utils*" -o $GCOVOUTPUT/abcd_test.info &> /dev/null
# we do not care about libs
lcov --remove $GCOVOUTPUT/abcd_test.info "*/abcd/lib/*" -o $GCOVOUTPUT/abcd_test.info &> /dev/null

lcov --extract $GCOVOUTPUT/abcd_test.info "*/abcd/*" -o $GCOVOUTPUT/abcd_test.info &> /dev/null
echo "Generating HTML output..."
genhtml --output-directory "$GCOVOUTPUT/cov_htmp_$d" $GCOVOUTPUT/abcd_test.info

# zero the counters
# lcov --directory $GCOV_PREFIX --zerocounters
