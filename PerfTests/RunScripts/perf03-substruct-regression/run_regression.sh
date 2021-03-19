#!/bin/bash

FORCE_COUT=0
FORCE_MAKE=0
FORCE_CKPTS=0
VERBOSE_LEVEL=2

DIR_BASE=/local/julia/performance_regression
APP_CKPTS=$DIR_BASE/work/neopz_trunk_build/PerfTests/RunScripts/setup-substruct-ckpts/setup-substruct-ckpts.sh 
APP_PTEST1=$DIR_BASE/work/neopz_trunk_build/PerfTests/RunScripts/perf02-substruct-assemble/run_cubo_threads.sh
APP_PTEST2=$DIR_BASE/work/neopz_trunk_build/PerfTests/RunScripts/perf02-substruct-assemble/pthread_vs_tbb-cubo_binckpt.sh
MAKE_ARGS="-j64"

function fail {
	echo "FAIL: $@" 
	exit 1
}

function verbose {
	if [ $1 -le $VERBOSE_LEVEL ]; then
		echo $2 
	fi
}

# Main
REVp=-1
CHANGE=1

# checking whether directory with source code exists  
if [ -d $DIR_BASE/work/neopz_trunk ]; then 
	# saving number of previous revision
	REVp=`svn info $DIR_BASE/work/neopz_trunk | grep Revision | cut -d ' ' -f2`  
fi

# getting number of current revision
REV=`svn info https://neopz.googlecode.com/svn/trunk/ | grep Revision | cut -d ' ' -f2`  

# checking whether revision changed since last checkout
verbose 2 "Checking whether there has been changes since last checkout."
if [ $REV -eq $REVp ]; then 
	verbose 1 "No changes since last checkout."
	CHANGE=0
else
	verbose 2 "There has been changes since last checkout."
fi

if [ $CHANGE -eq 1 -o $FORCE_COUT -eq 1 ] 
then
	# removing source if exists
	if [ -d $DIR_BASE/work/neopz_trunk ]; then rm -r $DIR_BASE/work/neopz_trunk || fail "could not remove directory."; fi

	# checking out current revision
	verbose 2 "Checking out revision $REV from https://neopz.googlecode.com/svn/trunk/"
	svn checkout https://neopz.googlecode.com/svn/trunk/ $DIR_BASE/work/neopz_trunk || fail "could not checkout a working copy."
	verbose 1 "Check out done."
fi

if [ $CHANGE -eq 1 -o $FORCE_MAKE -eq 1 ] 
then 
	# removing builded files if existing
	if [ -d $DIR_BASE/work/neopz_trunk_build ]; then rm -r $DIR_BASE/work/neopz_trunk_build || fail "could not remove directory."; fi
	
	# creating directory for builded files of svn
	mkdir $DIR_BASE/work/neopz_trunk_build && cd $DIR_BASE/work/neopz_trunk_build || fail "could not create directory."

	# configuring cmake
	verbose 2 "Configuring CMake..."
	cmake -DREAL_TYPE=double -DCMAKE_BUILD_TYPE=Release -DBUILD_PERF_TESTS=ON -DBUILD_PROJECTS=OFF -DBUILD_TUTORIAL=OFF -DBUILD_UNITTESTING=OFF \
		-DUSING_METIS=ON -DUSING_LOG4CXX=OFF -DUSING_BOOST=ON -DUSING_OPENSSL=ON -DUSING_TBB=ON \
			-DPERFTEST_LARGE_DATA_DIR:PATH=$DIR_BASE/checkpoints -DPERFTEST_APPS_DIR:PATH=$DIR_BASE/work/neopz_trunk_build/PerfTests \
				$DIR_BASE/work/neopz_trunk || fail "problems while configuring cmake."

	# building current revision
	verbose 2 "Building source files of revision $REV."
	make $MAKEARGS || fail "could not build files of $DIR_BASE/work/neopz_trunk."
	verbose 1 "Source builded."
else
	exit 2
fi

# creating subdirectory for current revision
mkdir -p $DIR_BASE/results/r$REV || fail "could not create directory for saving results."


cp -f $DIR_BASE/work/neopz_trunk/PerfTests/neopz_perf_data/SubStruct/inputs/cube1.txt $DIR_BASE/checkpoints/SubStruct/inputs || \
	fail "could not reach one or more input files for generating checkpoints."
cp -f $DIR_BASE/work/neopz_trunk/PerfTests/neopz_perf_data/SubStruct/inputs/8andares02.txt $DIR_BASE/checkpoints/SubStruct/inputs || \
	fail "could not reach input file for generating checkpoints."

if [ $FORCE_CKPTS -eq 1 ]; then
	[ -f $APP_CKPTS ] || fail "cannot generate checkpoints ( $APP_CKPTS is not a file)."
	$APP_CKPTS || fail "could not generate checkpoints." 
	verbose 1 "Check points generated."
fi

# running performance test
[ -f $APP_PTEST1 ] || fail "cannot run performance test ($APP_PTEST1 is not a file)."
cd $DIR_BASE/work/neopz_trunk_build/PerfTests/RunScripts/perf02-substruct-assemble && rm -f *.rdt || fail "directory does not exist."
verbose 1 "Starting performance test..."
for i in 1 2 3
do
	./run_cubo_threads.sh || fail "problems while running $APP_PTEST1" 
done
verbose 1 "Performance test went succesful."
mv -f *.rdt "$DIR_BASE/results/r$REV"

[ -f $APP_PTEST2 ] || fail "cannot run performance test ($APP_PTEST2 is not a file)."
cd $DIR_BASE/work/neopz_trunk_build/PerfTests/RunScripts/perf02-substruct-assemble && rm -f *.rdt || fail "directory does not exist."
verbose 1 "Starting performance test..."
for i in 1 2 3
do
        ./pthread_vs_tbb-cubo_binckpt.sh || fail "problems while running $APP_PTEST2"
done
verbose 1 "Performance test went succesful."
mv -f *.rdt "$DIR_BASE/results/r$REV"

