
VERBOSE_LEVEL=1

function fail
{
    echo "FAIL: $@"
    exit 1
}

function verbose
{
    local level=$1
    shift;
    if [ $VERBOSE_LEVEL -ge $level ]; then
	echo $@
    fi
}

APP="@PERFTEST_APPS_DIR@/Substruct/Perf-SubStruct"
CMD="$APP -mp @PERFTEST_DATA_DIR@/Substruct/input/8andares02.txt" 

# Main
verbose 1 "perf01 - substruct performance test."
verbose 1 "app: $APP"
verbose 1 "cmd: $CMD"

# Check if application binary exists
[ -f $APP ] || fail "Application $APP is not a file"
 
/usr/bin/time -l $CMD


