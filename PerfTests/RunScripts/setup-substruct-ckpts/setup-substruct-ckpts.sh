
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

# Main
verbose 1 "Generate the substruct checkpoints."
verbose 1 "app: $APP"
verbose 1 "cmd: $CMD"

# Check if application binary exists
[ -f $APP ] || fail "Application $APP is not a file"

DATADIR="@PERFTEST_DATA_DIR@"

CMD="$APP -mc $DATADIR/Substruct/inputs/cube1.txt" 

#for plevel in 1 2; do
for plevel in 1; do
  for ns in 1 2 4 8 16 32 64; do 

    echo "Generating cube chekpoints for plevel = $plevel and nsubs = $ns"

    CMD="$APP -mc \"$DATADIR/Substruct/inputs/cube1.txt\" \
        -dc1 \"$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt1\" \
        -dc2 \"$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt2\" \
        -dc3 \"$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt3\" "

    echo "CMD=$CMD"

    $APP -mc "$DATADIR/Substruct/inputs/cube1.txt" \
        -dc1 "$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt1" \
        -dc2 "$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt2" \
        -dc3 "$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt3"

    cp -f "$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt1" \
        "$DATADIR/Substruct/inputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt1"

    cp -f "$DATADIR/Substruct/outputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt2" \
        "$DATADIR/Substruct/inputs/cubo1.p$plevel.nsub$ns.t.@REAL_TYPE@.ckpt2"

  done   
done
