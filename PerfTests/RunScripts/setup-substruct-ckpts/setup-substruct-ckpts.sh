
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

APP="@PERFTEST_APPS_DIR@/SubStruct/Perf-SubStruct"

# Main
verbose 1 "Generate the substruct checkpoints."
verbose 1 "app: $APP"
verbose 1 "cmd: $CMD"

# Check if application binary exists
[ -f $APP ] || fail "Application $APP is not a file"

DATADIR="@PERFTEST_LARGE_DATA_DIR@"

function gen_ckpts
{
    INPUT=$1
    BASENAME=$2
    READOPT=$3
    CKPTYPE=$4
    EXTRAARGS=$5

    for plevel in 1 2; do
	for ns in 2 4 8 16 32 64; do # Does not work for ns == 1

	    echo "Generating $BASENAME chekpoints for plevel = $plevel and nsubs = $ns"

	    CMD="$APP $EXTRAARGS $READOPT \"$DATADIR/SubStruct/inputs/$INPUT\" \
              -gen_c1_md5 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt1.md5\" \
        	  -gen_c2_md5 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt2.md5\" \
        	  -gen_c3_md5 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt3.md5\" \
              -nsub $ns -p $plevel"

        if [ "$INPUT" == "8andares02.txt" ]; then
            CMD="$CMD  \
              -dc1 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt1\" \
        	  -dc2 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt2\" "
            # Skip ckpt3 for predio -- Too big
        else
            CMD="$CMD  \
              -dc1 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt1\" \
        	  -dc2 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt2\" \
        	  -dc3 \"$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt3\" "
        fi

	    eval $CMD &> /dev/null
	    RET=$?

	    if [ $RET != 0 ]; then
		    echo "ERROR when executing: $CMD"
		    echo "Return code = $RET"
	        echo "CMD=$CMD"
	    else
	        cp -f "$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt1" \
		        "$DATADIR/SubStruct/inputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt1"
	        
	        cp -f "$DATADIR/SubStruct/outputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt2" \
		        "$DATADIR/SubStruct/inputs/$BASENAME.p$plevel.nsub$ns.t.double.$CKPTYPE.ckpt2"
	    fi
	done   
    done
}

gen_ckpts cube1.txt cubo1 -mc txt 
gen_ckpts cube1.txt cubo1 -mc bin -bc
gen_ckpts 8andares02.txt predio -mp bin -bc
