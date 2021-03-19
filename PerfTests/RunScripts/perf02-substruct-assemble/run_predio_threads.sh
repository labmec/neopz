# Select the time command arguments
if /usr/bin/time -l ls &> /dev/null; then
    TIMEARGS="-l"
elif /usr/bin/time --verbose ls &> /dev/null; then
    TIMEARGS="--verbose"
fi
echo "TIMEARGS = $TIMEARGS"

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

function run_cfg
{
    #   -nt_a : number of threads to assemble each submesh
    local NTA=$1
    #   -nt_d : number of threads to decompose each submesh
    local NTD=$2
    #   -nt_m : number of threads to multiply
    local NTM=$3
    #   -nt_sm : number of threads to process the submeshes
    local NTSM=$4
    #   -p : plevel
    local P=$5

    for ns in 2 4 8 16 32 64; do

        BASEOUT="predio.@REAL_TYPE@.bin.ckpt1.p$P.nsub$ns.nt_a.$NTA.nt_d.$NTD.nt_m.$NTM.nt_sm.$NTSM"

        IF="predio.p$P.nsub$ns.t.@REAL_TYPE@.bin.ckpt1"
        MD5F="@PERFTEST_LARGE_DATA_DIR@/SubStruct/outputs/$OF.md5"
        CMD="$APP -bc -cf1 @PERFTEST_LARGE_DATA_DIR@/SubStruct/inputs/$IF -chk_c3_md5 $MD5F -st3 -ass_rdt $BASEOUT.ass.rdt -cre_rdt $BASEOUT.cre.rdt " 
        CMD="$CMD -nt_a $NTA -nt_d $NTD -nt_m $NTM -nt_sm $NTSM -p $P"
        verbose 1 "cmd: $CMD"
        
        rm -f "$OF"
        /usr/bin/time $TIMEARGS $CMD &> "$BASEOUT.output.txt"
        RET=$?

        if [ $RET != 0 ]; then
            echo "[ERROR]: $APP tool returned $RET."
            ERRORS=$((ERRORS+1))
        else
            echo "[OK]: MD5 for checkpoint 3 matches golden MD5."
            OKS=$((OKS+1))
        fi
    done
}

PLEVEL=2
VERBOSE_LEVEL=1
APP="@PERFTEST_APPS_DIR@/SubStruct/Perf-SubStruct"

# Main
verbose 1 "perf01 - substruct performance test: predio step."
verbose 1 "app: $APP"

# Check if application binary exists
[ -f $APP ] || fail "Application $APP is not a file"

ERRORS=0
OKS=0

#for PLEVEL in 1 2; do
for PLEVEL in 1; do
    echo "Start at checkpoint 1, dump checkpoint 3 and stop: plevel = $PLEVEL"

    # 1 thread for all the submeshs
    for t in 0 1 2 4 8 16 32 64; do
        run_cfg $t $t $t 0 $PLEVEL
    done
    
    # 2 threads for all the submeshs
    for t in 0 1 2 4 8 16 32 64; do
        run_cfg $t $t $t 2 $PLEVEL
    done

    # 4 threads for all the submeshs
    for t in 0 1 2 4 8 16 32 64; do
        run_cfg $t $t $t 4 $PLEVEL
    done

    # 8 threads for all the submeshs
    for t in 0 1 2 4 8 16 32 64; do
        run_cfg $t $t $t 8 $PLEVEL
    done

    # 16 threads for all the submeshs
    for t in 0 1 2 4 8 16 32 64; do
        run_cfg $t $t $t 8 $PLEVEL
    done

    # 32 threads for all the submeshs
    for t in 0 1 2 4 8 16 32 64; do
        run_cfg $t $t $t 8 $PLEVEL
    done

    # 64 threads for all the submeshs
    for t in 0 1 2 4 8 16 32 64; do
        run_cfg $t $t $t 8 $PLEVEL
    done

done

echo "# of ERROR: $ERRORS"
echo "# of OK   : $OKS"

