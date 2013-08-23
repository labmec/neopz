
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

# Select the time command arguments
if /usr/bin/time -l ls &> /dev/null; then
    TIMEARGS="-l"
elif /usr/bin/time --verbose ls &> /dev/null; then
    TIMEARGS="--verbose"
fi
echo "TIMEARGS = $TIMEARGS"

APP="@PERFTEST_APPS_DIR@/SubStruct/Perf-SubStruct"

# Main
verbose 1 "perf01 - substruct performance test: cubo1 p1 1st step."
verbose 1 "app: $APP"

# Check if application binary exists
[ -f $APP ] || fail "Application $APP is not a file"

ERRORS=0
OKS=0

echo "Start at checkpoint 1, dump checkpoint 3 and stop"
for ns in 1 2 4 8; do

  BASEOUT="cubo1.@REAL_TYPE@.txt.ckpt1.p1.nsub$ns"

  IF="cubo1.p1.nsub$ns.t.@REAL_TYPE@.txt.ckpt1"
  OF="cubo1.p1.nsub$ns.t.@REAL_TYPE@.txt.ckpt3"
  CMD="$APP -cf1 @PERFTEST_LARGE_DATA_DIR@/SubStruct/inputs/$IF -dc3 $OF -st3 -ass_rdt $BASEOUT.ass.rdt -cre_rdt $BASEOUT.cre.rdt " 

  verbose 1 "cmd: $CMD"

  rm -f "$OF"
  /usr/bin/time $TIMEARGS $CMD &> "$BASEOUT.output.txt"
  
  GOLDEN="@PERFTEST_LARGE_DATA_DIR@/SubStruct/outputs/$OF"

  # Side by side
  # DIFFOPTIONS="--suppress-common-lines -y -W 100"
  # Regular diff
  DIFFOPTIONS="--normal"

  if diff $DIFFOPTIONS "$OF" "$GOLDEN" > "$OF.diff" ; then
      echo "[OK]: $OF matches the golden output."
      OKS=$((OKS+1))
  else
      echo "[ERROR]: Checkpoint file \"$OF\" differs from golden. Differences at $OF.diff."
      ERRORS=$((ERRORS+1))
  fi

done


echo "# of ERROR: $ERRORS"
echo "# of OK   : $OKS"

