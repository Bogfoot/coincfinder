#!/usr/bin/env bash
set -euo pipefail

# User-tunable defaults
coinc=250          # coincidence half-window (ps)
delay_start=9      # ns
delay_end=12       # ns
delay_step=0.01    # ns
start=0            # start time (s); 0 means beginning
stop=600           # stop time (s); 0 means end

for f in 8hMeasurement/*.bin; do
  [ -e "$f" ] || { echo "No .bin files found in 8hMeasurement" >&2; exit 1; }
  stem=$(basename "$f" .bin)
  outdir="CoincEvents/$stem"
  mkdir -p "$outdir"

  echo "Running CoincPairs on $f"
  ./CoincPairs "$f" "$coinc" "$delay_start" "$delay_end" "$delay_step" "$start" "$stop" "$outdir/rate.csv" --dump-events
  echo "Done: $f"
done
