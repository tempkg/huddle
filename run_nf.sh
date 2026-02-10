#!/bin/bash

set -e
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p logs


nohup conda run -n nf-env nextflow run main.nf \
  -c nextflow.config \
  -name "run_${TIMESTAMP}" \
  -with-report "logs/report_${TIMESTAMP}.html" \
  -with-trace "logs/trace_${TIMESTAMP}.txt" \
  > "logs/nextflow_${TIMESTAMP}.log" 2>&1 &


# Save PID
echo $! > nextflow.pid
echo "Nextflow running with PID: $!"
echo "Check logs in: logs/"

