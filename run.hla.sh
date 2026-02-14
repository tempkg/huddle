#!/bin/bash

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p logs


nohup conda run -n nf-env nextflow run main.nf \
  -name "run_${TIMESTAMP}" \
  > "logs/nextflow_${TIMESTAMP}.log" 2>&1 &


# Save PID
echo $! > nextflow.pid
echo "Nextflow running with PID: $!"
echo "Check logs in: logs/"

