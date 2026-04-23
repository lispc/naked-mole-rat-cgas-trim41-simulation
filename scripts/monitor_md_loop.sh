#!/bin/bash
# Continuous MD monitoring - runs every hour

while true; do
    bash scripts/monitor_md.sh
    sleep 3600
done
