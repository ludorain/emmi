#!/usr/bin/env bash
set -euo pipefail

RUNNING_FILE="/tmp/EMMI.running"
CMD="/emmi/measure/measure_2/start.sh"

BOARD="31"
SENSOR="A1"

case "${SENSOR}" in
    *1) MODES="iv map vover temp" ;;
    *2) MODES="iv map" ;;
    *)  echo "Unknown sensor: $sensor"; exit 1 ;;
esac

### HACK
MODES="map vover temp iv"

wait_until_finished() {
    # Wait until the running file appears
    while [[ ! -e "$RUNNING_FILE" ]]; do
	echo "Wait until the running file appears"
        sleep 1
    done

    # Wait until the running file disappears
    while [[ -e "$RUNNING_FILE" ]]; do
	echo "Wait until the running file disappears"
        sleep 1
    done
}

for mode in ${MODES}; do
    echo "Starting measurement: $mode"

    ${CMD} ${BOARD} ${SENSOR} ${mode}

    echo "Waiting for completion..."
    wait_until_finished

    echo "Measurement $mode completed"
done

echo "All measurements completed"
