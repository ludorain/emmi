#! /usr/bin/env bash
set -x
set -u

plh250add="10.0.8.10"
kydir="/home/eic/CERNIA/keithley"

vmin=48 # S13360
vmax=60 # S13360
vstp=.1
nrep=10

TEMPS="18 20 22"
SLEEP=15m

### HACK

#TEMPS="20"

vmin=48
#vmax=54

RUNNAME=$(date +%Y%m%d-%H%M%S)

ky_setup() {
    ${kydir}/serial-device-cmd.py --device ttyS0 --list ${kydir}/setup.auto.cmd
    sleep 5
}

ky_read() {
    current=$(${kydir}/serial-device-cmd.py --device ttyS0 --cmd READ?)
    echo $current
}

plh250() {
    cmd=$1
    echo "${1}" | nc -W1 -w1 ${plh250add} 9221 | tr -d '\r'
}

ivread() {
    nread=$1
    voltage=$(vbias_get)
    for I in $(seq 1 ${nread}); do
	current=$(ky_read)
	echo ${voltage} ${current}
    done
}

vbias_on() {
    while [ $(plh250 "OP1?") != 1 ]; do
	plh250 "OP1 1; *OPC?"
	sleep 1
    done
}

vbias_off() {
    while [ $(plh250 "OP1?") != 0 ]; do
	plh250 "OP1 0; *OPC?"
	sleep 1
    done
}

vbias_set() {
    voltage=$1
    plh250 "V1 ${voltage}; *OPC?"
    sleep 1
}

vbias_get() {
    voltage=""
    while [ -z "$voltage" ]; do
	voltage=$(plh250 "V1?" | awk '{print $2}')
	sleep 1
    done
    echo ${voltage}
}

ivscan() {
    tag=$1
    rm -f ${tag}.IVscan.txt
    
    vbias_off
    vbias_set ${vmin}
    vbias_on
    sleep 1m
    
    for volt in $(seq ${vmin} ${vstp} ${vmax}); do
	vbias_set ${volt}
	ivread ${nrep} | tee -a ${tag}.IVscan.txt
    done
    
    vbias_set ${vmin}
    vbias_off
}

ivscan_loop() {
    ky_setup
    for tset in ${TEMPS}; do
	/emmi/system/peltier/update-tset.sh ${tset} && sleep ${SLEEP}
	tag="${tset}c"
	ivscan ${tag}
	root -b -q -l "/home/eic/CERNIA/measure/IVscan/IVscan.C(\"${tag}.IVscan.txt\")"
    done
}

### stop read_keithley
systemctl --user stop read_keithley.service

/home/eic/bin/telegram_message.sh "Started EMMI IV scan: ${RUNNAME}"

DDIR="/home/eic/DATA/EMMI/ivscan"
mkdir -p ${DDIR}/${RUNNAME}
cd ${DDIR}/${RUNNAME}
{ time -p ivscan_loop; } &>> ivscan.log
cd -

/home/eic/bin/telegram_message.sh "EMMI IV scan completed: ${RUNNAME}"

### start read_keithley
systemctl --user start read_keithley.service

