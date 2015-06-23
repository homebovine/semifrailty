#!/bin/bash
sleep 2

for sta in $(seq 1 20 1000)
do
    echo $sta
    sed "s/sta <- 1/sta <- $sta/" <~/semifrailty/trunk/semifrailty/code/estimation.r>~/semifrailty/trunk/semifrailty/code/estimation$sta.r
    sed "s/estimation/estimation$sta/" <~/cluster/semifrailty.pbs> ~/cluster/semifrailty$sta.pbs
    qsub  ~/cluster/semifrailty$sta.pbs
    sleep 1
done

