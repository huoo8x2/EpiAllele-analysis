#!/usr/bin/env bash
conda activate offinder

mkdir -p /data02/hukaijie/EpiAllele/result/off-target

cas-offinder "/data02/hukaijie/EpiAllele/script/offinder/mouse_input.txt" G0 "/data02/hukaijie/EpiAllele/result/off-target/mouse_offtarget.txt"
cas-offinder "/data02/hukaijie/EpiAllele/script/offinder/human_input.txt" G0 "/data02/hukaijie/EpiAllele/result/off-target/human_offtarget.txt"
