#!/bin/bash

dir="results/KINGCANYON/"
nsamples=$(ls -1 results/KINGSCANYON/NA_results*.txt | wc -l)

echo $nsamples

for i in $(ls results/KINGSCANYON/NA_results*.txt);
do
  ./bin/plot_inversion2.py -f $i -o file
done


for i in $(seq 1 $nsamples);
do
  mv results/KINGSCANYON/NA_results_sample$i.png results/KINGSCANYON/Inversion$i.png
done


for i in $(ls results/KINGSCANYON/NA_results*.txt);
do
  ./bin/plot_sampleDEM.py -f $i -o file
done


for i in $(seq 1 $nsamples);
do
  mv results/KINGSCANYON/NA_results_sample$i.png results/KINGSCANYON/Map$i.png
done

cd results/KINGSCANYON

for i in $(seq 1 $nsamples);
do
  montage Inversion$i.png Map$i.png -geometry +2+2 Sample$i.png
done

montage $(ls Sample*.png | sort -V -f) -geometry +2+2 AllInversions.png

rm Inversion*.png Map*.png Sample*.png
echo "Plots generated"
