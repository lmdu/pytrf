#!/bin/bash

#binary directory
bindir=/mnt/d/study/pytrf/benchmark/bins

#output directory
outdir=/mnt/d/study/pytrf/benchmark/results/atr

#summary output file
outfile=$outdir/summary.txt

#store benchmark time and memory
tempfile=time_mem.tmp

#record memory
memorys=()

#record elapsed time
times=()

#number of programs
num=-1

#input fasta files
gfiles=$@

measure_memory_time(){
	if [ $# -eq 1 ]; then
		/usr/bin/time -f "%e %M" -o $tempfile $1
	else
		{ /usr/bin/time -f "%e %M" -o $tempfile $1 > $2 ;}
	fi

	ret=$?

	let num++

	if [ $ret -eq 0 ]; then
		arr=($(cat $tempfile))
		times[$num]=${arr[0]}
		memorys[$num]=${arr[1]}

		#remove temp file
		if [ -e "$tempfile" ]; then
			rm $tempfile
		fi
	else
		times[$num]=0
		memorys[$num]=0
	fi
}

printf "genome\tsize\tcount\tMISA\t\tSciRoKoCo\t\tPERF\t\tGMATA\t\tPhobos\t\tKmer-SSR\t\tSSRIT\t\tTantan\t\tPytrf\t\t\n" > $outfile

for gfile in ${gfiles[@]}; do
	memorys=()
	times=()
	num=-1
	filename=$(basename $gfile)
	filename="${filename%.*}"

	#TRF
	measure_memory_time "trf $gfile 2 7 7 80 10 50 100 -h -d"
	mv $gfile.2.7.7.80.10.50.100.dat $outdir/$filename.trf.out

	#SciRoKoCo
	measure_memory_time "mono $bindir/SciRoKo/SciRoKoCo.exe -i $gfile -o $outdir/$filename.sciroko.td -mode mmfp -seedl 10"

	#Phobos
	measure_memory_time "$bindir/phobos/bin/phobos_64_libstdc++6 -M exact -l 13 -u 1 -U 6 --minPerfection 100 --outputFormat 3 --reportUnit 0 $gfile $outdir/$filename.phobos.out"

	#ULTRA
	measure_memory_time "$bindir/ultra --minunit 1 -o $outdir/$filename.ultra.out  $gfile"

	#Tantan
	measure_memory_time "$bindir/tantan -f4 -w100 -n3 $gfile" "$outdir/$filename.tantan.out"

	#PYTRF
	measure_memory_time "pytrf findatr -M 100 -o $outdir/$filename.pytrf.out $gfile"

	#get genome information
	array=($(python3 get_fasta_info.py $gfile))

	#genome size
	gsize=${array[0]}

	#sequence counts in genome
	seqcounts=${array[1]}

	#print result
	printf "%s\t%s\t%s" $filename $gsize $seqcounts >> $outfile
	for((i=0;i<=$num;i++)); do
		printf "\t%.2f\t%.2f" ${memorys[$i]} ${times[$i]} >> $outfile
	done
	printf "\n" >> $outfile
done
