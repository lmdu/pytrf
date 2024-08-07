#!/bin/bash

#binary directory
bindir=/mnt/d/study/pytrf/benchmark/bins

#output directory
outdir=/mnt/d/study/pytrf/benchmark/results/ssr

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

	#MISA
	measure_memory_time "perl $bindir/misa.pl $gfile"
	mv $gfile.misa $outdir/$filename.misa.out

	#SciRoKoCo
	measure_memory_time "mono $bindir/SciRoKo/SciRoKoCo.exe -i $gfile -o $outdir/$filename.sciroko.td -mode misa -m 13-7-5-4-3-3"

	#PERF
	measure_memory_time "PERF -i $gfile -m 1 -M 6 -u repeat_units.txt -o $outdir/$filename.perf.out"

	#RPTRF
	#measure_memory_time "$bindir/RPTRF -s $gfile -m 6 -t 13"
	#cat $PWD/result-*.txt > $outdir/$filename.rptrf.out
	#rm $PWD/result-*.txt

	#GMATA
	measure_memory_time "perl $bindir/GMATA/gssr.pl -i $gfile -o $outdir/$filename.gmata.out -m 1 -x 6 -r 3 -s 0"

	#Phobos
	measure_memory_time "$bindir/phobos/bin/phobos_64_libstdc++6 -M exact -l 13 -u 1 -U 6 --minPerfection 100 --outputFormat 3 --reportUnit 0 $gfile $outdir/$filename.phobos.out"

	#MREPS
	#measure_memory_time "$bindir/mreps -minsize 13 -minperiod 1 -maxperiod 6 -exp 3 -noprint -fasta $gfile" "$outdir/$filename.mreps.out"

	#dot2dot
	#measure_memory_time "$bindir/dot.linux -s $gfile -c dot.config -o $outdir/$filename.out.dot -l 1 -L 6"

	#KMERSSR
	measure_memory_time "$bindir/kmer-ssr -i $gfile -o $outdir/$filename.kmer.out -p 1-6 -r 3 -n 13"

	#SSRIT
	measure_memory_time "perl $bindir/ssrit.pl $gfile" "$outdir/$filename.ssrit.out"

	#Tantan
	measure_memory_time "$bindir/tantan -f4 -b0 -j0 -w6 -n3 $gfile" "$outdir/$filename.tantan.out"

	#PYTRF
	measure_memory_time "pytrf findstr -o $outdir/$filename.pytrf.out -r 13 7 5 4 3 3 $gfile"

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
