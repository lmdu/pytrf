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
	/usr/bin/time -f "%e %M" -o $tempfile $1 > /dev/null 2>&1

	let num++

	if [ ! ${memorys[$num]} > 0 ]; then
		memorys[$num]=0
		times[$num]=0
	fi

	arr=($(cat $tempfile))

	#clear temp file
	if [ -e "$tempfile" ]; then
		rm "$tempfile"
	fi

	times[$num]=${arr[0]}
	memorys[$num]=${arr[1]}
}

printf "genome\tsize\tcount\tMISA\t\tPERF\t\tRPTRF\t\tGMAT\t\tMREPS\t\tKmer-SSR\t\tPytrf\t\t\n"

for gfile in ${gfiles[@]}; do
	memorys=()
	times=()
	filename=$(basename $gfile)
	filename="${filename%.*}"

	#MISA
	measure_memory_time "perl bin/misa.pl $gfile"
	mv $gfile.misa.out results/$filename.misa.out

	#SciRoKoCo
	measure_memory_time "bin/SciRoKoCo.exe -i $gfile -mode misa -m 14-7-5-4-3-3"
	mv results.sciRo results/$filename.sciroko.out

	#PERF
	measure_memory_time "PERF -i $gfile -m 1 -M 6 -u repeat_units.txt -o results/$filename.perf.out"

	#RPTRF
	measure_memory_time "bin/RPTRF -s $gfile -m 6 -t 14"
	cat result-*.txt > results/$filename.rptrf.out
	rm result-*.txt

	#GMAT
	measure_memory_time "perl bin/gssr.pl -i $gfile -o results/$filename.gmat.out -m 1 -x 6 -r 3 -s 0"

	#MREPS
	#measure_memory_time "bin/mreps -fasta -minperiod 1 -maxperiod 6 -minsize 14 -exp 3 $gfile > results/$filename.mreps.out"

	#KMERSSR
	measure_memory_time "bin/kmer-ssr -i $gfile -o results/$filename.kmer.out -p 1-6 -r 3 -n 14"

	#PYTRF
	measure_memory_time "pytrf findstr -o results/$filename.pytrf.out -r 14 7 5 4 3 3 $gfile"

	#get genome information
	array=($(python3 get_fasta_info.py $gfile))

	#genome size
	gsize=${array[0]}

	#sequence counts in genome
	seqcounts=${array[1]}

	#print result
	printf "%s\t%s\t%s" $filename $gsize $seqcounts
	for((i=0;i<=$num;i++)); do
		printf "\t%.2f\t%.2f" ${memorys[$i]} ${times[$i]}
	done
	printf "\n"
done
