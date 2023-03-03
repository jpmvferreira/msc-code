for DATANAME in "LISA" "ET"; do
	if [[ "$DATANAME" == "LISA" ]]; then
		Nlist=$(seq 5 5 30)
		IDlist=$(seq 1 15)
		SAMPLES=7500
		WARMUP=2500
	elif [[ "$DATANAME" == "ET" ]]; then
		Nlist="100 250 500 750 1000"
		IDlist=$(seq 1 5)
		SAMPLES=2500
		WARMUP=1000
	fi
	
	for MODEL in "LCDM" "fotis-noOmegar"; do		
		for N in $Nlist; do
			for ID in $IDlist; do
				DATASET=${DATANAME}-N${N}-${ID}
				OUTPUT=${MODEL}_${DATASET}
				
				echo "Running model $MODEL with dataset $DATASET"
				echo "Start time: $(date +'%H:%M:%S, %d/%m/%y (%A)')"
				
				nice -n 10                         \
				smc-stan -m model/$MODEL.stan      \
				         -d data/Nvar/$DATASET.csv \
				         -o output/Nvar/$OUTPUT    \
				         -cfg config/$MODEL.yml    \
				         -s $SAMPLES               \
				         -w $WARMUP                \
				         --PSIS-LOO-CV             \
				         --WAIC                    \
				         -hp                       \
				&> STDOUTandSTDERR_$OUTPUT.log

		        echo "End time: $(date +'%H:%M:%S, %d/%m/%y (%A)')"
				echo ""
				 
				mv STDOUTandSTDERR_$OUTPUT.log          \
				output/Nvar/$OUTPUT/log/STDOUTandSTDERR.log  
			done
		done
	done
done
