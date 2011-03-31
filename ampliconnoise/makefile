all:
	(cd FCluster && make)
	(cd FClusterM && make)
	(cd FastaUnique && make)
	(cd NDist && make)
	(cd PyroDist && make)
	(cd PyroNoise && make)
	(cd SeqDist && make)
	(cd SeqDistM && make)
	(cd PyroNoiseM && make)
	(cd SeqNoise && make)
	(cd SeqNoiseM && make)
	(cd SplitClusterEven && make) 		
	(cd Perseus && make) 		

clean:
	(cd FCluster && make clean)
	(cd FClusterM && make clean)
	(cd FastaUnique && make clean)
	(cd NDist && make clean)
	(cd PyroDist && make clean)
	(cd PyroNoise && make clean)
	(cd PyroNoiseM && make clean)
	(cd SeqDist && make clean)
	(cd SeqNoise && make clean)
	(cd SeqDistM && make clean)
	(cd SeqNoiseM && make clean)
	(cd SplitClusterEven && make clean)
	(cd Perseus && make clean)
	(cd bin && rm -rf *)
install: bin
	cp FCluster/FCluster bin
	cp FClusterM/FClusterM bin
	cp FastaUnique/FastaUnique bin
	cp NDist/NDist bin
	cp PyroDist/PyroDist bin
	cp PyroNoise/PyroNoise bin
	cp PyroNoiseM/PyroNoiseM bin
	cp SeqDist/SeqDist bin
	cp SeqNoise/SeqNoise bin
	cp SeqDistM/SeqDistM bin
	cp SeqNoiseM/SeqNoiseM bin
	cp SplitClusterEven/SplitClusterEven bin
	cp Perseus/Perseus bin
bin:
	mkdir bin
