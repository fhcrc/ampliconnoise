all:
	(cd FCluster && make)
	(cd FastaUnique && make)
	(cd NDist && make)
	(cd PyroDist && make)
	(cd PyroNoise && make)
	(cd SeqDist && make)
	(cd PyroNoiseM && make)
	(cd PyroNoiseA && make)
	(cd SeqNoise && make)
	(cd SplitClusterEven && make)
	(cd SplitClusterClust && make)  		
	(cd Perseus && make) 		
	(cd PerseusD && make) 		
clean:
	(cd FCluster && make clean)
	(cd FastaUnique && make clean)
	(cd NDist && make clean)
	(cd PyroDist && make clean)
	(cd PyroNoise && make clean)
	(cd PyroNoiseM && make clean)
	(cd PyroNoiseA && make clean)
	(cd SeqDist && make clean)
	(cd SeqNoise && make clean)
	(cd SplitClusterEven && make clean)
	(cd SplitClusterClust && make clean) 
	(cd Perseus && make clean)	
	(cd PerseusD && make clean)
	(cd bin && rm -rf *)
install:
	cp FCluster/FCluster bin
	cp FastaUnique/FastaUnique bin
	cp NDist/NDist bin
	cp PyroDist/PyroDist bin
	cp PyroNoise/PyroNoise bin
	cp PyroNoiseM/PyroNoiseM bin
	cp PyroNoiseA/PyroNoiseA bin
	cp SeqDist/SeqDist bin
	cp SeqNoise/SeqNoise bin
	cp SplitClusterEven/SplitClusterEven bin
	cp SplitClusterClust/SplitClusterClust bin
	cp Perseus/Perseus bin
	cp PerseusD/PerseusD bin
