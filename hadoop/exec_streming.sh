hadoop jar /opt/cloudera/parcels/CDH-5.6.0-1.cdh5.6.0.p0.45/lib/hadoop-mapreduce/hadoop-streaming.jar \
        -file mapper.py    -mapper mapper.py \
        -file reducer.py   -reducer reducer.py \
        -input airports/input -output airports/output \
        -numReduceTasks 1 
