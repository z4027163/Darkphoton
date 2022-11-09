find 2018A/ -type f -size -500k >& tem2018A.log
find 2018B/ -type f -size -500k >& tem2018B.log
find 2018C/ -type f -size -500k >& tem2018C.log
find 2018D1/ -type f -size -500k >& tem2018D1.log
find 2018D2/ -type f -size -500k >& tem2018D2.log
python read2018.py >& list2018.txt
find 2017C/ -type f -size -500k >& tem2017C.log
find 2017D/ -type f -size -500k >& tem2017D.log
find 2017E/ -type f -size -500k >& tem2017E.log
find 2017F/ -type f -size -500k >& tem2017F.log
python read2017.py >& list2017.txt
find 2018*/ -type f -size -500k >& 2018step2.log
find 2017*/ -type f -size -500k >& 2017step2.log
python read2018_2.py >& list2018_2.txt
python read2017_2.py >& list2017_2.txt
