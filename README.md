# EarphoneBenchmarkLab v2.1
An open-source headphone testing/analyzing tool
## Note
- please run in Matlab R2018b or higher
- please follow the terms and conditions of GPLv3
- please use it at your own risk
## How to use
- HpTest_Raw.m: data analyzer/calculator
   - data will be written to the table defined in _resultFile_. please make sure labels are added to the file before running. 
   - target frequency response is defined in _curve_
   - please locate impulse response and other data in _readWav()_
- THDtest.m: non-linear distortion measuring tool
   - use system's default audio channels(windows only, further tests are needed if running on a different os)
   - data will be written to _path_
   - pre-existing files will not be replaced by default, remove those files manually if needed.
   
## Not familiar with Matlab？
I recommend [this book](https://www.amazon.com/Matlab-Practical-Introduction-Programming-Problem-ebook-dp-B00DG25ITW/dp/B00DG25ITW/ref=mt_kindle?_encoding=UTF8&me=&qid=)
