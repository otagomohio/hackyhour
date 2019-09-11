# Parallelizing for loops

Hopefully, you now see that for loops are super useful and can save you a lot of time! However, the loops we've been working with go through each item in the "list of things we are going to do stuff to" one at a time. One way to get a speed up is to do the stuff to all of the "list of things we are going to stuff to" all at once, i.e. making our loop run on "the things" in parallel.  

### Our wee little loop we will parallelize
To demonstrate this, I've got an example for loop below. This loop takes every 100th number from 0 to 100,000, then adds from 1 to 99 to each of the original numbers, and shoots all this into file called array_test.txt.
```
for i in `seq 0 100 100000`; # These are the numbers from 0 to 100,000 that we will do stuff to
  do for j in `seq 1 99`; # We are using this loop to do stuff 99 times in a row
    do (( i = i + 1 )); # Each of those 99 times, we are increasing the value of $i by 1
    echo $i >> array_test.txt; # And writing these new values of $i out to a file called array_test.txt
  done;
done;  
```
