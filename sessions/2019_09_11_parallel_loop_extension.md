# Parallelizing for loops

Hopefully, you now see that for loops are super useful and can save you a lot of time! However, the loops we've been working with go through each item in the "list of things we are going to do stuff to" one at a time. One way to get a speed up is to do the stuff to all of the "list of things we are going to stuff to" all at once, i.e. making our loop run on "the things" in parallel.  

### Our wee little loop we will parallelize
To demonstrate this, I've got an example for loop below. This loop takes every 100th number from 0 to 100,000, it then counts up from that number to 99+that number, and then counts down again to the original number. It shoots all these numbers into a file called array_test.txt.
```
for i in `seq 0 100 100000`; # These are the numbers from 0 to 100,000 that we will do stuff to
  do for j in `seq 1 99`; # We are using this loop to do stuff 99 times in a row
    do (( i = i + 1 )); # Each of those 99 times, we are increasing the value of $i by 1
    echo $i >> array_test.txt; # And writing these new values of $i out to a file called array_test.txt
  done;
  for j in `seq 1 99`; # We are using this loop to do stuff 99 times in a row
    do (( i = i - 1 )); # Each of those 99 times, we are decreasing the value of $i by 1
    echo $i >> array_test.txt; # And writing these new values of $i out to a file called array_test.txt
  done;
done;  
```
To look at the speeds ups we get from parallelizing, we first need to know how long this loop takes. We are going to figure this out by running the `date` command before and after this loop to see how long it takes. Here's a version you can copy and paste that has the `date` command on each end (with comments removed, and using ';' instead of breaking the loop across lines:
```
date; for i in `seq 0 100 100000`;   do for j in `seq 1 99`;     do (( i = i + 1 ));     echo $i >> array_test.txt;   done;   for j in `seq 1 99`;     do (( i = i - 1 ));     echo $i >> array_test.txt;   done; done; date
```
Here's the results of that:
```
Thu 12 Sep 2019 11:13:42 NZST
Thu 12 Sep 2019 11:14:11 NZST
```
Our original loop took 29 s. Our file we've created, array_test.txt looks like this:
```
1
2
3
...
100002
100001
100000
```
OK, so that's the loop we'll try and speed up!

### First stop in parallelizing: creating an array
For the parallelization to work, we need to feed the loop an array. The first thing we need to do, therefore, is define our array: this is the thing we are going to do stuff to. In this case it is a list of numbers from 0 to 100000, sampling every 100 but it could be a list of files too. The key is to define the array by putting whatever group of things we want inside of the parentheses.
```
array=(`seq 0 100 100000`)
```
Let's double check that our array does in fact go from 0 to 100,000 in increments of 100. To echo out the results of an array you have to do things a little differently to a normal variable. You have to stick its name followed by [@] (this just means all values in the array e.g. if you just wanted the 10th value, you could stick a 10 in there i.e. array[10]). You then need enclose the name+[@] in curly braces, stick a $ in front of it (just like we normally have to do for variables in bash) and then enclose the whole thing in quotes...phew!
```
echo "${array[@]}"
```

Cool, OK. So the next thing we need to do is to define our "task". This is part of the for loop in the original example above, but when working with an array we need to pull it out of the for loop and turn it into a function. Let's call this function, 'thing_we_want_to_do':
```
thing_we_want_to_do() {
		for j in `seq 1 99`;
		  do (( i = i + 1 ));
		  echo $i >> array_test.txt;
		done
    for j in `seq 1 99`;
      do (( i = i - 1 ));
      echo $i >> array_test.txt;
    done;
}
```
Here's a version you can just paste into the command line:
```
thing_we_want_to_do() { for j in `seq 1 99`; do (( i = i + 1 )); echo $i >> array_test.txt; done;  for j in `seq 1 99`; do (( i = i - 1 ));  echo $i >> array_test.txt; done; }
```

Now we are going to combine our array, and our function. To match our function 'thing_we_want_to_do' above, we are going to use the $i variable to represent the various values in our array. In this first example we are going to run a 'normal', non-parallelized loop. To see how long it  takes we are going to run the 'date' command before and after. However, before we run it, we need to remove the array_test.txt from our simple for loop above!
```
rm array_test.txt
date
for i in "${array[@]}"; 
  do thing_we_want_to_do $i;
done;
date
```
Here's a version of that you can copy and paste:
```
rm array_test.txt; date; for i in "${array[@]}";    do thing_we_want_to_do $i; done; date
```
OK, so that took 17s. It was a wee bit quicker than our original loop, because we'd already saved the `seq 0 100 100000` as an array, so the computer didn't have to do this part as well as the add 1...99 and take away 1...99.
```
Thu 12 Sep 2019 11:30:10 NZST
Thu 12 Sep 2019 11:30:37 NZST
```
However, this isn't parallelizing it. But the next step will!

### Second stop in parallelizing: for reals, here we go!

