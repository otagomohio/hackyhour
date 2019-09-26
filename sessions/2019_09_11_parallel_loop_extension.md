# Parallelizing for loops

Hopefully, [after the last hacky hour in for loops](https://github.com/otagomohio/hackyhour/blob/master/sessions/2019_09_11.md) you now see that for loops are super useful and can save you a lot of time! To illustrate this, let's start off with a simple loop:
```
for i in `seq 0 100 100000`;
do echo $i
done
```
Each time this loop runs through it takes a new number going from 0 to 100000 by every 100th number (the `seq 0 100 100000` part), and then echoes this number out (the `echo $i` part). This shows how good computers are at repeating things. However, the loop above goes through each item in the "list of things we are going to do stuff to" one at a time. One way to get a speed up is to do the stuff to all of the "list of things we are going to stuff to" all at once, i.e. making our loop run on "the things" in parallel.  

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
    done;
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
Let's try the function by running it:

```
rm array_test.txt
thing_we_want_to_do
```

You should see a new file array_test.txt:

```
1
2
3
...
96
97
98
99
...
```

Now we are going to combine our array, and our function. To match our function 'thing_we_want_to_do' above, we are going to use the $i variable to represent the various values in our array. In this first example we are going to run a 'normal', non-parallelized loop. To see how long it  takes we are going to run the 'date' command before and after. However, before we run it, we need to remove the array_test.txt from our simple for loop above or we will just keep expanding that old file!
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
OK, so that took 27s. It was a wee bit quicker than our original loop, because we'd already saved the `seq 0 100 100000` as an array, so the computer didn't have to do this part as well as the add 1...99 and take away 1...99.
```
Thu 12 Sep 2019 11:30:10 NZST
Thu 12 Sep 2019 11:30:37 NZST
```
However, this isn't parallelizing it. Our array_test.txt file has still been written in sequential order:
```
1
2
3
...
100002
100001
100000
```
This is because our loop has run one after the other using one thread. However the next step we do will parallelize this whole thing!

### Second stop in parallelizing: for reals, here we go!
Are you ready for this? To parallelize the code it is one small sneaky change: we just replace the semicolon before done with an '&'. Instead of waiting for each $i to be done, it just starts it and then starts the next one without waiting for it to be done.

```
rm array_test.txt
date
for i in "${array[@]}"; 
  do thing_we_want_to_do $i & done
date
```
Here's a version you can copy and paste:
```
rm array_test.txt; date; for i in "${array[@]}";    do thing_we_want_to_do $i & done; date
```

7s! Phew, much faster! My machine has 4 cores, so by going parallel it has allowed us about a 4x speed up
```
Thu 12 Sep 2019 11:41:23 NZST
Thu 12 Sep 2019 11:41:30 NZST
```

Apart from the speed, there is another way we can tell that we've got multiple threads going at the same time. If we look at the contents of array_test.txt this time, you'll see something like this:
```
201
202
401
...
96902
96901
96900
```
All the values are there, but they are not in order, because they correspond to different "i" array input values (e.g. 200, 400, 96000). This is because the code is running in PARALLEL (which is also why it was speedier: 7 s vs 27 s,  because many different array values were being run at the same time. This obvioulsy isn't a good idea if you want output that completes in a sequential order (but there are ways of getting around this e.g. writing out separate log files and then combining and sorting them at the end), but it can be a great idea if you've got a bunch of independent tasks you want done on similar things (e.g trimming adaptor off sequencing reads on a bunch of different files etc etc).

The one thing you need to be careful about is that this parallel loop will take up all the CPUs available to it. We can limit it using another way of achieving parallelization: xargs. 

### Playing nice and not stealing all the cores: using xargs for parallelization
To use xargs we 'll need to tweak our function slightly so that it knows it is taking its arguments from the command line (the i=$1 bit)
```
thing_we_want_to_do() {
    i=$1;
    for j in `seq 1 99`;
      do (( i = i + 1 ));
      echo $i >> array_test.txt;
    done;
    for j in `seq 1 99`;
      do (( i = i - 1 ));
      echo $i >> array_test.txt;
    done;
}
```
Version you can copy and paste:
```
thing_we_want_to_do() { i=$1; for j in `seq 1 99`; do (( i = i + 1 )); echo $i >> array_test.txt; done;  for j in `seq 1 99`; do (( i = i - 1 ));  echo $i >> array_test.txt; done; }
```
Let's try the function buy running it from 100:

```
rm array_test.txt
thing_we_want_to_do
```

We then need to export our thing_we_want_to_do function to the current shell so that xargs will be able to see it.
```
export -f thing_we_want_to_do
```
We'll also remove the version of array_test.txt we have from last time:
```
rm array_test.txt
```
OK, so xargs. Xargs would rather take our arguments one at a time than as an array, so we are going to take the original seq statement we used to make our array and feed it directly to xargs (if you were processing a bunch of files or something similar, you could feed it the `ls` command you used to identify those files). Other things going on here: 
* -n 1 means take the arguments one at a time 
* -P 1 means use just one processor (so this should be about the speed of our original one by one for loop). 
* -I starting_i means, take the value that has just been fed to xargs and call it 'starting_i'. The starting_i in the bash -c 'thing_we_want_to_do starting_i' string is then substituted with the actual number value read in from the seq 0 100 100000 statement. Then our little function then takes that number and calls it $i for its calcualtions (adding 1 through 99, taking away 1 through 99)
```
date
seq 0 100 100000 | xargs -n 1 -P 1 -I starting_i bash -c 'thing_we_want_to_do starting_i'
date
```
Here's a version you can copy and paste:
```
date; seq 0 100 100000 | xargs -n 1 -P 1 -I starting_i bash -c 'thing_we_want_to_do starting_i'; date
```

As expected, with one processor is slower than using all four cpus. This took 37 s which is also a little slower than even our original loop! The overhead for managing the threads through xargs likely makes it a little slower than just a normal for loop in this case.
```
Thu 12 Sep 2019 11:54:17 NZST
Thu 12 Sep 2019 11:54:53 NZST
```
The output (array_test.txt) is also sequential because it is just one processor ticking through. 
```
1
2
3
...
100002
100001
100000
```

Let's remove that (rm array_test.txt) and try again with 4 processors (the max on my machine):
```
rm array_test.txt
date
seq 0 100 100000 | xargs -n 1 -P 4 -I starting_i bash -c 'thing_we_want_to_do starting_i'
date
```
Here's a version you can copy and paste:
```
rm array_test.txt; date; seq 0 100 100000 | xargs -n 1 -P 4 -I starting_i bash -c 'thing_we_want_to_do starting_i'; date
```

18 seconds - speedier, but not as speedy as doing the original parallel loop using the array. Again, this is likely due to the overhead for managing threads but if you were in a situation where it would be a party foul to steal all the cpus (say, on a shared machine for example), then xargs would be your friend!
```
Thu 12 Sep 2019 12:03:41 NZST
Thu 12 Sep 2019 12:03:59 NZST
```
As for our original parallel array, the output captured in array_test.txt is not sequential:
```
101
1
102
...
100002
100001
100000
```
However, if instead of writing to a single file, we were trimming sequencing reads and writing them out to separate files, you could see how parallel processing could be really speedy and useful!
