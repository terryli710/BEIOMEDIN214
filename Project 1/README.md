# Project 1

 

## Important Information

1. Project 1 is due by **5:00:00****pm PDT on October 2, 2020**. Projects turned in within 24 hours after the deadline will receive a 10% deduction. Projects turned in within 48 hours after the deadline will receive a 20% deduction. Projects turned in within 72 hours after the deadline will receive a 30% deduction. **No submissions will be accepted more than 72 hours after the deadline.** 

   Note: we have given a substantial amount of time to complete this particular project due to its challenging nature. We recommend starting this assignment early, as troubleshooting can be quite time consuming.

2. Programming projects must be completed individually. You may discuss algorithms with others, but the coding should be done alone. You must explicitly name everyone with whom you discussed this project in the header comments of your code and in the collaboration attestation in the project quiz. Students must abide by the terms of the Stanford Honor Code.

3. **Remember to consult Piazza,** as many common questions will be asked and answered there. Remember that Piazza participation (asking and answering questions) makes up the attendance portion of the grade.

   

4. Prior to beginning the assignment please read the [Code Policy](https://canvas.stanford.edu/courses/122868/pages/code-policy). Because the skeleton code includes a lot of static typing to ensure implementational consistency and easier debugging, **Python 3.5 or higher is required for this course**.

   

5. Source code must run **exactly** as specified below. We will be running your code on test alignments on Gradescope.

# Introduction 

We often want to compare the DNA or protein sequences from different species to examine evolutionary relationships or to look at shared functionality. In order to make these comparisons, we need to examine which portions of these sequences align. During lecture, we covered the Smith-Waterman and Needleman-Wunsch algorithms for sequence alignment using dynamic programming. We will be implementing this algorithm using the affine gap penalty, as described in class and in [Durbin et al](https://canvas.stanford.edu/courses/122868/files/6438133/download?wrap=1)[![Preview the document](https://canvas.stanford.edu/images/preview.png)](https://canvas.stanford.edu/courses/122868/files/6438133/download?wrap=1). (equations 2.16). Your program should be implemented such that it can be used for either global ends-free or local alignment (the input file has an indicator on which alignment it wants).

Once you have implemented this algorithm, we will use sequence alignment to help identify the location of a disease-relevant RNA, look at sequence homology, and help make decisions about using a therapeutic source for a protein.

**Learning Goals**

- Understand and implement sequence alignment algorithms (N-W and S-W) with the affine gap penalty

- Implement, understand, and debug a dynamic programming based algorithm.

- Application/Interpretation:

- - Know when to use sequence alignment, understand reasons for different inputs (nucleic acid vs. protein), types of alignment (local vs. global)
  - Understand how gap parameters affect sequence alignments

 

# Implementation Notes

Your code should follow exactly the equations on page 29 of the [Durbin et al](https://canvas.stanford.edu/courses/122868/files/6438133/download?wrap=1)[![Preview the document](https://canvas.stanford.edu/images/preview.png)](https://canvas.stanford.edu/courses/122868/files/6438133/download?wrap=1) textbook, i.e., there should be 3 scoring matrices, M, Ix, and Iy. A few things to note:

 

1. According to Durbin, sequence x (or sequence 1) is on the top of the scoring matrix, while sequence y (or sequence 2) is on the left side of the scoring matrix. Therefore, the i-th character in x corresponds to the i-th *column* in the scoring matrix, while the j-th character in y corresponds to the j-th *row* in the scoring matrix. In the equations on page 29 (equations 2.16), the first index is the *column* index, and the second index is the *row* index. ***\*This is opposite to what we usually see in programming languages and to what Russ did in class.\**** **We recommend that you follow the convention that Russ described in class!** While it will work either way, the standard convention will make it easier for you to debug, and easier for you to receive partial credit.***\*
   \****

2. The Durbin et al textbook assumes the same gap penalty for both sequences. However, we have different values for the two sequences. Therefore, in the Ix equation, **d should be dy, and e should be ey, because Ix means inserting a gap in sequence y. Similarly, in the Iy equation, d should be dx, and e should be ex.**

3. Be sure that you have a strong conceptual understanding of the algorithm before you attempt to encode it. It is wise to begin with pseudocode for this project.

4. Your program should be able to trace **all** the alignments that give the best alignment score. When you need to report more than one equivalent alignment, you should leave a blank line and then print the alignments in a similar manner, as in the example. The ordering of the alignments will not affect the autograder.

5. You should be careful with rounding errors. If you are getting rounding errors when using the == operator to compare numbers, note that in general, two floats will only satisfy == if they were produced from the same operations. Instead of **x == y** use **abs((x - y)) < delta** (we provide the **fuzzy_equals** function in the starter code). If you do not use this, you may wind up missing some alignments.

6. Be aware that we may test on an asymmetric score matrix.

7. The policy is to have 

   no end-gap penalties

    imposed, on either end. 

   Thus the ends of the sequences that do not participate in the alignment do not count towards the score. The example output above has truncated the end gaps in the output. E.g. for the above example, the first alignment shown is

   ```
   ATGC
   AGGC
   ```

    

   with the leading A in sequence 1 excluded because it is aligned with an end gap. This is still a global alignment! An alternative representation of the first alignment is

   ```
   AATGC 
   _AGGC
   ```

    

   which is equivalent, and still has a score of 3.0. Your program should produce **only the first type of output**, with the end gaps omitted, and end gaps are not penalized. This requires you to change a little bit in the global alignment algorithm. Think about the initialization of the first row and first column of the score matrix.

    

Hints for local alignment vs. global alignment:

1. Local never records a score in the score matrix below 0.

   

2. When looking for best match, local looks for the best score anywhere in the matrix, while global looks only in last row and last column.

   

3. Local requires negative scores in the match matrix (note: the "match matrix" is the matrix which specifies what the score is for matching one amino acid or nucleotide with another, this is NOT the same as the "score matrix" which stores the scores for best alignments at any position; the "score matrix" is the one you are building up from top left to bottom right). Without negative scores for mismatches an alignment's score would never decrease and the local alignment would never terminate.

   

4. Remember that with ends-free alignment, you exclude all leading and trailing insertions or deletions. For example, lets say the alignment for AAAACCGTAC and CCGTAC would look as follows:

 

```
AAAACCGTAC
____CCGTAC
```

 

The leading gaps would ***NOT\*** be included in the score, and the alignment should be reported as follows:

```
CCGTAC
CCGTAC
```

  \5. During local alignment, if you trace back to a cell that contains pointers to a zero in the M matrix and a pointer to a zero in the Ix or Iy matrix, you should only follow the pointer to the zero in the M matrix and terminate your traceback there only. This will prevent you from having alignments that are right-sided substrings.

#  

# Program specifications

### Input

The input to the program is a file with lines containing the following information.

1. A sequence of letters indicating sequence A
2. A sequence of letters indicating sequence B
3. An indication of whether local (1) or global (0) alignment is sought.
4. A set of gap penalties for introducing gaps into A or B.
5. The symbol alphabets to be used (e.g. ATGC for DNA strands and 21 single-letter abbreviations for the proteins, but there could be any alpha-numeric symbol (A-Z, a-z, 0-9))
6. Lines showing the score between an element in A and one in B.

An [example input file](https://canvas.stanford.edu/courses/122868/files/6438099/download?wrap=1):

```
AATGC
AGGC
0
0 0 0 0
4
ATGC
4
ATGC
1 1 A A 1
1 2 A T 0
1 3 A G 0
1 4 A C 0
2 1 T A 0
2 2 T T 1
2 3 T G 0
2 4 T C 0
3 1 G A 0
3 2 G T 0
3 3 G G 1
3 4 G C 0
4 1 C A 0
4 2 C T 0
4 3 C G 0
4 4 C C 1
```

And an annotated version of example input file, both of which can be found here in the [P1 student files directory.](https://canvas.stanford.edu/courses/122868/files/6502046/download) (This is just for your reference; your program will not need to read such files.)

 

```
AATGC          ;sequence A (length = # of columns)
AGGC           ;sequence B (length = # of rows)
0              ;0 = global,  1 = local
0 0 0 0        ;gap open penalty for A (dx), gap extension for A (ex), open for B (dy), extend B (ey)
4              ;number of letters in alphabet for sequence A = NA
ATGC           ;letters in alphabet for sequence A
4              ;number of letters in alphabet for sequence B = NB
ATGC           ;letters in alphabet for sequence B
1 1 A A 1      ;
1 2 A T 0      ;These are the entries in the match matrix between
1 3 A G 0      ; alphabet for A and alphabet for B.  Letters appear
1 4 A C 0      ; in same order as lines 6 and 8.
2 1 T A 0      ; There are NA * NB entries in this matrix, and each
2 2 T T 1      ;    entry is a row.
2 3 T G 0      ; First col is row number in match matrix
2 4 T C 0      ; Second col is col number in match matrix
3 1 G A 0      ; Third col is letter from alphabet A  
3 2 G T 0      ; Fourth col is letter from alphabet B
3 3 G G 1      ; Fifth col is match score between the two letters
3 4 G C 0      ; 
4 1 C A 0      ; Thus, this matrix gives positive score only for
4 2 C T 0      ; exact matches.  There are many other matrices
4 3 C G 0      ; also possible.
4 4 C C 1
```

 

For simplicity, you may assume that the input sequences will never be longer than 300 letters and the symbol alphabet size is never larger than 26.

 

### Output

The output of the program is a file with lines containing the following information. Please read carefully, as failure to comply with formatting requests will result in lost points.

 

1. Score for the best alignment (rounded to the first decimal place i.e. 3.1415->3.1).
2. All alignments which achieve the best score, with the format:
   1. empty line
   2. sequence A (with necessary gaps)
   3. sequence B (with necessary gaps) such that aligned characters from A and B are in the same column.

- Remember not to include the gaps at the very start and end of the alignment.
- Please use the underscore ("_") symbol for a gap as in the example output below.
- The order in which the alignments are output will not affect the autograder.

The resulting [output file](https://canvas.stanford.edu/courses/122868/files/6438094/download?wrap=1) for the above example:

```
3.0

ATGC
AGGC

AATG_C
A__GGC

ATG_C
A_GGC

AATGC
AG_GC

AATGC
A_GGC
```

### Command Line Arguments

Your program must be named align.py and take exactly 2 arguments: input_file output_file **exactly** as follows:

```
python3 align.py input_file output_file
```

 

**Failure for code to execute in exactly this manner will result in a grade of zero for the code portion of this project. Please note that you must execute this command from the terminal. Excessive print() statements will result in failure of code to run efficiently and will result in a grade of zero for the code portion of this project.**

 

[Links to an external site.](http://bmi214.stanford.edu/index.php?n=Site.Project1?action=sourceblock&num=7)

# Starter Code

Download the zipped code and example directory [**here**](https://canvas.stanford.edu/courses/122868/files/6502046/download). Starter code is provided in **align_skeleton.py**. This code will form the basis for you to build off of in your implementation and provides us with a framework to test functions and give you partial credit. Areas where there is flexibility in implementational decisions will be noted in comments and docstrings. Comfort with object-oriented programming and familiarity with dynamic programming will be necessary for this assignment. 

Make a copy of this called **align.py** and fill in the sections marked “FILL IN”.  

We have also provided a set of functions that may help you answer the quiz questions in **align_quiz_functions.py**.

*Unit Tests*

When coding up an algorithm, it is important to have test examples so you can check the behavior of the functions you write.

We have provided you with the skeleton of testing library to test your code in **align_test.py**. A couple unit tests are included to test your functions for reading in input, to test the match matrix class, and to test Ix updates. Fill in the unit tests, and use these to test that your program is working! Also think about how you might write additional unit tests; write one or two. What would be useful to test? Why? 

To run **align_test.py**, make sure it is in the same directory as your align.py file and the test_example.input are in the same directory. 

```
python3 align_test.py
```

It will tell you if the tests fail or pass.

*'**Note - you do not have to implement unit tests, but we strongly recommend you do. For grading purposes, we will be using these and similar tests to give partial credit. See note on grading at the end.***



*Debugging advice:*



1.   Write methods to print the score matrices (scores AND pointers) and the traceback. But make sure to turn these off at submit time.
2.   Work through a simple example by hand, make sure your code matches this!
3.   Write and use unit tests!



# Running the Alignments

You should test your program on several alignments, several of which are provided in the zipped file above. For provided input files in the **examples** directory, test your program and compare your output against the provided corresponding sample output files. **DO NOT MODIFY the example alignment files.** 

## Submission Instructions

1. Take the quiz on Canvas. We recommend you take the quiz after having written the code and use your code, but you can take the quiz *without* having finished. It involves understanding the algorithm and how it works on different inputs. The quiz will take some time to complete, please budget accordingly. 

   

2. Upload all source code to[ Gradescope (Links to an external site.)](https://www.gradescope.com/courses/186205) (your **align.py** and **align_test.py** files) as a zipped folder. Do not include irrelevant files.

 

# Rubric

| **Total 129 points**                                         |                                                              |                                                              |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Coding Score (total 85 points):                              | In order to give more partial credit, your coding score will be composed of two separate scores |                                                              |
|                                                              | 1. Raw accuracy score (Maximum number of points = 63)        | We will run your align program on 14 different inputs (some of these are held out for the set of inputs you were given) For each input you can receive a max of 4.5 points: 1 point for the correct score 1 point for correctly formatted score 1 point for the correct number of alignments 1 point for all correct alignments 0.5 point for some correct alignments |
|                                                              | 2. Implementation score (Maximum number of points = 22)      | Filled in local alignment score matrices correctly - 2 points **Correct score matrix updates - 8 points (3 each for Ix, Iy, 2 for M) **Make sure your end gap penalties are correct!Started the traceback in the correct place in global alignment - 4 points **Started the traceback in the correct place in local alignment - 4 points **Parameters are loaded correctly - 4 points** Hint: Think about how these will be tested with variations of the unit tests! Make sure yours work :) |
| Quiz questions score: 34 points                              |                                                              |                                                              |
| Style score: 10 points (see [Code Policy](https://canvas.stanford.edu/courses/122868/pages/code-policy) page) |                                                              |                                                              |

 

## Helpful Conceptual References

The following illustration is very helpful for gaining some visual intuition for the recursion formulas covered in class. Make sure you understand conceptually what's going on here and look at which arrow is pointing where and why.

![alignments](https://canvas.stanford.edu/courses/122868/files/6502010/download?wrap=1)

Additionally, here are two demos that allow you to play with some alignments: [Demo1 (Links to an external site.)](https://gtuckerkellogg.github.io/pairwise/demo/) and [Demo2 (Links to an external site.)](http://experiments.mostafa.io/public/needleman-wunsch/)

Finally, if you're struggling with the dynamic programming portion, the following blog post might be helpful: [Blog (Links to an external site.)](https://www.freecodecamp.org/news/demystifying-dynamic-programming-3efafb8d4296/)

# Key Concepts

Application

- Global vs local alignment
- End gap penalties
- Uses of sequence alignment

Implementation

- Ix and Iy:

- - Initialization of Ix, Iy
  - Arrow pointing from Ix to Iy

- Using appropriate end gap penalties

- Data structure use:

- - Pointers, matrices
  - Keeping track of match list

- Traceback

- - How/where to start the traceback
  - How to end the traceback
  - How to branch, update?

- - - Keep track of alignment so far, coords of next steps, go thru next steps
    - Keep track of completed alignments in some sort of data structure
    - Make sure to traceback through ties

- Right-sided substrings
- Global vs local alignment
- Initializing first row/col, terminating traceback in local alignment

# References

- [Mount Readings (Links to an external site.)](http://www.bioinformaticsonline.org/): Chapter 3 & pp 122-147, 240-259
- [Needleman-Wunsch global alignment (1970)](https://canvas.stanford.edu/courses/122868/files/6438131/download)[![Preview the document](https://canvas.stanford.edu/images/preview.png)](https://canvas.stanford.edu/courses/122868/files/6438131/download)
- [Gotoh alignment speedup article (1982)](https://canvas.stanford.edu/courses/122868/files/6438134/download)[![Preview the document](https://canvas.stanford.edu/images/preview.png)](https://canvas.stanford.edu/courses/122868/files/6438134/download)
- [Gribskov and Devereux's Sequence Analysis Primer (1994)](https://canvas.stanford.edu/courses/122868/files/6438136/download?wrap=1)[![Preview the document](https://canvas.stanford.edu/images/preview.png)](https://canvas.stanford.edu/courses/122868/files/6438136/download?wrap=1). 
- [Alignment algorithms from Durbin, et al.](https://canvas.stanford.edu/courses/122868/files/6438133/download?wrap=1)[![Preview the document](https://canvas.stanford.edu/images/preview.png)](https://canvas.stanford.edu/courses/122868/files/6438133/download?wrap=1)