# Introduction
This project implements the Burrows-Wheeler Transform (BWT) and uses it for exact pattern matching. The BWT is a reversible, block-sorting transformation. A special end-of-string character (`$`), which is lexicographically smaller than every other character, is appended to the text, and all cyclic rotations of the string are generated and sorted lexicographically to form the Burrows-Wheeler Matrix (BWM). The last column of the BWM is the BWT. It is a permutation of the original text that tends to group identical characters into long runs, which makes the text much easier to compress.

The BWT also underlies modern sequence alignment. Nearly all established short-read aligners use an index built on it (the FM-index) to map reads to a reference genome. Bowtie, for example, aligns short sequences to genome-sized references quickly while keeping a small memory footprint.

# Pseudocode
```
def BWT(string):
    IF string does not end with $:
        append $ to string
    FOR i from 0 to length of string - 1:
        build rotation = string[i:] + string[:i]     (move the first i chars to the back)
        store the pair (rotation, i)
    SORT the pairs alphabetically by rotation, keeping each i with its row
    FOR each (rotation, i) in sorted order:
        add the last char of rotation to the BWT
        add i to positions                            (this is the suffix array)
    RETURN BWT, positions


def suffix_array(string):
    FOR i from 0 to length of string - 1:
        store the pair (string[i:], i)                (suffix starting at i)
    SORT the pairs alphabetically by suffix
    RETURN the starting positions, in sorted order


def BWT_from_suffix_array(text, suffix_positions):
    start with an empty list for the last column
    FOR each pos in suffix_positions:
        IF pos == 0:
            add the LAST char of text ($)             (nothing comes before it, so wrap around)
        ELSE:
            add text[pos - 1]                         (the char just before this suffix)
    RETURN the list joined into a string


def cal_count(string):
    count how many times each char appears (freq)
    set cumulative = 0
    FOR each char in sorted order:
        IF char is $:
            skip it                                   (so count of the first letter is 0)
        count[char] = cumulative                      (number of letters smaller than char)
        cumulative += freq[char]
    RETURN count


def cal_occur(bwt_string):
    FOR each distinct char:
        make a list of 0s, one per position
    FOR each position i and its char current:
        FOR each char in the dict:
            IF i > 0:
                copy occur[char][i - 1] into occur[char][i]
        add 1 to occur[current][i]
    RETURN occur                                      (occur[char][i] counts bwt[0..i], inclusive)


def update_range(lower, upper, count, occur, a):
    IF lower == 0:
        a_before_lower = 0
    ELSE:
        a_before_lower = occur[a][lower - 1]          (a's before the range)
    a_through_upper = occur[a][upper]                 (a's through the end of the range)
    lower_new = count[a] + a_before_lower + 1
    upper_new = count[a] + a_through_upper
    RETURN lower_new, upper_new


def find_match(query, reference):
    IF query is empty:
        RETURN empty list
    IF reference does not end with $:
        append $ to reference
    suffix_positions = suffix_array(reference)
    transformed = BWT_from_suffix_array(reference, suffix_positions)
    count = cal_count(reference)
    occur = cal_occur(transformed)
    set lower = 0, upper = length of transformed - 1  (every row matches the empty string)
    FOR each char in query, from last to first:
        IF char is not in count:
            RETURN empty list                         (char never appears in reference)
        lower, upper = update_range(lower, upper, count, occur, char)
        IF lower > upper:
            RETURN empty list                         (no match)
    collect suffix_positions[row] for every row from lower to upper
    SORT the positions                                (rows are alphabetical, not text order)
    RETURN positions


def run_length_encode(bwt_string):
    start with empty encoded text
    IF bwt_string is empty:
        RETURN empty encoded text
    current_char = first char, run_length = 1
    FOR each remaining char:
        IF char == current_char:
            run_length += 1
        ELSE:
            add current_char + run_length to encoded text
            current_char = char, run_length = 1
    add the final current_char + run_length to encoded text
    RETURN encoded text


def run_length_decode(encoded):
    start with empty decoded text, i = 0
    WHILE i < length of encoded:
        symbol = encoded[i], move i forward one
        read all following digits into a number
        add symbol repeated that many times to decoded text
    RETURN decoded text
```

# Successes
**Collaboaration and Communication:** One of the major things that made this project successful was our ability and willingness to meet sychronously a few times to talk through the parts we weren't understanding and bounce ideas off of each other. We all contributed different ideas for implementations, and worked together through various versions of our code to make sure we all agreed on which implementation strategy made the most sense and which would work best for us. 


# Struggles
**Initial Setup:** We had some struggles with the `def BWT()` function. Initially, the code included `BWT += r[0][-1]` which read the first column instead of the last one. We also initially struggled with understanding how to utilize the two count dictionaries in the suggested `def cal_count()` default pseudocode. After looking into the `Counter()` dictionary import, it made sense. 

**Pattern Matching:** This seemed to be the hardest part for our group. We initially had some confusion differentiating the `occur` and `count` dictionaries. We had to go back into some class notes as well as the HTML to get visuals and explanations on what they each were supposed to contain.

# Personal Reflections
## Yulia
This project was difficult to wrap my head around conceptually, but when my group members and I met to go through the html examples we were able to get a good grasp of what should be happening step-by-step, which made writing the pseudocode much easier. When reading the default pseudocode given by the assignment, we were confused on some aspects and wondered if they may be intentional errors. Using these inconsistencies, we came up with different ideas, not all of which ended up being implemented. Also, this assignment did not have a big driver program calling all the functions together, so we opted to include test code chunks to make sure it is running as expected. Overall we made good progress and the concepts made sense to me!


## Kailey
Understanding the BWT algorithm was difficult conceptually. The toughest part for me to comprehend was pattern matching. A key aspect of this project was reviewing and comprehending the concept with my group. We met with the time we all could meet to create a consensus on our understanding of the algorithm and how we should approach this project. We were able to walk-through the entirety of the .html provided to us and ask questions when there was confusion and plan how we would approach constructing the code. In addition to our own pseudocode, having the pseudocode provided, created a supplementary learning material that would provide the basis for our code. I felt like I was overly ambitious in the beginning of the project, thinking we would require complex data structures like NumPy arrays to configure the data; however after taking some time to process the material, I was able to realize that this step was possibly not needed for this project. With reading more about NumPy and its capabilities and thinking about larger datasets and scalability, using arrays would be helpful in this aspect. 

This project has helped me develop a new learning approach to tackling these projects. I first learned that I need to comprehend the algorithm and perform an example by hand if possible. Then write out the pseudocode. From the pseudocode, I can go line by line constructing code and use test code to check my progress! Finally, I try to include biological considerations if I am able to code such complexities. 

## Maggie
As expected, we had to spend a lot of time working together to conceptually understand this project. I think we were all initilially on the same page in our understanding: we understood the beginning part of this program, but were confused on the pattern matching aspect. We spent a few meetings just oing over pseudocode and the HTML provided to try to wrap our heads around it. Like Yulia mentioned, we had some trouble understanding the initial pseudocode that was given, and had to parse through that piece by peice. We ended up using the given pseudocode more as a guide for our new ideas, which seemed to work pretty well. OUr collaboratin was vital to me understanding the concepts here, and I feel like we ended in a good spot of all of us understanding the concept and the code we wrote.

# Generative AI Appendix
Claude Opus 5.5 was used to reformat/correct our original pseudocode to match the code changes.
