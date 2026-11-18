
## Project 8: Viterbi Decoding for Hidden Markov Models

# Introduction
This project implements the Viterbi algorithm for decoding hidden Markov models (HMMs) in Python. Given an HMM and a sequence of observations, it returns the most probable sequence of hidden states and that path's log probability. The implementation works for any HMM whose parameters are supplied by the user and does not depend on the example data.

In an HMM, the states are hidden and only the symbols they emit are observed. A gene finder never sees the labels "exon" and "intron," only the nucleotides, so the labels must be inferred (Eddy, 2004).

An HMM is specified by:
- a set of **hidden states**;
- an **alphabet** of observable symbols;
- **initial probabilities**: the probability of starting in each state;
- **transition probabilities**: the probability of moving from each state to each state (a Markov chain on the hidden states);
- **emission probabilities**: the probability of each state emitting each symbol.


**Running example: CpG islands**
In vertebrate genomes, the dinucleotide CG is rarer than base composition alone would predict, except in short stretches called CpG islands, which often mark promoters. A two-state HMM with an "island" state and a "genome background" state, each with its own nucleotide composition, is the classic teaching model for segmenting DNA (Durbin et al., 1998). The example data in this repository follow that model.


**The Viterbi algorithm**

Andrew Viterbi (1967) devised his algorithm to decode convolutional codes in digital communication. Forney (1973) later recognized it as a general solution for finding the most probable path through a hidden Markov process. It answers the question:

> Given the observations, which single sequence of hidden states has the highest joint probability with them?

Viterbi decoding is how an HMM turns a long, unlabelled sequence into an annotation. GENSCAN models human gene structure (exons, introns, splice sites, untranslated regions) as a probabilistic state machine and reports the most probable gene structure for a stretch of genomic DNA (Burge & Karlin, 1997). ChromHMM segments the genome into chromatin states, such as "active promoter" or "repressed," from combinations of histone marks (Ernst & Kellis, 2010, 2012). In both cases, the labelled genome browser tracks that biologists rely on come from a decoding step like the one implemented here.


**Features**
- **Optimal path finding:** returns the single most probable hidden-state sequence.
- **Dynamic programming:** fills a T × N score table (T positions, N states) in O(T·N²) time and O(T·N) space, avoiding redundant computation.
- **Log-space computation:** prevents numerical underflow on long sequences.
- **Traceback:** reconstructs the optimal path from stored back-pointers.
- **Input validation:** checks that states are consistent across all parameters and that every probability distribution sums to 1.


**Why log space?**
A path probability is a product of many numbers less than 1. Multiplying many of these numbers cause the result to drop below the smallest possible 64-bit number (underflow), and causes every path to incorrectly look equally impossible. Working with logarithms turns products into sums, which keep them in range for comparison. 


**Algorithm overview**
1. **Initialization**
   - Build the log-probability matrix from the log initial probabilities plus the log emission probability of the first observation.
   - Initialize a traceback matrix for path reconstruction.
2. **Recursion**
   - For each position and each state, find the previous state that maximizes the score, then add the log transition and log emission terms.
   - Store both the scores and the traceback pointers.
3. **Termination and traceback**
   - Select the final state with the highest log probability.
   - Follow the traceback pointers backward to reconstruct the optimal path.
  
**Our main dataset used**
We implemented Viterbi by attempting to use the entire chromosome 22 sequence, consisting of 50,818,468 bases. Using a UCSC derived file, cpgIslandExt.txt, we extracted the CpG regions of chr22 and used this to define our hidden states (CpG vs Genomic) and calculate its initial, emission, and transition probabilities to run the chr22 sequence with. 

# Pseudocode
`train_cpg_hmm` estimates the HMM parameters from a chromosome FASTA file and the UCSC CpG island annotation. The two states are `I` (CpG island) and `G` (genomic background).

`_checksum` is a static method used to make sure all the probabilities (ie each transition probability) adds up to 1. It will raise a `ValueError` if they do not.

`__init__` initializes the HMM object and sets parameters. There is input erro checking in here, and `ValueError` will be raised if there is an error.

`_viterbi_recursion` uses the current state to find the previous state that maximizes the score, and stores them as pointers for traceback.

`viterbi` sets up the probability matrix with initial probabilties and the first observation, and initializes a traceback matrix.

The pseudocode for this project is as follows:

```
FUNCTION train_cpg_hmm(fasta_path, cpg_path, chrom, pseudocount)
    # load sequence
    seq <- concatenate all non-header lines of fasta_path, uppercased
    bases <- seq as an array of single characters

    # label each position
    is_island <- array of False, one entry per base
    for each line in cpg_path:
        if line's chromosome = chrom:
            is_island[start : end] <- True

    valid_base <- positions where base is A, C, G, or T (excludes N)
    mask["I"] <- valid_base AND is_island
    mask["G"] <- valid_base AND NOT is_island

    NORMALIZE(counts): return counts divided by sum of counts

    # initial probabilities
    for each state: state_counts[state] <- number of positions in mask[state]
    init_probs <- NORMALIZE(state_counts)

    # emission probabilities
    for each state:
        for each base in A, C, G, T:
            base_counts[base] <- count of base within mask[state] + pseudocount
        emit_probs[state] <- NORMALIZE(base_counts)

    # transition probabilities
    both_valid <- positions t where base[t] and base[t+1] are both valid
    for each such t:
        record (is_island[t], is_island[t+1]) as a transition
    count G->G, G->I, I->G, I->I, adding pseudocount to each
    trans_probs["I"] <- NORMALIZE({I: I->I, G: I->G})
    trans_probs["G"] <- NORMALIZE({I: G->I, G: G->G})

    return init_probs, trans_probs, emit_probs


CLASS HiddenMarkovModel
FUNCTION _checksum(name, prob)
    total <- sum of all values in prob
    if total is not close to 1.0 (tolerance 1e-9):
        raise error "name must sum to 1, got total"


FUNCTION __init__(init_probs, trans_probs, emit_probs)
    store init_probs, trans_probs, emit_probs
    states <- list of keys of emit_probs

    # validate inputs
    if init_probs, trans_probs, or emit_probs is empty:
        raise error

    if the state sets of init_probs, trans_probs, and emit_probs differ:
        raise error

    _checksum("init_probs", init_probs)
    for each (state, row) in trans_probs:
        _checksum("trans_probs[state]", row)
    for each (state, row) in emit_probs:
        _checksum("emit_probs[state]", row)

    for each (state, row) in trans_probs:
        if keys of row do not match states:
            raise error "trans_probs[state] references unknown state"


FUNCTION _viterbi_recursion(seq, prob_matrix, traceback)
    next_pos <- number of rows in prob_matrix

    # base case
    if next_pos = length of seq:
        return prob_matrix, traceback

    prev <- prob_matrix[next_pos - 1]
    probs <- empty map
    pointers <- empty map

    for each state s in states:
        # find the best previous state
        best_prev <- None
        best_score <- -infinity
        for each state p in states:
            score <- prev[p] + log(trans_probs[p][s])
            if score > best_score:
                best_prev <- p
                best_score <- score

        # apply transition and emission
        probs[s] <- prev[best_prev]
                    + log(trans_probs[best_prev][s])
                    + log(emit_probs[s][seq[next_pos]])
        pointers[s] <- best_prev

    append probs to prob_matrix
    append pointers to traceback

    return _viterbi_recursion(seq, prob_matrix, traceback)


FUNCTION _viterbi_traceback(traceback, last_state)
    path <- [last_state]
    for t from (length of traceback - 1) down to 1:
        append traceback[t][last element of path] to path
    reverse path
    return path joined into a string


FUNCTION viterbi(seq)
    # initialization
    first_row <- empty map
    for each state s in states:
        first_row[s] <- log(init_probs[s]) + log(emit_probs[s][seq[0]])
    prob_matrix <- [first_row]
    traceback <- [ map of each state to None ]

    # recursion
    prob_matrix, traceback <- VITERBI_RECURSION(seq, prob_matrix, traceback)

    # termination
    last_row <- last row of prob_matrix
    last_state <- state s in states with the highest last_row[s]
    log_prob <- last_row[last_state]

    # traceback
    path <- VITERBI_TRACEBACK(traceback, last_state)

    return path, log_prob
```

Our test code follows a toy dataset we created:
```
# Example initial probabilities (probability of starting in each state)
init_probs = {
    "I": 0.2,
    "G": 0.8
}

# Example transition probabilities (probability of moving from one state to another)
trans_probs = {
    "I": {"I": 0.7, "G": 0.3},
    "G": {"I": 0.1, "G": 0.9}
}

# Example emission probabilities (probability of observing a symbol in a given state)
emit_probs = {
    "I": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
    "G": {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
}

toy_hmm = HiddenMarkovModel(init_probs, trans_probs, emit_probs)
print(f"State path for input 'GGCACTGAA': {toy_hmm.viterbi("GGCACTGAA")}")
print(f"State path for input 'ATATCGCGCGCGAT': {toy_hmm.viterbi("ATATCGCGCGCGAT")}")
```
With the expected output:
```
State path for input 'GGCACTGAA': ('GGGGGGGGG', -13.929108456051065)
State path for input 'ATATCGCGCGCGAT': ('GGGGIIIIIIIIGG', -21.20203079978548)
```


# Successes
**Class construction:** We successfully created the class `HiddenMarkovModel` and were able to implement different functions and attributes using `self`

**Recursion:** We were able to implement recursion in our `_viterbi_recursion()` function. THis helped us find the position of the previous state with the maximum score, and store both scores as traceback pointers. 

**Log-space computation:** We successfully implemented log-space calculations to avoid the problem of underflow in out calculations.

**Communication and workflow:** As a team, we were able to use bith synchronous and asynchronous work sessions effectively. We used our synchrinous team meetings to understand the algorithm, work through pseudocode, and coordinate schedules and timing for peices of the project. We then were all ab le to effectively work asynchronously through GitHub and kept communcation through Teams messages

# Struggles

First error we encountered below when running the toy dataset, the error check loops over the keys only, so it crashes. We modified the code to use sets: `if set(row) != set(self.states):`
```
ValueError                                Traceback (most recent call last)
Cell In[7], line 21
     17     "I": {"A": 0.1, "C": 0.4, "G": 0.4, "T": 0.1},
     18     "G": {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
     19 }
     20 
---> 21 toy_hmm = HiddenMarkovModel(init_probs, trans_probs, emit_probs)
     22 print(f"State path for input 'GGCACTGAA': {toy_hmm.viterbi("GGCACTGAA")}")
     23 print(f"State path for input 'ATATCGCGCGCGAT': {toy_hmm.viterbi("ATATCGCGCGCGAT")}")

Cell In[6], line 51, in HiddenMarkovModel.__init__(self, init_probs, trans_probs, emit_probs)
     47         for state, row in emit_probs.items():
     48             self._checksum(f"emit_probs[{state!r}]", row)
     49 
     50         # trans_probs references unknown state
---> 51         for state, row in trans_probs:
     52             if row.key() != self.states:
     53                 raise ValueError("trans_probs[{state}] references unknown state")
     54         #-----------------------------------

ValueError: not enough values to unpack (expected 2, got 1)
```

The next major error message we got was in actually running the Viterbi on the chr22 sequence, after successfully loading the model with calculated initial, emission, and transition probabilities for chr22. It raises a key error upon encountering "N" bases in the sequence string, so we needed to account for that possibility in the Viterbi as well. Our next steps are to optimize to handle that case by writing a message in the error checking section notifying when the input sequence has bases other than A C G and T, and in the same check run Viterbi with the non-N chunks of sequence, retaining the correct positions.
```
{
	"name": "KeyError",
	"message": "'N'",
	"stack": "---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
Cell In[13], line 28
     24     sequence = \"\".join(line.strip() for line in fasta_file if not line.startswith(\">\")).upper()
     25 print(f\"\
Loaded chr22: {len(sequence):,} bases\")
     26 
     27 # Run Viterbi on the whole chromosome
---> 28 path, log_prob = chr22_hmm.viterbi(sequence)

Cell In[8], line 105, in HiddenMarkovModel.viterbi(self, seq)
    101         # Set up probability matrix using initial probabilities and the first observation
    102         # Initialize a traceback matrix for path reconstruction
    103         first_row = {}
    104         for s in self.states:
--> 105             first_row[s] = math.log(self.init_probs[s]) + math.log(self.emit_probs[s][seq[0]])
    106         prob_matrix = [first_row]
    107 
    108         traceback = [{s: None for s in self.states}]

KeyError: 'N'"
}
```


# Personal Reflections
## Yulia
Starting this notebook with my group members from scratch with no given datasets was initially challenging, as we had to chat through all the components that we would need. First we started with conceptually understanding the Viterbi lecture notes. Then, understanding what types of biological datasets could work well with the Viterbi model, and realizing that we had to make the Viterbi class very generalizable to handle any sorts of input data and N number of hidden states. We decided to implement a toy dataset, and a chromosome 22 dataset using the CpG and Genomic states. I primarily worked on the function to incorporate the chr22 sequence, define the states and the probabilities to pass it into the model. It took a while for me to understand how to calculate emission and transition probabilities, so I had to use outside resources for that. We were not able to run the model on the chr22 data since we did not handle the N stretches in the sequence itself, although we did in calculating probabilities, so that will be our next step. I look forward to implementing the next steps of the Hidden Markov Model we will learn about, and working with other biological datasets that fit the Viterbi model- so that we can further optimize our code to be generalizable and handle different types of cases. 

## Maggie
Initially, I was feeling a little overwhelmed since we had no staring notebook to give us an outline of what this project should look like. It took some time of all three of us sitting down and working through the algorithm for me to feel a little less overwhelemed. We took a good amount of time in the beginning just sitting down with our notes and trying to comprehend the algorithm. When putting it together, Yulia found a dataset from chromosome 22 that we were able to implement and sue for running the program. We ended up not being able to completely run the dataset, because we ran into a problem with stretches of character `N` that we did not have time to fix. This would be something that we as a roup would like to work on handling for the next assignments. Overall we worked well together, and made sure that ecveryone was understanding the different aspects of this project.

## Other member
Other members' reflections on the project

# References
Durbin, R., Eddy, S. R., Krogh, A., & Mitchison, G. (1998). Biological sequence analysis: Probabilistic models of proteins and nucleic acids. Cambridge University Press. https://doi.org/10.1017/CBO9780511790492

Eddy, S. R. (2004). What is a hidden Markov model? Nature Biotechnology, 22(10), 1315–1316. https://doi.org/10.1038/nbt1004-1315

Ernst, J., & Kellis, M. (2010). Discovery and characterization of chromatin states for systematic annotation of the human genome. Nature Biotechnology, 28(8), 817–825. https://doi.org/10.1038/nbt.1662

Ernst, J., & Kellis, M. (2012). ChromHMM: Automating chromatin-state discovery and characterization. Nature Methods, 9(3), 215–216. https://doi.org/10.1038/nmeth.1906

Forney, G. D., Jr. (1973). The Viterbi algorithm. Proceedings of the IEEE, 61(3), 268–278. https://doi.org/10.1109/PROC.1973.9030 

Viterbi, A. (1967). Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. IEEE Transactions on Information Theory, 13(2), 260–269. https://doi.org/10.1109/TIT.1967.1054010

# Generative AI Appendix
Claude Opus 5.5 was used to help construct the `train_cpg_hmm()` function, with input questions: "What is the best way to process and store sequence data for use in hidden states and calculating probabilities?" "How do I calculate emission probabilities?" "How do I calculate transition probabilities?". From these results, suggestions to store the sequence in Numpy array as bytes was implemented, alongside adding in pseudocounts to prevent any probability from being 0 which would fail in Viterbi log calculation, and calculating transition "codes". 
