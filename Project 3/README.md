# Project 3 Instructions

### Important Information

1. Project 3 is due by **5:00pm on October 30, 2020**. Late policy applies up to 3 days as specified in the syllabus.

1. Programming projects must be completed individually. You may discuss algorithms with others, but the coding should be done alone. You must explicitly name everyone with whom you discussed this project in the header comments of your code and in the collaboration attestation in the project quiz. Students must abide by the terms of the Stanford Honor Code.
2. Remember to consult Piazza, as many common questions will be asked and answered there.
3. Prior to beginning the assignment, please read the [Code Policy](https://canvas.stanford.edu/courses/122868/pages/code-policy).
4. Read this project page carefully; your source code must run **exactly** as specified below. We will be unit testing your code on Gradescope.

### Key Dependencies

You can download all files required for the project [here](https://canvas.stanford.edu/courses/122868/files/folder/p3).

**PyRosetta**: Many of the utility functions in this project require PyRosetta, a Python wrapper for the Rosetta software suite. You will need to request an academic license for free at this link: 

[https://els.comotion.uw.edu/express_license_technologies/pyrosetta (Links to an external site.)](https://els.comotion.uw.edu/express_license_technologies/pyrosetta).

Once approved, download the latest PyRosetta from [http://www.pyrosetta.org/dow (Links to an external site.)](http://www.pyrosetta.org/dow). Make sure you select the right distribution for your machine (Linux/MacOS vs. Windows). Follow the appropriate installation instructions at the bottom of the page.

**Make sure to do this step as early as possible—installing PyRosetta may take a while\**

 

**Biopython:** The other dependency you will need is [Biopython (Links to an external site.)](https://biopython.org/wiki/Documentation), which is useful for interacting with PDB files. This can be installed using `pip install biopython`.

### Background

Proteins are extraordinary molecular machines that perform and mediate fundamental functions in living organisms, including molecular signaling, reaction catalysis, and energy transformation. This remarkable diversity in function of natural proteins is made possible by variations in their three-dimensional structures, each of which has been optimized through evolution to perform its particular function. It is also known that a protein’s structure is determined by its amino acid sequence, and that the folded conformation of a protein corresponds to the minimum-energy state for its sequence. Thus, the ability to predict the structure of a protein from its amino acid sequence alone would be very useful for addressing many biological problems, such as predicting protein function from genomic data and rationally engineering novel proteins to perform a desired function. For this reason, the problem of **protein structure prediction** has been a holy grail in structural bioinformatics for decades.

We will focus on *ab initio* structure prediction, which aims to fold a protein from sequence alone without using a template structure. In this project, you will implement a fragment assembly algorithm, one of the most widely-used and successful approaches for *ab initio* protein structure prediction. This is a simplified version of the Rosetta protein folding algorithm described in lecture, which uses a Monte Carlo simulation with simulated annealing to iteratively assemble fragments from known protein structures into a novel structure that matches the input sequence.

### Project Overview

In this assignment, you will fold some small target (or "native") proteins. We give you data for three proteins:

1. Helix – this is a simple 17-residue alpha helix. We **strongly** recommend you use this structure for debugging and testing your code, as it should take less than 5 minutes to fold.
2. 1FW4 – this is the crystal structure of the C-terminal domain of calmodulin, an important calcium binding signaling protein.
3. 1UBQ – this is the crystal structure of ubiquitin, a small regulatory protein that is present in almost all eukaryotic tissues.

Implement the fragment assembly protocol for *ab initio* protein folding by following the implementation details below. In each simulation, you will fold the protein from an extended configuration using two stages: (1) a coarse assembly stage, which samples fragments of length 9, and (2) a refinement stage, which samples fragments of length 3. After the full simulation, you will perform an energy minimization and report (1) the final energy of the folded conformation and (2) the RMSD to the target structure (in Angstroms).

First, make sure your program works consistently on *helix.pdb*. The total simulation should take 5 minutes or less on your local machine and achieve a final RMSD of less than 1Å to the target structure.

Then, run 10 simulations each for 1FW4 and 1UBQ. Since these proteins are slightly larger, each simulation could take 5-10 minutes. **With this in mind, make sure you leave enough time to run all of these simulations.** If it is taking longer than 10 minutes per simulation, think about how you could implement your code more efficiently—the most common reason for inefficient simulations is unnecessarily repeating operations during sampling.

Finally, visualize the lowest-RMSD structures in PyMol (if you haven't used this before, see Appendix).

### Data Formats

We include data files in the following formats:

**.fasta files:** These files contain protein sequences, represented as a string of 1-letter amino acid codes. The lines starting with the symbol “>” simply denote the protein and chain that the sequence comes from.

**.pdb files:** This file format is used to store protein structures derived from X-ray crystallography data. These files are fixed-width and contain metadata information as well as the identity and xyz-coordinates of every atom in the protein. More information on the PDB file format is available on [Wikipedia (Links to an external site.)](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)), but this is not necessary to complete the project; we provide functions to interact with PDB files as necessary.

**.frag files:** The fragment assembly protocol you will be implementing requires sampling from a “fragment library”, which contains a set of candidate fragments for each sequence position in the input protein. These fragment libraries are generated by comparing short windows in the input sequence with known protein structures using PSI-BLAST (a variation of BLAST for protein sequences) and secondary structure prediction. Based on these comparisons, the most similar fragments (typically the top 200) are generated for each position. We have pre-computed 3- and 9-residue fragment libraries and provide them to you in Rosetta fragment file format (see [here (Links to an external site.)](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/fragment-file) for details). The only important data are the **phi** and **psi** columns, but it is useful to be familiar with the other information contained in these files.

**.rmsd files:** Typically, in *ab initio* folding you do not have access to the target structure (for obvious reasons), so you must sample randomly from the fragment library at each position. However, for this assignment you will be implementing *biased* forward folding, in which we reduce our sample space of fragments to only those which are closest in structure to the target fold. We have pre-calculated the root-mean-square deviation (RMSD) between each fragment in the fragment library and its corresponding fragment in the native structure. These files are tab-delimited with the following three columns:

```
position          fragment          rmsd
```

You do not need to modify these files, although the script used to generate the .rmsd files is given to you in *code/get_lowRMS_frags.py* if you are interested.

 

### Implementation Overview

1. ### Monte Carlo sampling with simulated annealing

**Monte Carlo fragment sampling**

This is probably the most widely used and successful strategy for folding a protein *ab initio*. This approach is so effective because it provides a **fast**, **knowledge-based** strategy for searching the space of possible structures for a protein. This conformational space is intractably large, especially for longer sequences, so a naïve sampling approach that permutes the protein structure at random is likely to take a very long time to reach a “good” structure. However, we can leverage the more than 150,000 known protein structures in the PDB to selectively sample regions of the space that we predict to be close to the target structure. This results in **fragment libraries** that are customized for the protein of interest, containing short sections of structures that have similar sequences to short windows of the input sequence. By sampling these fragments, we are able to effectively “jump” rapidly between different regions of the conformational space in search of low-energy configurations. In short, the procedure consists of a series of steps, or moves, where each move consists of two stages: (1) a fragment-based move in the conformational space and (2) an evaluation of the favorability of the move.

A typical *ab initio* folding protocol, such as Rosetta, starts with a round of coarse fragment assembly using 9-residue fragments. Each Monte Carlo move proceeds as follows:

1. Sample a random position in the protein chain. Note that this position corresponds to the **left-hand side** of the fragment window. Think carefully about which positions in the chain you can actually sample from; there are restrictions imposed by the fragment length and the definition of a torsional angle (see illustration [here](https://ww2.chemistry.gatech.edu/~lw26/structure/protein/secondary_structure/phi_psi/right.html)).
2. Retrieve list of candidate fragments for that position (see below).
3. Sample a random fragment from this list of candidates
4. Insert the sampled fragment into the current protein chain at the sampled position **in torsion space.** In other words, replace each torsion angle (![LaTeX: \phi](https://canvas.stanford.edu/equation_images/%255Cphi)and ![LaTeX: \psi](https://canvas.stanford.edu/equation_images/%255Cpsi)) in the current protein chain with corresponding torsion angles from the selected fragment.
5. Score the resulting conformation and accept or reject using the Metropolis criterion (see below).

After this coarse 9-mer assembly, we perform a refinement stage that uses the same Monte Carlo procedure to sample from 3-mer fragments. This enables a finer adjustment of the torsions in small regions of the protein, which provides more precision in regions such as small loops.

**Choice of scoring function**

After making a move, we then need to evaluate whether the move was beneficial. To do this, we must define a scoring function. Since it has been shown that a protein sequence will fold into the conformation with the lowest free energy, the scoring function should represent the energy of the structure as well as possible; we call this the **energy function.**

For maximum accuracy, we could calculate energy using every atom in the simulation. However, this is extremely expensive and the full-atom conformational space is not very smooth, making Monte Carlo more difficult. To alleviate these issues, rather than calculating energy over every atom in the protein, we use a **centroid** representation of each amino acid. In a centroid representation, the backbone atoms (N, C, C, O) remain the same, but the side chain is simplified to a single pseudo-“atom”, whose position is determined by the centroid of the side chain and whose radius and atomic properties (charge, polarity, etc.) are determined by the residue’s identity. This not only makes energy calculation faster, but also smooths out the energy landscape to aid in sampling. In this project, we use the Rosetta ‘score3’ centroid energy function, which is largely knowledge-based (i.e. derived from statistics calculated over known structures, rather than physics equations). See Table 1 in this [paper (Links to an external site.)](https://www-sciencedirect-com.stanford.idm.oclc.org/science/article/pii/S0076687904830040) for a detailed description of each energy term. You will not be asked to calculate any of these terms, but you will need to have a conceptual understanding of their meaning for the quiz.

![Screen Shot 2019-10-15 at 5.29.19 PM.png](https://canvas.stanford.edu/users/186210/files/5031348/preview?verifier=PV5dngSkwh4cmk63suOefb6raOf7bnWbN3ZEo7Vp)

Visualization of all-atom vs. centroid representations. ([Kmiecik et al., 2016](http://pubs.acs.org/doi/pdf/10.1021/acs.chemrev.6b00163))

**Metropolis criterion**

After making a single move (i.e. permuting a fragment) and measuring the energy of the resulting structure, you need a way to decide whether or not to accept the move. Since we are seeking the global energy minimum, we generally want to accept moves that decrease energy and reject moves that increase it. However, a naïve solution that simply rejects all moves that increase energy will by construction force the simulation to find the nearest energy minimum, regardless of whether it is a global minimum. This results in sub-optimal structures for most simulations. We can relax this criterion by accepting some moves that increase energy, allowing the simulation to escape local minima and explore the energy landscape. Specifically, we accept moves with a probability that depends on how much the energy increases using the **Metropolis criterion:**

![equ1](equ1.svg)

where ![LaTeX: \Delta E=E_{after}-E_{before}](https://canvas.stanford.edu/equation_images/%255CDelta%2520E%253DE_%257Bafter%257D-E_%257Bbefore%257D) is the change in energy produced by the move, k is the Boltzmann constant (k = 1 for this project), and T is the temperature of the system. Using this criterion, we guarantee that if we run the simulation for long enough (i.e. as time goes to infinity), the probability of observing any particular configuration is given by the Boltzmann distribution

![equ2](E:\GoogleDrive\Courses\BIOMEDIN214\Project 3\equ2.svg)

thus ensuring that our sampling is physically accurate.

**Simulated annealing**

An important consideration when evaluating moves with the Metropolis criterion is the temperature parameter T. High temperatures result in almost all moves being accepted, while low temperatures ensure that only decreases or *very small* increases in energy are accepted. An effective way to perform energy minimization is by starting at a high temperature and gradually reducing the temperature during the course of the simulation. This procedure is known as **simulated annealing**.

In this assignment, you will use an exponential annealing schedule with an annealing rate of 0.999:

![LaTeX: T_t=T_{start}\cdot0.999^t](https://canvas.stanford.edu/equation_images/T_t%253DT_%257Bstart%257D%255Ccdot0.999%255Et)

so that the temperature at each step t is given by

![LaTeX: T_t=0.999\:\cdot\:T_{t-1}](https://canvas.stanford.edu/equation_images/T_t%253D0.999%255C%253A%255Ccdot%255C%253AT_%257Bt-1%257D)

Note that the temperature is annealed after every move that is **accepted**, not every move attempted. We also ignore the Boltzmann constant k for simplicity (i.e. k = 1). Here, we will use the following annealing schedules:

- For the 9-mer fragment assembly stage, anneal from ![LaTeX: T_{start}=100](https://canvas.stanford.edu/equation_images/T_%257Bstart%257D%253D100) to ![LaTeX: T_{end}=1](https://canvas.stanford.edu/equation_images/T_%257Bend%257D%253D1)
- For the 3-mer refinement stage, anneal from ![LaTeX: T_{start}=1](https://canvas.stanford.edu/equation_images/T_%257Bstart%257D%253D1)to ![LaTeX: T_{end}=0.1](https://canvas.stanford.edu/equation_images/T_%257Bend%257D%253D0.1)

**Energy minimization (relaxation)**

The result of the first two stages is a folded protein, but it is essentially a complex combination of a large number of fragments from other proteins and is not required to be globally consistent. For that reason, we perform energy minimization on the output structure. To do this, we convert the protein structure back to a full-atom representation and run Rosetta’s [FastRelax](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/relax) minimization protocol. FastRelax efficiently optimizes the side chain configurations of the protein while allowing small changes to the overall backbone in order to find a minimum-energy structure. You do not need to implement this; we provide a **relax** function in **utils.py** that does this for you.

**Evaluation**

After minimization, the main evaluation metrics for the final folded protein are **energy** and **root-mean-square deviation (RMSD)** from the target structure. To calculate RMSD, we superimpose the predicted structure onto the target structure and measure RMSD between only the C atoms of the two structures:

![LaTeX: RMSD\left(A,\:B\right)=\sqrt{\frac{1}{N}\sum_i\left(C\alpha_{i,A}-C\alpha_{i,B}\right)^2}](https://canvas.stanford.edu/equation_images/RMSD%255Cleft(A%252C%255C%253AB%255Cright)%253D%255Csqrt%257B%255Cfrac%257B1%257D%257BN%257D%255Csum_i%255Cleft(C%255Calpha_%257Bi%252CA%257D-C%255Calpha_%257Bi%252CB%257D%255Cright)%255E2%257D)

where A and B are the structure to be compared, and N is the number of residues in each structure. The superimposition and RMSD calculation is also implemented for you in the **utils.relax** function.

Example superposition of a protein before and after minimization

![Screen Shot 2019-10-15 at 12.46.49 PM.png](https://canvas.stanford.edu/users/186210/files/5031399/preview?verifier=b3NT15EBpP8k6ZY3EtgeVJqFjBJUJ7L5B6gBfdRr)

 

2. ### Implementation Instructions

Putting everything together, you will implement a Monte Carlo fragment search with simulated annealing as follows, using the class provided in **FragmentSampler.py**:

1. Assemble 9-mers using MC procedure from T = 100 to T =1
   1. At each step:
      1. Sample random 9-residue window in sequence
      2. Get list of **top N** candidate fragments at that position
      3. Sample random candidate fragment
      4. Replace torsion angles of selected 9-residue window with torsions from selected fragment. Note that you will have to copy the protein here, because if the insertion is rejected it must be returned to its previous state.
      5. Measure energy and accept/reject using Metropolis criterion
      6. If accept, anneal T and move to new protein position. If reject, sample new fragment.
   2. Return the **best** (lowest-energy) structure produced during this simulation. **This may or may not be the last structure generated.**
   3. During the simulation, write the iteration, temperature, and energy to a log file for later reference
2. Assemble 3-mers using MC procedure from T = 1 to T = 0.1, starting with the best structure from stage 1.
   1. Do the same as above, but using 3-mer fragments
   2. **Again, keep track of the best (lowest-energy) structure produced during this simulation.** Save this in a PDB file called *best.**pdb*.
   3. During the simulation, write the iteration, temperature, and energy to a log file for later reference
3. Perform energy minimization on the best structure from stage 2. After each simulation, make sure to log the energy and the RMSD returned after relaxation. The relax function also saves a few intermediate structures in PDB files for visualization: if the best structure from your simulation is *best.pdb*, then *best_fast_relax.pdb* is the relaxed structure, and *best_fast_relax_aligned.pdb* is the relaxed structure superimposed into the same coordinate space as the native structure.

**Protein class**

We use the Protein class in **Protein.py** to store all information about our current conformation.

We store protein conformations using PyRosetta *Pose objects*. This makes our lives much easier by taking care of a lot of the complicated mechanics under the hood. You don’t need to worry about this too much, except realize that the *self.pose* attribute in **Protein.py** is a Pose object. **Do not try to change this to a different type of object, it is critical for the simulation to work.**

Since we work in torsion space for all fragment sampling, the only interactions you need to do with the pose objects is (1) access the ![LaTeX: \phi](https://canvas.stanford.edu/equation_images/%255Cphi)ϕ and ![LaTeX: \psi](https://canvas.stanford.edu/equation_images/%255Cpsi)ψ angles at a position, and (2) set the ![LaTeX: \phi](https://canvas.stanford.edu/equation_images/%255Cphi)ϕ and ![LaTeX: \psi](https://canvas.stanford.edu/equation_images/%255Cpsi)ψ angles at a position to new values. See the following tutorial for information on how to do so: [https://www.cse.huji.ac.il/~fora/81855/exercises/ex4.pdf](https://www.cse.huji.ac.il/~fora/81855/exercises/ex4.pdf)

Another important note about Pose objects: to copy a pose, you **cannot** simply define `pose_copy = pose`. This will simply create a new pointer to the same object in memory, so any changes made to `pose_copy` will also be made to `pose`. Instead, define a new pose and use `pose_copy.assign(pose)`, as shown in lines 34–35 of **Protein.py**.

 

**Fragment set class**

The pre-calculated fragment libraries for each input protein are given in the *.frag files, and the RMSDs of each fragment relative to the target structure are given in the *.rmsd files. In **FragmentSet.py**, implement a class that holds this fragment set information so that it can be accessed by your simulation.

In particular, in order to reduce the size of the sampling space, we do not sample from all fragments in the library at a given position. Instead, we *bias* the simulator such that it only samples from the top-N closest fragments in RMSD to the target structure at each position. To do this, the class must do the following:

- Parse the .frag and .rmsd file for a structure and store them in some data structure (think about what is most efficient for fast access). Since we do all fragment replacements in torsion space, we only need to store the (![LaTeX: \phi](https://canvas.stanford.edu/equation_images/%255Cphi)ϕ, ![LaTeX: \psi](https://canvas.stanford.edu/equation_images/%255Cpsi)ψ) angles for each fragment at each position.
- **Note:** the "position" we refer to here is the position in the protein you are folding. This is the position denoted in the header lines of the fragment files (see screenshot below). The position in the third column represents the position in the protein from which the fragment is derived, and can be ignored. Beware: the position listed after the literal "P" in the columns does not always correspond to the position in the header, so do not use this either!
- The fragment index is given in the last column. This corresponds to the fragment number in the .rmsd file.

![Screen Shot 2019-10-24 at 9.39.46 PM-1.png](https://canvas.stanford.edu/users/186210/files/5081057/preview?verifier=b1X4f19j8D1mcud3Mu4unBLNsqfZNetu0pTxPmY2)

- Implement the function `get_lowRMS_frags`, which returns the top N fragments as a list of lists. Each sub-list represents a fragment, and each element is a (![LaTeX: \phi](https://canvas.stanford.edu/equation_images/%255Cphi), ![LaTeX: \psi](https://canvas.stanford.edu/equation_images/%255Cpsi)) tuple representing the torsion at that fragment position in degrees. For example, a 3-mer fragment could be represented as the following: [(-60.892, 142.456), (-72.281, 128.933), (-132.337, -175.477)]. **The autograder will test this function, so do not change the output format.**

 

3. ### Main program specifications

Since the fragment assembly procedure is stochastic and is not guaranteed to find the minimum-energy conformation, we need to perform many simulations in order to find one that folds into the correct structure. In a real-life structure prediction setting, you would run on the order of 10,000 simulations to find the best folded structure. However, by using small, simple proteins and biased fragment sampling that considers only the top 3 fragments at each position, we do not have to do nearly as many to get a decent result.

Using the file **main.py**, write a program that runs **10** fragment-based structure prediction simulations using the procedure outlined above. Each simulation should use a different random seed, and thus should produce a different final structure.

After you have run all 10 simulations, find the simulation with the **lowest RMSD** using your log file *sim_log.txt*. You should be able to achieve a minimum RMSD of **less than 5 Å**. Finally, use PyMol to visualize *best_fast_relax_aligned.pdb* from the most successful simulation in the same window as *<target>.pdb*. The two structures should overlap closely. Take a screenshot of your best result for 1FW4 and 1UBQ.

**Inputs**

Your program **main.py** should take the following command-line inputs. The argparse library in Python makes this very straightforward. For code development, simply run the simulations using the default values, using helix.fasta and only one simulation. Run 10 simulations for each of 1FW4.fasta and 1UBQ.fasta; you will be asked to experiment with the other parameters for the quiz. Example command:

```
python main.py --fasta helix.fasta --logdir helix_log
```

| **Parameter** | **Description**                                      | **Default** |
| ------------- | ---------------------------------------------------- | ----------- |
| --fasta       | .fasta file containing sequence                      | None        |
| --logdir      | directory to save all log files                      | '.'         |
| --nsims       | number of simulations                                | 1           |
| --nfrags      | number of fragments to sample from at each iteration | 3           |
| --anneal_rate | temperature annealing parameter                      | 0.999       |

**Outputs**

**Log files:** You should save the following **tab-delimited** log files:

- For each input sequence, log a summary of all the simulations you have run, consisting of three columns, corresponding to the simulation number, the energy after relaxation, and the RMSD to target after relaxation:

```
   sim_number    energy    rmsd
```

- - **This is the most important log file for the autograder. It must be saved in the main log directory (specified by --logdir) as *simulation_summary.txt*.**
- For each simulation, log *at least* the following quantities throughout the simulation:

```
   iteration    temperature    energy
```

- - You can save other quantities if you to keep track of them, but you will need the above columns for the quiz.
  - You can log assembly and refinement stages to separate log files or the same one, it doesn't matter.

- We highly recommend creating a separate log directory for each simulation in order to store the PDB files and log files for each simulation. One possible directory structure is as follows (for the helix structure, with only 5 simulations). You don't need to save the initial or target structures, but it can be helpful in visualization.

  ​      ./helix_log

  ​         |__ /sim_01

  ​              |___ initial.pdb

  ​              |___ <target>.pdb

  ​              |___ best.pdb

  ​              |___ best_fast_relax.pdb

  ​              |___ best_fast_relax_aligned.pdb

  ​              |___ sim_01_log.txt

  ​         |__ /sim_02

  ​         |__ /sim_03

  ​         |__ /sim_04

  ​         |__ /sim_05

  ​         |__ simulation_summary.txt

**PDB files:** Each simulation should save at least *best.pdb*; after relaxation you should also get PDB files corresponding to relaxed and aligned structures.

**Helpful hints:**

- Make sure you are not repeating operations — this will result in inefficient code. This includes sampling positions and calculating the fragment candidates at that position.
- A "step" corresponds to an **accepted** move, not an **attempted** move. Make sure you only anneal temperature and move to the next step after accepting a move.
- Keep track of the fragments you have sampled at each position. If you sample every available fragment and do not accept any of them, you have likely already found the optimal fragment for the current configuration so you can move on to a new position.
- You can also speed things up by choosing smart data structures: lookups are much faster in hash-based objects like sets and dictionaries.
- It may be helpful to log certain quantities during the simulation to help you debug.
- When sampling positions, be very careful about which positions you are sampling. Also, note that **numpy.random.randint** samples from [low, high) *exclusive of the high value,* while **random.randint** samples from [low, high] *inclusive of the high value*. Make sure you are using the correct upper limit!
- The files given to you follow a consistent naming convention – take advantage of this and **do not hardcode the file names!** This will cause you to fail the autograder.
- To copy a Pose object**, do not** simply use `pose_copy = pose`. This is a shallow copy (i.e. will simply create a new pointer to the same object in memory, so any changes made to `pose_copy` will also be made to `pose`). Instead, define a new pose and use `pose_copy.assign(pose)`, as shown in lines 34–35 of **Protein.py**.

 

Submission Instructions

1. Take the quiz on Canvas. We recommend you take the quiz after having written the code and use your code, but you can take the quiz *without* having finished. It involves understanding the algorithm and how it works on different inputs. The quiz will take some time to complete, please budget accordingly. 
2. Upload all source code to [Gradescope (Links to an external site.)](https://www.gradescope.com/courses/49640)as a zipped folder. This must include **main.py, FragmentSampler.py**, **FragmentSet.py**, **Protein.py**, and **utils.py** in order for the autograder to work properly. Do not include irrelevant files.

References

- [Structure prediction with Rosetta (Links to an external site.)](https://www-sciencedirect-com.stanford.idm.oclc.org/science/article/pii/S0076687904830040)
- [Rosetta centroid energy function (Links to an external site.)](https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/centroid-score-terms)

 

### Appendix: How to use PyMol for visualization

Arguably the most interesting and satisfying part of any protein simulation is the visualization. Visually comparing the final folded structure of our protein to its native crystal structure allows us to better understand where our approach succeeds and where it struggles. It also provides a visual sanity check that a structure with a low RMSD and low energy actually does look similar to the real protein. Finally, visualizing how the structure changes over the course of the simulation helps us to understand how our simulation procedure works to find low-energy conformations.

One of the most popular software programs for visualizing protein structures is PyMol, which is free for students and available on all platforms (Linux, Mac, and Windows). You can download PyMol from [https://pymol.org. (Links to an external site.)](https://pymol.org/)

The PyMol interface consists of a large display screen, a console at the top (boxed in red), and an object panel on the right side (boxed in green).

![Screen Shot 2019-10-14 at 8.27.54 PM.png](https://canvas.stanford.edu/users/186210/files/5031416/preview?verifier=WkGsIBStnMkJN1tX9gJyW8d5hQ2T1ajphe7hvxhO)

The two main ways to interact with PyMol are by typing commands into the console and by using the buttons in the object panel. You will need to visualize some of your outputs for the project quiz. In addition to the tips listed below, there is a PyMol cheat sheet posted in the Project 3 files on Canvas.

Basic trackpad commands (for Mac, these may differ slightly on Windows):

- Rotate view: click and drag
- Zoom: two-finger pinch
- Pan: Option + click and drag

Some important console commands:

- To load a structure:

```
    load /path/to/structure.pdb
```

- To load a structure into a named object MyProtein (in PyMol syntax, commands are structured as 'command, object'):

```
    load /path/to/structure.pdb, MyProtein
```

- Loading multiple PDB files into one object results in **frames** (shown here as 1/100). If your files are labeled consistently, you can load multiple at once using the following command:

```
    for i in range(1,1001):cmd.load('pose_{}.pdb'.format(i), 'pose')
```

- These frames can be played through one-at-a-time or in a continuous loop using the play buttons in the bottom right of the screen. You can control framerate and other options under the 'Movie' menu. Try saving PDB files for the first 100-1000 iterations of a simulation, and then load them into PyMol and play through the frames to watch your protein fold!

![Screen Shot 2019-10-14 at 8.38.56 PM.png](https://canvas.stanford.edu/users/186210/files/5031422/preview?verifier=cfZSFpxh5STcYhcbshYS3yezOPEKyVM7EsCUbYjk)

- To delete an object:

```
    delete MyProtein
```

 

Some important notes about the object panel:

- Clicking the object name will toggle the object on/off
- The buttons ‘A’ (Action), ‘S’ (Show), ‘H’ (Hide), ‘L’ (Label), ‘C’ (Color) each reveal a menu that controls how the object is displayed. Play around with these menus a bit to get a sense for what each of them does.
- For the purposes of this assignment, we recommend visualizing the proteins in “cartoon” representation, which makes it easy to see the overall secondary structure of the protein (helices, sheets, and loops). Do this using ‘S’ > ‘as’ > ‘cartoon’.

![Screen Shot 2019-10-14 at 8.48.10 PM.png](https://canvas.stanford.edu/users/186210/files/5031423/preview?verifier=qqhGOyO2ouGuLDO1oQYQ3ZOHpcruMzMZVEkohGKJ)    ![Screen Shot 2019-10-14 at 8.53.41 PM.png](https://canvas.stanford.edu/users/186210/files/5031424/preview?verifier=fWk4pYLgKqOL2VYnWlcvvzpigV736KkffcsA3NTF)

- Color is controlled using the ‘C’ menu. It can be useful to color by secondary structure (‘C’ > ‘by ss’) to visualize high-level fold, as shown below.

![Screen Shot 2019-10-14 at 8.47.52 PM.png](https://canvas.stanford.edu/users/186210/files/5031427/preview?verifier=gFzGVVDGuahh0YK5JpcJWINuoDxbEHFtiuBl2glF)

- To select residues in the structure, you can simply click on the location you want to select. To view the corresponding position in the sequence or select residues by position, use the command

```
   set seq_view
```

- That command will result in the sequence appearing as shown below: you can select the amino acids by letter/position directly. This will be useful for the quiz.

![Screen Shot 2019-10-22 at 10.15.35 PM.png](https://canvas.stanford.edu/users/186210/files/5070801/preview?verifier=vWczhT8aWz8jxboad1kYtWjZMe7eaTOwq0jsnfG2)