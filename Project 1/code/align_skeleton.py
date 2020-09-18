"""
This file provides the skeleton code for align.py.

Locations with "FILL IN" in comments are where you need to add code. Implementational notes and hints
are incorporated throughout. DO NOT modify the static types that are pre-defined for each function unless
explicitly mentioned. Changing expected behavior of functions will likely result in loss of points by the
autograder.

Usage: python align.py input_file output_file
"""


import sys
from typing import Set, Tuple  # NOTE: You may need to "pip install typing" locally if this import gives you errors


#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a: float, b: float) -> bool:
    """
    Checks if two floating point numbers are equivalent.
    :param a: first number
    :param b: second number
    :return: True if a and be are within epsilon = 10^-6
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

#### ------- CLASSES ------- ####
class MatchMatrix(object):
    """
    A class representation of the match matrix, S. Stores the scores of matches between characters.
    """
    def __init__(self):
        """
        Initialize MatchMatrix class.

        NOTE: There is a bit of freedom in how to implement this.
        """
        ### FILL IN ###
        pass

    def set_score(self, a: str, b: str, score: float) -> None:
        """
        Updates or adds a score for a specified match

        :param a: the character from sequence A
        :param b: the character from sequence B
        :param score: the score to set the match M(a,b)
        """
        ### FILL IN ###
        pass

    def get_score(self, a: str, b: str) -> float:
        """
        Returns the score for a particular match, where a is the
        character from sequence A and b is from sequence B.

        :param a: the character from sequence A
        :param b: the character from sequence B
        :return: the score of that match, M(a,b)
        """
        ### FILL IN ###
        pass


class ScoreMatrix(object):
    """
    A class representation of the score matrices (M, Ix, Iy), which will be dynamically updated.
    The score matrix consists of a 2-D array of ScoreEntries that are updated during alignment
    and used to output the maximum alignment.
    """

    def __init__(self, name: str, nrow: int, ncol: int) -> None:
        """
        Initialize ScoreMatrix class.

        :param name: identifier for the score matrix, should be in {Ix, Iy, M}
        :param nrow: number of rows for ScoreMatrix
        :param ncol: number of columns for ScoreMatrix
        """
        self.name = name
        self.nrow = nrow
        self.ncol = ncol

        self.score_matrix  = ...# FILL IN
        # you need to figure out a way to represent this and how to initialize
        # Hint: it may be helpful to create a ScoreEntry class so that each entry in the
        # ScoreMatrix is an instance of this ScoreEntry class. This will help you keep track of all the
        # information that must be stored in each entry of the score matrix

    def get_score(self, row: int, col: int) -> int:
        """
        Return the current score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :return: current score at position (row, col)
        """

        ### FILL IN ###
        pass

    def set_score(self, row: int, col: int, score: int) -> None:
        """
        Set the score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :param score: score to set at position (row, col)
        """

        ### FILL IN ###
        pass

    def get_pointers(self, row: int, col: int) -> Set[Tuple[int, int]]:
        """
        Return the indices of the entries being pointed to. Remember, these are essential for the final traceback.

        :param row: row index
        :param col: column index
        :return: a set of indices (represented as tuples) corresponding to other entries being pointed to for traceback
        ex: {(1,0), (1,1)}
        """

        ### FILL IN ###
        pass

    def set_pointers(self, row: int, col: int, pointers: Set[Tuple[int, int]]) -> None:
        """
        Add pointers to each entry in the score matrix.

        :param row: row index
        :param col: column index
        :param pointers: set of pointers to add to your entry
        """

        ### FILL IN ###
        pass

    def print_scores(self) -> str:
        """
        Returns a nicely formatted string containing the scores in the score matrix. This function is OPTIONAL
        (i.e. will not be checked by autograder) but will be extremely helpful for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        :return: a string representation of the scores in the score matrix
        """

        ### FILL IN ###
        pass

    def print_pointers(self) -> str:
        """
        Returns a nicely formatted string containing the pointers in the score matrix. This function is OPTIONAL
        (i.e. will not be checked by autograder) but will be extremely helpful for debugging!
        """

        ### FILL IN ###
        pass


class AlignmentParameters(object):
    """
    A class to hold the alignment parameters.
    """

    def __init__(self) -> None:
        """
        Initialize AlignmentParameters object with default parameters.
        """

        # The definition for all of these class variables are documented in the annotated version of the input file
        # on the P1 page on Canvas.
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False 
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.len_alphabet_a = 0
        self.alphabet_a = ""
        self.len_alphabet_b = 0
        self.alphabet_b = ""
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file: str) -> None:
        """
        Read the parameters from an input file and update the alignment parameters accordingly

        :param input_file: path to the alignment input file (whose structure is defined on the project page)
        """

        ### FILL IN ###
        pass


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by calling "align()"

    NOTE ON IMPLEMENTATION: You might find it helpful for code structure/organization to create additional functions
    that handle subtasks during the alignment. This is totally acceptable. However, you MUST implement the following
    functions with the expected behavior as outlined in the docstrings, which we will check in order to give
    at least partial credit for your P1 implementations:

    - populate_score_matrices
    - update_m, update_ix, and update_iy
    - find_traceback_start

    NOTE 2: Don't forget about that fuzzy_equals function at the top.

    """

    def __init__(self, input_file: str, output_file: str) -> None:
        """
        Initialize Align object.

        :param input_file: alignment input file path
        :param output_file: file path to write the output alignment
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters()

        ### FILL IN ###
        # NOTE: be careful about how you initialize these!
        self.m_matrix = ...
        self.ix_matrix = ...
        self.iy_matrix = ...

    def align(self):
        """
        Main method for running the alignment.

        Note: there is no static typing on the method, as you can choose to return arbitrary
        intermediates/output if it's helpful for debugging. The essential minimal functionality that this
        method must have is to write the resulting alignments to the output file in the format specified in the
        project page
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # perform a traceback and write the output to an output file
        ### FILL IN ###

    def populate_score_matrices(self) -> None:
        """
        Populate the score matrices based on the data in align_params. Should call update(i,j) for each entry
        in the score matrices.
        """

        ### FILL IN ###
        pass

    def update(self, row: int, col: int) -> None:
        """
        Update all matrices at a given row and column index.

        :param row: index of row to update
        :param col: index of column to update
        """
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row: int, col: int) -> None:
        """
        Update matrix M.

        :param row: row index
        :param col: column index
        """

        ### FILL IN ###
        pass

    def update_ix(self, row: int, col: int) -> None:
        """
        Update matrix Ix.

        :param row: row index
        :param col: column index
        """

        ### FILL IN ###
        pass

    def update_iy(self, row: int, col: int) -> None:
        """
        Update matrix Iy.

        :param row: row index
        :param col: column index
        """

        ### FILL IN ###
        pass

    def find_traceback_start(self) -> Tuple[float, Set[Tuple[int, int]]]:
        """
        Find the location(s) to start the traceback and the corresponding best score.
        NOTE: Think carefully about how to set this up for local alignment.

        :return: The value of the best score and the location(s) to start the traceback to produce this score.
        The expected format is (max_val, max_loc), where max_val is the best score and max_loc is a set containing
        the positions to start the traceback that produce the best score, where each position is represented by a tuple
        [ex. (5.5, {(1,2), (3,4)}) ].
        """

        ### FILL IN ###
        pass

    def traceback(self): ### FILL IN additional arguments ###
        """
        Perform a traceback.

        NOTE: There is no static typing on the method, as you have freedom to choose how you'd like to implement
        this method, including which arguments to include and what to return.

        HINT: It is extremely helpful for debugging to include a way to print the traceback path.
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)
        """

        ### FILL IN ###
        pass

    def write_output(self) -> None:
        """
        Write the output of an alignment to the output file.
        """

        ### FILL IN ###
        pass


# DO NOT MODIFY THE CODE BELOW!
def main():
    """
    Run the align function from command line, passing the input and output paths as arguments.
    """

    # Check that arguments are in present in the command line as expected
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__=="__main__":
    main()
