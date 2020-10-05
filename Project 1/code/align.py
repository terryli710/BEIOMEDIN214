"""
This file provides the skeleton code for align.py.

Locations with "FILL IN" in comments are where you need to add code. Implementational notes and hints
are incorporated throughout. DO NOT modify the static types that are pre-defined for each function unless
explicitly mentioned. Changing expected behavior of functions will likely result in loss of points by the
autograder.

Usage: python align.py input_file output_file
"""

import sys
from typing import Set, Tuple, \
    List  # NOTE: You may need to "pip install typing" locally if this import gives you errors
import numpy as np


#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a: float, b: float) -> bool:
    """
    Checks if two floating point numbers are equivalent.
    :param a: first number
    :param b: second number
    :return: True if a and be are within epsilon = 10^-6
    """
    epsilon = 10 ** (-6)
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
        self.data = {}
        pass

    def set_score(self, a: str, b: str, score: float) -> None:
        """
        Updates or adds a score for a specified match

        :param a: the character from sequence A
        :param b: the character from sequence B
        :param score: the score to set the match M(a,b)
        """
        self.data[(a, b)] = score
        pass

    def get_score(self, a: str, b: str) -> float:
        """
        Returns the score for a particular match, where a is the
        character from sequence A and b is from sequence B.

        :param a: the character from sequence A
        :param b: the character from sequence B
        :return: the score of that match, M(a,b)
        """
        try:
            return self.data[(a, b)]
        except KeyError:
            print("Key error for {}->{}".format(a, b))
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

        class ScoreEntry(object):
            '''
            To store entries of score matrix.

            score: score of this entry
            pointer: (set of tuples) record the path for traceback
            '''

            def __init__(self):
                self.score = 0
                self.pointer = set()
        # set up matrix with instances of scoreEntry
        self.score_matrix = np.array([[ScoreEntry() for j in range(self.ncol)] for i in range(self.nrow)])
        # you need to figure out a way to represent this and how to initialize
        # Hint: it may be helpful to create a ScoreEntry class so that each entry in the
        # ScoreMatrix is an instance of this ScoreEntry class. This will help you keep track of all the
        # information that must be stored in each entry of the score matrix

    def __checkIndex__(self, row, col):
        if 0 <= row < self.nrow and 0 <= col < self.ncol:
            return True
        else:
            print("IndexError:")
            print("Score matrix dimension [{},{}], but requested position [{},{}]".format(self.nrow, self.ncol, row, col))
            return False

    def get_score(self, row: int, col: int) -> float:
        """
        Return the current score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :return: current score at position (row, col)
        """
        if self.__checkIndex__(row, col):
            entry = self.score_matrix[row, col]
            return entry.score
        else:
            return 0
        pass

    def set_score(self, row: int, col: int, score: float) -> None:
        """
        Set the score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :param score: score to set at position (row, col)
        """
        if self.__checkIndex__(row, col):
            self.score_matrix[row, col].score = round(score, 3)
        pass

    def get_pointers(self, row: int, col: int) -> Set[Tuple[int, int, str]]:
        """
        Return the indices of the entries being pointed to. Remember, these are essential for the final traceback.

        :param row: row index
        :param col: column index
        :return: a set of indices (represented as tuples) corresponding to other entries being pointed to for traceback
        ex: {(1,0), (1,1)}
        """

        if self.__checkIndex__(row, col):
            entry = self.score_matrix[row, col]
            return entry.pointer
        pass

    def set_pointers(self, row: int, col: int, pointers: Set[Tuple[int, int, str]]) -> None:
        """
        Add pointers to each entry in the score matrix.

        :param row: row index
        :param col: column index
        :param pointers: set of pointers to add to your entry
        """

        if self.__checkIndex__(row, col):
            self.score_matrix[row, col].pointer.update(pointers)
        pass

    def __formatTable__(self, data: np.array) -> str:
        """
        This function formats the table to make it nicely looked
        :param data: data to be formatted
        :return: a tring for print
        """
        # convert to string
        data = data.astype(str)
        # padding
        arrayLen = np.vectorize(len)
        col_width = np.max(arrayLen(data)) + 2
        output = []
        for row in data:
            output.append("".join(entry.ljust(col_width) for entry in row))
        return "\n".join(output)

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
        data = np.zeros((self.nrow, self.ncol))
        for j in range(self.ncol):
            for i in range(self.nrow):
                score = self.score_matrix[i, j].score
                data[i, j] = str(score)
        return self.__formatTable__(data)
        pass

    def print_pointers(self) -> str:
        """
        Returns a nicely formatted string containing the pointers in the score matrix. This function is OPTIONAL
        (i.e. will not be checked by autograder) but will be extremely helpful for debugging!
        """
        data = np.zeros((self.nrow, self.ncol)).astype(str)
        for j in range(self.ncol):
            for i in range(self.nrow):
                pointer = self.score_matrix[i, j].pointer
                if pointer:
                    data[i, j] = str(pointer)
                else:
                    data[i, j] = '{}'

        return self.__formatTable__(data)
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
        # read file and dividing values
        try:
            with open(input_file) as f:
                content = f.read()
        except FileNotFoundError:
            print("File name \"{}\" invalid.".format(input_file))
            return
        lines = [line.split() for line in content.split('\n')]
        self.seq_a = lines[0][0]
        self.seq_b = lines[1][0]
        self.global_alignment = not int(lines[2][0])
        self.dx = float(lines[3][0])
        self.ex = float(lines[3][1])
        self.dy = float(lines[3][2])
        self.ey = float(lines[3][3])
        self.len_alphabet_a = int(lines[4][0])
        self.alphabet_a = lines[5][0]
        self.len_alphabet_b = int(lines[6][0])
        self.alphabet_b = lines[7][0]
        for line in lines[8:]:
            # except for empty line
            if line == [''] or line == []: continue
            # set score iteratively
            self.match_matrix.set_score(line[2], line[3], float(line[4]))
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
        # loading parameters
        self.align_params.load_params_from_file(input_file)
        # seq_a as rows, seq_b as cols, set matrices to appropriate dims
        self.nrow, self.ncol = len(self.align_params.seq_a) + 1, len(self.align_params.seq_b) + 1
        # best score when A[i] matches with B[j]
        self.m_matrix = ScoreMatrix('M', self.nrow, self.ncol)
        # best score when A[i] matches with '_'
        self.ix_matrix = ScoreMatrix('Ix', self.nrow, self.ncol)
        # best score when B[j] matches with '_'
        self.iy_matrix = ScoreMatrix('Iy', self.nrow, self.ncol)

    def align(self, verbose=False):
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
        score, start_points = self.find_traceback_start()
        traces = []
        for start_point in start_points:
            traces += self.traceback(start_point)
        # print the results
        # for trace in traces: print(self.visualizePaths(trace))
        # store the alignments
        if verbose:
            for trace in traces: print(self.__printOutput__(trace), '\n')
        score_content = str(round(score, 1)) + '\n\n'
        self.alignment_result = score_content + '\n'.join(self.__printOutput__(trace) + '\n' for trace in traces)
        self.write_output()

    def populate_score_matrices(self) -> None:
        """
        Populate the score matrices based on the data in align_params. Should call update(i,j) for each entry
        in the score matrices.
        """
        # for initial entries, no end gap, set as 0
        for row in range(self.nrow):
            self.m_matrix.set_score(row, 0, 0)
            self.ix_matrix.set_score(row, 0, 0)
            self.iy_matrix.set_score(row, 0, 0)
        for col in range(self.ncol):
            self.m_matrix.set_score(0, col, 0)
            self.ix_matrix.set_score(0, col, 0)
            self.iy_matrix.set_score(0, col, 0)
        # for other entries, update accordingly
        for row in range(1, self.nrow):
            for col in range(1, self.ncol):
                self.update(row, col)
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

    def __getMax__(self, values: List[float], coord: List[Tuple[int, int, str]]):
        """
        Get the maximum score from provided values and their coordinates
        :param values: list of float value compare
        :param coord: coordinates as form of (row, col, matrix_name)
        :return: (max_value, {coordinates1, coordinates2})
        """
        # assert lengths are the same
        assert len(values) == len(coord), "value({}) and coordinates({}) length differ".format(len(values), len(coord))
        max_value = max(values)
        result = set()
        # iteratively find max values and positions
        for i in range(len(values)):
            if fuzzy_equals(max_value, values[i]): result.update([coord[i]])
        return tuple((max_value, result))

    def update_m(self, row: int, col: int) -> None:
        """
        Update matrix M.

        :param row: row index
        :param col: column index
        """

        # get s(i, j) from match matrix
        letter_a, letter_b = self.align_params.seq_a[row - 1], self.align_params.seq_b[col - 1]
        score_ij = self.align_params.match_matrix.get_score(letter_a, letter_b)
        # choices from which we choose the maximum
        values = [self.m_matrix.get_score(row - 1, col - 1) + score_ij,
                  self.ix_matrix.get_score(row - 1, col - 1) + score_ij,
                  self.iy_matrix.get_score(row - 1, col - 1) + score_ij]
        coords = [(row - 1, col - 1, self.m_matrix.name),
                  (row - 1, col - 1, self.ix_matrix.name),
                  (row - 1, col - 1, self.iy_matrix.name)]
        # get score and matrix
        max_score, pointers = self.__getMax__(values, coords)
        # set values
        # set negative values as 0 for local alignment
        if (not self.align_params.global_alignment) and (max_score < 0):
            self.m_matrix.set_score(row, col, 0)
        else:
            self.m_matrix.set_score(row, col, max_score)
        self.m_matrix.set_pointers(row, col, pointers)
        pass

    def update_ix(self, row: int, col: int) -> None:
        """
        Update matrix Ix.
        Iy records the best score if seq_a matches with '_'
        :param row: row index
        :param col: column index
        """

        # choices from which we choose the maximum
        pre_pos = (row - 1, col)
        values = [self.m_matrix.get_score(pre_pos[0], pre_pos[1]) - self.align_params.dy,
                  self.ix_matrix.get_score(pre_pos[0], pre_pos[1]) - self.align_params.ey]
        coords = [(pre_pos[0], pre_pos[1], self.m_matrix.name),
                  (pre_pos[0], pre_pos[1], self.ix_matrix.name)]
        # score and pointers
        max_score, pointers = self.__getMax__(values, coords)
        # set values
        # set negative values as 0 for local alignment
        if (not self.align_params.global_alignment) and (max_score < 0):
            self.ix_matrix.set_score(row, col, 0)
        else:
            self.ix_matrix.set_score(row, col, max_score)
        self.ix_matrix.set_pointers(row, col, pointers)
        pass

    def update_iy(self, row: int, col: int) -> None:
        """
        Update matrix Iy.
        Iy records the best score if seq_b matches with '_'
        :param row: row index
        :param col: column index
        """

        # choices from which we choose the maximum
        pre_pos = (row, col - 1)
        values = [self.m_matrix.get_score(pre_pos[0], pre_pos[1]) - self.align_params.dx,
                  self.iy_matrix.get_score(pre_pos[0], pre_pos[1]) - self.align_params.ex]
        coords = [(pre_pos[0], pre_pos[1], self.m_matrix.name),
                  (pre_pos[0], pre_pos[1], self.iy_matrix.name)]
        # score and pointer
        max_score, pointers = self.__getMax__(values, coords)
        # set values
        # set negative values as 0 for local alignment
        if (not self.align_params.global_alignment) and (max_score < 0):
            self.iy_matrix.set_score(row, col, 0)
        else:
            self.iy_matrix.set_score(row, col, max_score)
        self.iy_matrix.set_pointers(row, col, pointers)
        pass

    def find_traceback_start(self) -> Tuple[float, Set[Tuple[int, int, str]]]:
        """
        Find the location(s) to start the traceback and the corresponding best score.
        NOTE: Think carefully about how to set this up for local alignment.

        :return: The value of the best score and the location(s) to start the traceback to produce this score.
        The expected format is (max_val, max_loc), where max_val is the best score and max_loc is a set containing
        the positions to start the traceback that produce the best score, where each position is represented by a tuple
        [ex. (5.5, {(1,2), (3,4)}) ].
        """

        # if global alignment (no end gap penalty)
        if self.align_params.global_alignment:
            # search in the last row/col in matrices
            # last row + last col
            values = [self.m_matrix.get_score(self.nrow - 1, j) for j in range(self.ncol)] + \
                     [self.ix_matrix.get_score(self.nrow - 1, j) for j in range(self.ncol)] + \
                     [self.iy_matrix.get_score(self.nrow - 1, j) for j in range(self.ncol)] + \
                     [self.m_matrix.get_score(i, self.ncol - 1) for i in range(self.nrow)] + \
                     [self.ix_matrix.get_score(i, self.ncol - 1) for i in range(self.nrow)] + \
                     [self.iy_matrix.get_score(i, self.ncol - 1) for i in range(self.nrow)]
            coords = [(self.nrow - 1, j, self.m_matrix.name) for j in range(self.ncol)] + \
                     [(self.nrow - 1, j, self.ix_matrix.name) for j in range(self.ncol)] + \
                     [(self.nrow - 1, j, self.iy_matrix.name) for j in range(self.ncol)] + \
                     [(i, self.ncol - 1, self.m_matrix.name) for i in range(self.nrow)] + \
                     [(i, self.ncol - 1, self.ix_matrix.name) for i in range(self.nrow)] + \
                     [(i, self.ncol - 1, self.iy_matrix.name) for i in range(self.nrow)]
            result = self.__getMax__(values, coords)
        # if local alignment
        else:
            # search in the whole matrices
            values = []
            coords = []
            for matrix in [self.m_matrix, self.ix_matrix, self.iy_matrix]:
                for col in range(self.ncol):
                    for row in range(self.nrow):
                        values.append(matrix.get_score(row, col))
                        coords.append((row, col, matrix.name))
            result = self.__getMax__(values, coords)
        return result

    def __findMatrix__(self, name: str) -> ScoreMatrix:
        '''
        Map from matrix name to matrix instance.
        :param name: name of matrix
        :return: matrix instance
        '''
        matrix_dict = {matrix.name: matrix for matrix in [self.m_matrix, self.ix_matrix, self.iy_matrix]}
        return matrix_dict[name]

    def __tracebackRecursive__(self, position: Tuple[int, int, str],
                               traces: List[List[Tuple[int, int, str]]]) -> List[List[Tuple[int, int, str]]]:
        '''
        Recursively traceback, records all possible paths
        :param position: current position
        :param traces: recorded paths tracked
        :return: all recorded paths
        '''

        row, col, name = position
        # if end of matrix, return
        if row == 0 or col == 0:
            return traces
        # if not global and score is zero, return without the current one
        elif (not self.align_params.global_alignment) and fuzzy_equals(self.__findMatrix__(name).get_score(row,col),0):
            return [trace[:-1] for trace in traces]
        else:
            # find pointers
            matrix = self.__findMatrix__(name)
            pointers = matrix.get_pointers(row, col)
            results = []
            for pointer in pointers:
                # add new point in traces
                new_traces = [trace + [pointer] for trace in traces]
                # Run recursive operation!
                pointer_results = self.__tracebackRecursive__(pointer, new_traces)
                results += pointer_results
            return results

    def visualizePaths(self, path: List[Tuple[int, int, str]]) -> str:
        '''
        Visualize a path of traceback
        :param path: path input
        '''

        return '->'.join([self.__visualizePosition__(pos) for pos in path])

    def __visualizePosition__(self, position: Tuple[int, int, str]) -> str:
        '''
        Return a string to print a position
        :param position: position to be printed
        :return: string
        '''

        return '{}({},{})'.format(position[2], position[0], position[1])

    def __deleteDuplicatePaths__(self, traces: List[List[Tuple[int, int, str]]]) -> List[List[Tuple[int, int, str]]]:
        '''
        Drop duplicate paths (regardless of the last border position's difference)
        :param traces: list of paths
        :return: list of paths containing unique paths
        '''
        traces_wo_tail = [tuple(trace[:-1]) for trace in traces]
        trace_set = set()
        unique_index = []
        for i in range(len(traces_wo_tail)):
            trace = traces_wo_tail[i]
            if trace not in trace_set:
                unique_index.append(i)
                trace_set.update([trace])
        return [traces[i] for i in unique_index]

    def traceback(self, start_point: Tuple[int, int, str]) -> List[List[Tuple[int, int, str]]]:
        """
        Perform a traceback.

        NOTE: There is no static typing on the method, as you have freedom to choose how you'd like to implement
        this method, including which arguments to include and what to return.

        HINT: It is extremely helpful for debugging to include a way to print the traceback path.
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)
        """
        traces = self.__tracebackRecursive__(position=start_point, traces=[[start_point]])
        unique_traces = self.__deleteDuplicatePaths__(traces)
        return unique_traces

    def __omitEndGap__(self, alignment_a: str, alignment_b: str) -> Tuple[str, str]:
        return self.__recursiveDelete__(alignment_a, alignment_b)

    def __recursiveDelete__(self, str_a: str, str_b: str, pattern='_') -> Tuple[str, str]:
        if str_a.startswith(pattern) or str_b.startswith(pattern):
            return self.__recursiveDelete__(str_a[1:], str_b[1:])
        elif str_a.endswith(pattern) or str_b.endswith(pattern):
            return self.__recursiveDelete__(str_a[:-1], str_b[:-1])
        else: return str_a, str_b

    def __printOutput__(self, trace: List[Tuple[int, int, str]]) -> str:
        align_a, align_b = '', ''
        for element in trace:
            row, col, name = element
            if row == 0 or col == 0: break
            # if m_matrix, matched
            if name == self.m_matrix.name:
                align_a = self.align_params.seq_a[row - 1] + align_a
                align_b = self.align_params.seq_b[col - 1] + align_b
            # if ix_matrix, seq_a is matching with gap '_'
            elif name == self.ix_matrix.name:
                align_a = self.align_params.seq_a[row - 1] + align_a
                align_b = '_' + align_b
            # if iy_matrix, seq_b is matching with gap '_'
            elif name == self.iy_matrix.name:
                align_a = '_' + align_a
                align_b = self.align_params.seq_b[col - 1] + align_b
        return '\n'.join([align_a, align_b])

    def write_output(self) -> None:
        """
        Write the output of an alignment to the output file.
        """
        # if output_file is not specified, then save nothing
        import os
        if (self.output_file=='') or (not os.path.isdir(os.path.dirname(self.output_file))):
            print('Didn\'t save file')
            return
        with open(self.output_file, "w") as f:
            f.write(self.alignment_result)
        print('Saved to', self.output_file)
        return


# DO NOT MODIFY THE CODE BELOW!
def main():
    """
    Run the align function from command line, passing the input and output paths as arguments.
    """

    # Check that arguments are in present in the command line as expected
    if (len(sys.argv) != 3):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__ == "__main__":
    main()
