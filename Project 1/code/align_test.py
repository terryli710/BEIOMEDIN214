"""
Unit tests provide a way of testing individual components of your program.

This will help you with debugging and making sure you don't break your code! Ensuring these tests pass will be
extremely important before submitting to the autograder, which will permit only a limited number of submissions.

Here, we provide some unit tests and a skeleton for a few others.
Note that you do not have to follow these exactly, they are designed to help you.

What other unit tests might you want to write?
 - Think about the traceback and writing the output file. 
 - Try to write at least one or two additional tests as you think of them.

To run:
  python align_test.py

Make sure align.py is located in the same directory, and the test_example.input file is present!
"""

import unittest
import os
from align import *
from align_quiz_functions import *

TEST_INPUT_FILE = "test_example.input"


class TestAlignmentClasses(unittest.TestCase):

    def test_match_matrix(self):
        """
        Tests match matrix object
        """
        match_matrix = MatchMatrix()
        match_matrix.set_score("A", "C", 5)
        self.assertEqual(match_matrix.get_score("A", "C"), 5)


    def test_score_matrix_score(self):
        """
        Tests score matrix object score set + get methods
        """
        score_matrix = ScoreMatrix('test', 3, 4)
        score_matrix.set_score(1,2,10)
        self.assertEqual(score_matrix.get_score(1,2), 10)
        # this should be very similar to test match matrix
        return


    def test_score_matrix_pointers(self):
        """
        Tests score matrix object pointer set + get methods
        """
        score_matrix = ScoreMatrix('test', 4, 5)
        score_matrix.set_pointers(2, 3, {(0,1), (1,0)})
        self.assertEqual(score_matrix.get_pointers(2,3), {(0,1), (1,0)})
        return


    def test_score_matrix_print(self):
        score_matrix = ScoreMatrix('test', 4, 5)
        score_matrix.set_score(1, 2, 10)
        score_matrix.set_pointers(2, 3, {(0, 1), (1, 0)})
        scores = score_matrix.print_scores()
        pointers = score_matrix.print_pointers()
        print(scores)
        print(pointers)


    def test_param_loading(self):
        """
        Tests AlignmentParameters "load_params_from_file()" function
        """
        align_params = AlignmentParameters()
        align_params.load_params_from_file(TEST_INPUT_FILE)
        self.assertEqual(align_params.seq_a, "AATGC")
        self.assertEqual(align_params.seq_b, "AGGC")
        self.assertTrue(align_params.global_alignment)
        self.assertEqual(align_params.dx, 0.1)
        self.assertEqual(align_params.ex, 0.5)
        self.assertEqual(align_params.dy, 0.6)
        self.assertEqual(align_params.ey, 0.3)
        self.assertEqual(align_params.alphabet_a, "ATGC")
        self.assertEqual(align_params.alphabet_b, "ATGCX")
        self.assertEqual(align_params.len_alphabet_a, 4)
        self.assertEqual(align_params.len_alphabet_b, 5)

        # test that the match match is set up correctly
        #  if this fails, make sure you are loading the asymmetric matrix properly!
        match_mat = align_params.match_matrix
        self.assertEqual(match_mat.get_score("A", "X"), 0.3)
        self.assertEqual(match_mat.get_score("C", "G"), -0.3)
        self.assertEqual(match_mat.get_score("G", "C"), 0)


    def test_update_ix(self):
        """
        Test AlignmentAlgorithm's update Ix
        """

        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dy = 1
        align_params.ey = 0.5

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.m_matrix.set_score(2,2, 3)
        align.ix_matrix.set_score(2,2, 2.5)

        # run the method!
        align.update_ix(3, 2)

        score = align.ix_matrix.get_score(3,2)
        self.assertEqual(score, 2)

        pointers = align.ix_matrix.get_pointers(3,2)
        names = set([pointer[1] for pointer in pointers])
        self.assertTrue(set(("M", "Ix")).issubset(names))
        # note - in this example, it should point to M -AND- Ix
        # check by hand!


    def test_update_m(self):
        """
        Test AlignmentAlgorithm's update M
        """
        # configure alignment params
        align_params = AlignmentParameters()
        align_params.match_matrix.set_score("a","a", 1)
        align_params.seq_a = "aaaa"
        align_params.seq_b = "aaa"
        # create an alignment object
        align = Align("", "")
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2, 2, 3)
        align.ix_matrix.set_score(2, 2, 2.5)
        align.iy_matrix.set_score(2, 2, 5)


        # run the method!
        align.update_m(3, 3)

        score = align.m_matrix.get_score(3, 3)
        self.assertEqual(score, 6)

        pointers = align.m_matrix.get_pointers(3, 3)
        names = set([pointer[2] for pointer in pointers])
        self.assertTrue(set(("Iy")).issubset(names))
        return


    def test_update_iy(self):
        """
        Test AlignmentAlgorithm's update Iy
        """
        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dx = 1
        align_params.ex = 0.5

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2, 2, 3)
        align.iy_matrix.set_score(2, 2, 2.5)

        # run the method!
        align.update_iy(3, 2)

        score = align.iy_matrix.get_score(3, 2)
        self.assertEqual(score, 2)

        pointers = align.iy_matrix.get_pointers(3, 2)
        names = set([pointer[1] for pointer in pointers])
        self.assertTrue(set(("M", "Iy")).issubset(names))
        return


    def test_traceback_start(self):
        """
        Tests that the traceback finds the correct start
        Should test local and global alignment!
        """
        # Given some matrices, should find the correct start
        # configure alignment params
        align_params = AlignmentParameters()

        align = Align("", "")
        align.align_params = align_params

        # len(A) = len(B) = 2
        align.nrow, align.ncol = 3, 3
        align.m_matrix = ScoreMatrix("M", 3, 3)
        align.ix_matrix = ScoreMatrix("Ix", 3, 3)
        align.iy_matrix = ScoreMatrix("Iy", 3, 3)
        # set arbitrary scores
        align.m_matrix.set_score(1, 1, 5)
        align.ix_matrix.set_score(2, 1, 3)

        # if global
        align.align_params.global_alignment = True
        traceback_start_global = align.find_traceback_start()
        self.assertEqual(traceback_start_global, (3, {(2, 1, 'Ix')}))
        # if local
        align.align_params.global_alignment = False
        traceback_start_local = align.find_traceback_start()
        self.assertEqual(traceback_start_local, (5, {(1, 1, 'M')}))
        return


    def test_traceback(self):
        '''
        Test traceback function
        '''
        align = Align("test_example.input", "test_example.output")
        align.populate_score_matrices()
        traces = align.traceback(start_point=(5, 4, 'M'))
        for trace in traces:
            print(align.visualizePaths(trace))
        return

    def __testPaths__(self):
        # get project path
        file_path = os.path.dirname(os.path.abspath(__file__))
        project_path = os.path.abspath(os.path.join(file_path, os.path.pardir))
        # get example paths
        example_dir = os.path.join(project_path, 'examples')
        output_dir = os.path.join(project_path, 'test')
        return example_dir, output_dir

    def test_align(self):
        example_dir, output_dir = self.__testPaths__()
        input = os.path.join(example_dir, 'alignment_example4.input')
        output = os.path.join(example_dir, 'alignment_example4.output')
        input_path = os.path.join(example_dir, input)
        solution_path = os.path.join(example_dir, output)
        output_path = os.path.join(output_dir, output)
        align_test = Align(input_path, output_path)
        align_test.populate_score_matrices()
        align_test.align()
        output_align, solution_align = read_seq_pairs(output_path), read_seq_pairs(solution_path)
        self.assertTrue(compare_alignments(output_path, solution_path))
        return

    def __diffAlign__(self, align_a: List[List[str]], align_b: List[List[str]]):
        '''
        Find the part of alignments that only appear in one of the alignments
        :param align_a: alignment of a from function "read_seq_pairs"
        :param align_b: alignment of b from function "read_seq_pairs"
        :return: tuple of diff
        '''
        diff_a = [i for i in align_a if i not in align_b]
        diff_b = [i for i in align_b if i not in align_a]
        return (diff_a, diff_b)

    def test_examples(self):
        # get project path
        example_dir, output_dir = self.__testPaths__()
        examples_names = os.listdir(example_dir)
        example_input = [example for example in examples_names if example.endswith('.input')]
        example_output = [example.replace('.input', '.output') for example in example_input]
        test_success = True
        for input, output in zip(example_input, example_output):
            input_path = os.path.join(example_dir, input)
            solution_path = os.path.join(example_dir, output)
            output_path = os.path.join(output_dir, output)
            align_test = Align(input_path, output_path)
            align_test.populate_score_matrices()
            align_test.align()
            output_align, solution_align = read_seq_pairs(output_path), read_seq_pairs(solution_path)
            if not compare_alignments(output_path, solution_path):
                # print(self.__diffAlign__(output_align, solution_align))
                test_success = False
        self.assertTrue(test_success)

if __name__ == '__main__':
    unittest.main()