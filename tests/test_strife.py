import pyximport; pyximport.install()
import unittest
import strife.strife.strife
import scipy
from scipy.weave import inline

class TestStrife(unittest.TestCase):
    def setUp(self):
        self.a_game = strife.strife.strife.Strife()
        board_a = scipy.array([[ True,  True, False,  True,  True],
                               [False,  True, False, False, False],
                               [False, False,  True, False,  True],
                               [ True,  True, False,  True, False],
                               [ True,  True,  True,  True,  True]], dtype=bool)
        board_b = scipy.array([[False,  True,  True, False,  True],
                               [ True, False, False,  True, False],
                               [False, False,  True,  True,  True],
                               [ True,  True, False, False,  True],
                               [False, False, False, False,  True]], dtype=bool)
        board_c = scipy.array([[ True, False, False,  True, False],
                               [False, False, False, False,  True],
                               [ True, False, False,  True,  True],
                               [ True,  True, False, False,  True],
                               [ True,  True,  True,  True, False]], dtype=bool)
        self.a_game.board = scipy.array([board_a, board_b, board_c]).transpose(1,2,0)

class TestBoardShape(TestStrife):
    def test_board_shape(self):
        self.assertEqual(self.a_game.board.shape, (5,5,3))

class TestTruth(unittest.TestCase):#Strife):
    def test_inline_truth(self):
        self.assertEqual(scipy.array(True), True)
        inline_code = r"""
            signed long int q = {0};
            printf("1: %d; ", q);
            printf("%d == 1: ", T1(0));
            printf("%d\n", T1(0) == q);
            return_val = (T1(0) == q);"""
        t = scipy.bool_(False)
        self.assertEqual(inline(inline_code.format(0), ['t']), 1)
        self.assertEqual(inline(inline_code.format(1), ['t']), 0)
        t = scipy.bool_(True)
        self.assertEqual(inline(inline_code.format(0), ['t']), 0)
        self.assertEqual(inline(inline_code.format(1), ['t']), 1)

class TestCount(unittest.TestCase):
    def test_count(self):
        a_game = strife.strife.strife.Strife()
        a_game.board = scipy.array([[[1, 0, 0, 1],
                                     [1, 1, 0, 1],
                                     [0, 1, 1, 0],
                                     [1, 1, 1, 1]]], dtype='int64').transpose(1,2,0)
        '''
        with rows = [1,3), cols = [1,2)

        at center = (1,1) and radius = 1
        1 0 0
        1 1 0 == 5
        0 1 1

        center = (2,1)
        1 1 0
        0 1 1 == 7
        1 1 1
        '''
        sum_board = a_game.count_neighbors(rows   = scipy.array([1, 3]),
                                           cols   = scipy.array([1, 2]),
                                           gene   = 0,
                                           allele = 1,
                                           radius = 1)
        print("sumboard, array", sum_board, scipy.array([[5], [7]]))
        self.assertEqual((scipy.array([[5], [7]]) == sum_board).all(), True)

class TestDiffuse2by2(TestStrife):
    def test_diffuse(self):
        board_a = scipy.arange(2*2).reshape(2,2)
        board_b = scipy.arange(2*2, 2*2*2).reshape(2,2)
        board_c = scipy.arange(2*2*2, 3*2*2).reshape(2,2)
        board_0 = scipy.array([board_a, board_b, board_c]).transpose(1,2,0)
        self.a_game.board = board_0.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((0, 0)))
        board_90 = self.a_game.board.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((0, 0)))
        board_180 = self.a_game.board.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((0, 0)))
        board_270 = self.a_game.board.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((0, 0)))
        board_360 = self.a_game.board.copy()
        priori_board_90_a = scipy.array([[1,3],
                                         [0,2]])
        priori_board_90_b = scipy.array([[5,7],
                                         [4,6]])
        priori_board_90_c = scipy.array([[9,11],
                                         [8,10]])
        priori_board_90 = scipy.array([priori_board_90_a,
                                       priori_board_90_b,
                                       priori_board_90_c]).transpose(1,2,0)
        priori_board_180_a = scipy.array([[3,2],
                                         [1,0]])
        priori_board_180_b = scipy.array([[7,6],
                                         [5,4]])
        priori_board_180_c = scipy.array([[11,10],
                                         [9,8]])
        priori_board_180 = scipy.array([priori_board_180_a,
                                        priori_board_180_b,
                                        priori_board_180_c]).transpose(1,2,0)
        priori_board_270_a = scipy.array([[2,0],
                                          [3,1]])
        priori_board_270_b = scipy.array([[6,4],
                                          [7,5]])
        priori_board_270_c = scipy.array([[10,8],
                                          [11,9]])
        priori_board_270 = scipy.array([priori_board_270_a,
                                        priori_board_270_b,
                                        priori_board_270_c]).transpose(1,2,0)
        priori_board_360_a = scipy.array([[0,1],
                                          [2,3]])
        priori_board_360_b = scipy.array([[4,5],
                                          [6,7]])
        priori_board_360_c = scipy.array([[8,9],
                                          [10,11]])
        priori_board_360 = scipy.array([priori_board_360_a,
                                        priori_board_360_b,
                                        priori_board_360_c]).transpose(1,2,0)
        self.assertEqual((board_90  == priori_board_90 ).all(), True)
        self.assertEqual((board_180 == priori_board_180).all(), True)
        self.assertEqual((board_270 == priori_board_270).all(), True)
        self.assertEqual((board_360 == priori_board_360).all(), True)

class TestDiffuse2by2Modulus(TestStrife):
    def test_diffuse(self):
        board_a = scipy.arange(2*2).reshape(2,2)
        board_b = scipy.arange(2*2, 2*2*2).reshape(2,2)
        board_c = scipy.arange(2*2*2, 3*2*2).reshape(2,2)
        board_0 = scipy.array([board_a, board_b, board_c]).transpose(1,2,0)
        self.a_game.board = board_0.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((1, 1)))
        board_90 = self.a_game.board.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((1, 1)))
        board_180 = self.a_game.board.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((1, 1)))
        board_270 = self.a_game.board.copy()
        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((1, 1)))
        board_360 = self.a_game.board.copy()
        priori_board_90_a = scipy.array([[1,3],
                                         [0,2]])
        priori_board_90_b = scipy.array([[5,7],
                                         [4,6]])
        priori_board_90_c = scipy.array([[9,11],
                                         [8,10]])
        priori_board_90 = scipy.array([priori_board_90_a,
                                       priori_board_90_b,
                                       priori_board_90_c]).transpose(1,2,0)
        priori_board_180_a = scipy.array([[3,2],
                                         [1,0]])
        priori_board_180_b = scipy.array([[7,6],
                                         [5,4]])
        priori_board_180_c = scipy.array([[11,10],
                                         [9,8]])
        priori_board_180 = scipy.array([priori_board_180_a,
                                        priori_board_180_b,
                                        priori_board_180_c]).transpose(1,2,0)
        priori_board_270_a = scipy.array([[2,0],
                                          [3,1]])
        priori_board_270_b = scipy.array([[6,4],
                                          [7,5]])
        priori_board_270_c = scipy.array([[10,8],
                                          [11,9]])
        priori_board_270 = scipy.array([priori_board_270_a,
                                        priori_board_270_b,
                                        priori_board_270_c]).transpose(1,2,0)
        priori_board_360_a = scipy.array([[0,1],
                                          [2,3]])
        priori_board_360_b = scipy.array([[4,5],
                                          [6,7]])
        priori_board_360_c = scipy.array([[8,9],
                                          [10,11]])
        priori_board_360 = scipy.array([priori_board_360_a,
                                        priori_board_360_b,
                                        priori_board_360_c]).transpose(1,2,0)
        self.assertEqual((board_90  == priori_board_90 ).all(), True)
        self.assertEqual((board_180 == priori_board_180).all(), True)
        self.assertEqual((board_270 == priori_board_270).all(), True)
        self.assertEqual((board_360 == priori_board_360).all(), True)

class TestDiffuse3by3Modulus(TestStrife):
    def test_diffuse(self):
        board_a = scipy.arange(4*4).reshape(4,4)
        board_b = scipy.arange(4*4, 4*4*2).reshape(4,4)
        board_c = scipy.arange(4*4*2, 4*4*3).reshape(4,4)

        board_0 = scipy.array([board_a, board_b, board_c]).transpose(1,2,0)
        self.a_game.board = board_0.copy()

        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((3, 3)))
        board_90 = self.a_game.board.copy()

        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((3, 3)))
        board_180 = self.a_game.board.copy()

        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((3, 3)))
        board_270 = self.a_game.board.copy()

        strife.strife.strife.rotquad90(self.a_game.board, 0, scipy.array((3, 3)))
        board_360 = self.a_game.board.copy()

        priori_board_90 = scipy.array([[[ 3,  1,  2, 15],
                                        [ 4,  5,  6,  7],
                                        [ 8,  9, 10, 11],
                                        [ 0, 13, 14, 12]],
                                       [[19, 17, 18, 31],
                                        [20, 21, 22, 23],
                                        [24, 25, 26, 27],
                                        [16, 29, 30, 28]],
                                       [[35, 33, 34, 47],
                                        [36, 37, 38, 39],
                                        [40, 41, 42, 43],
                                        [32, 45, 46, 44]]]).transpose(1,2,0)

        priori_board_180 = scipy.array([[[15,  1,  2, 12],
                                         [ 4,  5,  6,  7],
                                         [ 8,  9, 10, 11],
                                         [ 3, 13, 14,  0]],
                                        [[31, 17, 18, 28],
                                         [20, 21, 22, 23],
                                         [24, 25, 26, 27],
                                         [19, 29, 30, 16]],
                                        [[47, 33, 34, 44],
                                         [36, 37, 38, 39],
                                         [40, 41, 42, 43],
                                         [35, 45, 46, 32]]]).transpose(1,2,0)

        priori_board_270 = scipy.array([[[12,  1,  2,  0],
                                         [ 4,  5,  6,  7],
                                         [ 8,  9, 10, 11],
                                         [15, 13, 14,  3]],
                                        [[28, 17, 18, 16],
                                         [20, 21, 22, 23],
                                         [24, 25, 26, 27],
                                         [31, 29, 30, 19]],
                                        [[44, 33, 34, 32],
                                         [36, 37, 38, 39],
                                         [40, 41, 42, 43],
                                         [47, 45, 46, 35]]]).transpose(1,2,0)

        priori_board_360 = scipy.array([[[ 0,  1,  2,  3],
                                         [ 4,  5,  6,  7],
                                         [ 8,  9, 10, 11],
                                         [12, 13, 14, 15]],
                                        [[16, 17, 18, 19],
                                         [20, 21, 22, 23],
                                         [24, 25, 26, 27],
                                         [28, 29, 30, 31]],
                                        [[32, 33, 34, 35],
                                         [36, 37, 38, 39],
                                         [40, 41, 42, 43],
                                         [44, 45, 46, 47]]]).transpose(1,2,0)

        self.assertEqual((board_90  == priori_board_90 ).all(), True)
        self.assertEqual((board_180 == priori_board_180).all(), True)
        self.assertEqual((board_270 == priori_board_270).all(), True)
        self.assertEqual((board_360 == priori_board_360).all(), True)

'''suite = unittest.TestLoader().loadTestsFromTestCase(TestTrueEq1)
unittest.TextTestRunner(verbosity=2).run(suite)
'''

#if __name__ == '__main__':
#    unittest.main()
