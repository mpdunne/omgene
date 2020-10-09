from utils.misc import anyperm, del_indices, grab_lines

from  unittest import TestCase


class TestMisc(TestCase):

    def test_grab_lines(self):
        lines = grab_lines('echo "hello\nthese\nare\nsome\nlines"')
        self.assertListEqual(['hello', 'these', 'are', 'some', 'lines', ''], lines)

    def test_del_indices(self):
        input_list = [10, 11, 12, 13, 14, 15, 16, 17, 18]

        check1 = del_indices(input_list, [])
        self.assertEqual([10, 11, 12, 13, 14, 15, 16, 17, 18], check1)

        check2 = del_indices(input_list, [1, 3, 6])
        self.assertEqual([10, 12, 14, 15, 17, 18], check2)

        check3 = del_indices(input_list, [6, 1, 3])
        self.assertEqual([10, 12, 14, 15, 17, 18], check3)

        check4 = del_indices(input_list, [6, 1, 1, 1, 1])
        self.assertEqual([10, 12, 13, 14, 15, 17, 18], check4)

        check5 = del_indices(input_list, [6, 100, 1])
        self.assertEqual([10, 12, 13, 14, 15, 17, 18], check5)

    def test_anyperm(self):
        check1 = anyperm([1, 2, 3], [[1, 3, 2]])
        self.assertTrue(check1)

        check2 = anyperm([1, 2], [[1, 3, 2]])
        self.assertFalse(check2)

        check3 = anyperm([1, 2, 3], [[1, 3, 2], [1, 3, 5]])
        self.assertTrue(check3)

