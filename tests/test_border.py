# tests/test_border.py
import unittest
from border.border import Border


class TestBorder(unittest.TestCase):
    def setUp(self):
        self.border = Border(30, 44, 68, 70, 100, 4, 10)

    def test_get_upper(self):
        self.assertEqual(self.border.get_upper(65), 4)

    def test_get_lower(self):
        self.assertAlmostEqual(self.border.get_lower(65), -5.74549562499999, places=1)

    def test_get_distance(self):
        self.assertAlmostEqual(self.border.get_distance(65), 9.74549562499999, places=1)


if __name__ == "__main__":
    unittest.main()
