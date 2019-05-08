from unittest import TestCase

import knot_pull

class TestDowker(TestCase):
    def test_dowker_is_string(self):
        s = "nothing for now"
        self.assertTrue(isinstance(s, basestring))