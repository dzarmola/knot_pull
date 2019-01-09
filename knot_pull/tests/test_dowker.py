from unittest import TestCase

import knotPull

class TestDowker(TestCase):
    def dowker_is_string(self):
        s = "nothing for now"
        self.assertTrue(isinstance(s, basestring))