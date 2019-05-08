from unittest import TestCase

import knot_pull.downloader as dwn

class TestConnection(TestCase):
    def test_urlopen_works(self):
        self.assertRaises(dwn.URLError, dwn._urlopen,"http://www.stolbeznog.tv")
        self.assertIsNotNone(dwn._urlopen, "https://www.google.com")


