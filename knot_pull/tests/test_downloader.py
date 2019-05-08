from unittest import TestCase

import knot_pull.downloader as dwn

class TestConnection(TestCase):
    def test_urlopen_works(self):
        self.assertRaises(dwn.URLError, dwn._urlopen,"http://www.stolbeznog.tv")
        self.assertIsNotNone(dwn._urlopen, "https://www.google.com")


class TestPDBreading(TestCase):

    def setUp(self):
        self.pdbid = '1uak'

    def test_access(self):
        self.assertTrue(dwn.check_if_in_PDB(self.pdbid,"A"))
        self.assertFalse(dwn.check_if_in_PDB(self.pdbid,"C"))
        self.assertEqual(dwn.get_ligands(self.pdbid,"A"),["SAM"])
        self.assertEqual(len(list(dwn.get_particular_chain(self.pdbid,"A"))),243)
        self.assertEqual(list(dwn.get_particular_chain(self.pdbid,"A")),dwn.get_all_chains(self.pdbid))
        self.assertEqual(dwn.get_chain_list(self.pdbid),["A"])
