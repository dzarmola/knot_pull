from unittest import TestCase

import knot_pull.downloader as dwn

class TestConnection(TestCase):
    def test_urlopen_works(self):
        self.assertRaises(dwn.URLError, dwn._urlopen,"http://www.stolbeznog.tv")
        self.assertIsNotNone(dwn._urlopen, "https://www.google.com")


class TestPDBreading(TestCase):

    def setUp(self):
        self.pdbid = '1uak'
        self.chain = 'A'
        self.bad_chain = 'B'
        self.ligand = ['SAM']

    def test_access(self):
        self.assertTrue(dwn.check_if_in_PDB(self.pdbid,"A"))
        self.assertFalse(dwn.check_if_in_PDB(self.pdbid,self.bad_chain))
        self.assertEqual(dwn.get_ligands(self.pdbid,self.chain),self.ligand)
        self.assertEqual(len(list(dwn.get_particular_chain(self.pdbid,self.chain)[0])),243)
        self.assertEqual(list(dwn.get_particular_chain(self.pdbid,self.chain)[0]),dwn.get_all_chains(self.pdbid)[0])
        self.assertEqual(dwn.get_chain_list(self.pdbid),[self.chain])

