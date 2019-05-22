from unittest import TestCase
import os
path = os.path.dirname(os.path.abspath(__file__))

import knot_pull.reader as rdr

class TestOnline(TestCase):
    def setUp(self):
        self.pdbid_onechain = '1uak'
        self.chain_onechain = 'A'
        self.pdbid_multi = '4mcb'

    def test_err_works(self):
        self.assertTrue(rdr.read_from_web(self.pdbid_onechain,'C')[2])
        self.assertFalse(rdr.read_from_web(self.pdbid_onechain,'C')[0])
        self.assertFalse(rdr.read_from_web(self.pdbid_onechain)[2])

    def test_be_works(self):
        self.assertEqual(str(rdr.read_from_web(self.pdbid_multi,'A')),
                         str(rdr.read_from_web(self.pdbid_multi,begin='1',end='219')))
        self.assertEqual(str(rdr.read_from_web(self.pdbid_onechain)),
                         str(rdr.read_from_web(self.pdbid_onechain, begin='i-1', end='i250')))


class TestPDB(TestCase):
    def setUp(self):
        self.pdb_file = '{}/2efv.pdb'.format(path)

    def test_be_works(self):
        self.assertTrue(rdr.read_from_pdb(self.pdb_file,"Z")[2])
        self.assertEqual(str(rdr.read_from_pdb(self.pdb_file)), str(rdr.read_from_pdb(self.pdb_file,"A")))
        self.assertTrue(rdr.read_from_pdb(self.pdb_file, begin='50',end='i12')[2])
        self.assertEqual(str(rdr.read_from_pdb(self.pdb_file)),str(rdr.read_from_pdb(self.pdb_file,begin='i6',end='82')))
        self.assertEqual(str(rdr.read_from_pdb(self.pdb_file)),str(rdr.read_from_pdb(self.pdb_file,begin='1',end='i87')))

class TestXYZ(TestCase):
    def setUp(self):
        self.xyz_file = '{}/05_123.xyz'.format(path)

    def test_be_works(self):
        self.assertTrue(rdr.read_from_xyz(self.xyz_file, begin='13',end='107')[2])
        self.assertNotEqual(str(rdr.read_from_xyz(self.xyz_file)),
                            str(rdr.read_from_xyz(self.xyz_file, begin='1', end='i22')))
        self.assertEqual(str(rdr.read_from_xyz(self.xyz_file,begin='1',end='i22')),
                         str(rdr.read_from_xyz(self.xyz_file,begin='i11',end='12')))


