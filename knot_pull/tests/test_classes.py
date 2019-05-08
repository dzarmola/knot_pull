from unittest import TestCase
from random import random

import knot_pull.kpclasses as kpc


class TestBead(TestCase):

    def setUp(self):
        self.bead = kpc.Bead([1,2,3])
        self.beadx = kpc.Bead([1,1,3],"X")
        self.bead.setNhand(self.beadx)
        self.beadx.setNhand(None)

    def test_bead_creation(self):
        self.assertTrue(self.bead.y == 2)
        self.assertTrue(str(self.bead) == "1,2,3")
        self.assertTrue(self.bead.isCa())
        self.assertFalse(self.beadx.isCa())

    def test_bead_connection(self):
        self.assertTrue(self.bead.Ndist == 1.)
        self.assertIsNone(self.beadx.Ndist)


class TestChainCopy(TestCase):

    def setUp(self):
        self.atom_chain = [kpc.Bead([random() for _ in range(3)]) for _ in range(10)]

    def test_chain_copy(self):
        _copy = kpc.chainDeepCopy(self.atom_chain)
        self.assertTrue(_copy[3].Nhand.Chand.y == _copy[3].y)
        self.assertTrue(_copy[2].Chand.y == _copy[3].y)
        self.assertTrue(self.atom_chain[3].y == _copy[3].y)
        self.atom_chain[3].y /= 2.
        self.assertFalse(self.atom_chain[3].y == _copy[3].y)
