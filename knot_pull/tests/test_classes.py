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


class TestLine(TestCase):

    def setUp(self):
        self.line1 = kpc.Line(5,1)
        self.line2 = kpc.Line(6,1)

    def test_line_behaviour(self):
        self.assertEqual(self.line1.topsign(),"")
        self.assertEqual(self.line2.topsign(),"-")
        self.assertTrue(self.line1<self.line2)
        self.assertTrue(self.line2>10)
        self.assertTrue(self.line2+10 == 16)
        self.assertTrue(self.line2-2 == 4)
        self.assertTrue(int(self.line2) == 6)
        self.assertTrue(self.line1%2)


class TestCrossing(TestCase):

    def setUp(self):
        self.line1 = kpc.Line(5,1)
        self.line2 = kpc.Line(4,0)
        self.line3 = kpc.Line(6,1)
        self.line4 = kpc.Line(9,0)
        self.cr1 = kpc.Crossing(5,1,4,0)
        self.cr2 = kpc.Crossing(9,0,6,1)

    def test_crossing_funcs(self):
        self.assertEqual(self.cr1[0],self.line1)
        self.assertEqual(self.cr1.l1.top,self.line1.top)

        self.cr1.reverse_topo()
        self.assertTrue(self.cr1[1].top == 1)

        self.assertEqual(self.cr2.min(),6)
        self.assertEqual(self.cr2.max(),9)
        self.assertEqual(self.cr2.uneven(),9)

        self.cr2.remove_previous(2)
        self.assertEqual(self.cr2.min(), 5)
        self.assertEqual(self.cr2.max(), 8)
        self.assertEqual(self.cr2.even(),8)
        self.assertEqual(self.cr2.values(), [8,5])
        self.assertEqual(self.cr2.to_code(), (5,8))
        self.assertEqual(self.cr2.to_mod_code(), (5,8))
        self.assertEqual(self.cr2.len_loop(), 3)
        self.assertEqual(self.cr2.sum_loop(), 13)
        self.assertTrue(self.cr2.has_values([5,8]))
        self.assertFalse(self.cr2.has_values([6,9]))


class TestCode(TestCase):

    def setUp(self):
        _code = [(1,-4),(3,-6),(5,-2),(6,-7)]
        self.code1 = kpc.Code()
        self.code1.read_in(_code)

        self.code2 = kpc.Code()
        for _c1,_c2 in _code:
            self.code2.add(kpc.Crossing(abs(_c1), _c1 < 0, abs(_c2), _c2 < 0))

        self.code_bad1 = kpc.Code()
        self.code_bad1.read_in([(1,-4),(3,-5)])

        self.code_bad2 = kpc.Code()
        self.code_bad2.read_in([(1, 4), (3, -5)])

    def test_code_behavior(self):
        self.assertEqual(self.code1[2], self.code2[2])
        self.assertIsNone(self.code1.check_yo())
        self.assertRaises(ValueError,self.code_bad1.check_yo)
        self.assertRaises(ValueError,self.code_bad2.check_yo)

        c = self.code1[0]
        self.code1.start_later_by(3)
        self.assertIs(c,self.code1[0])

        self.assertEqual(2,self.code2.index(self.code2[2]))

        self.code_bad1.make_way(5)
        self.assertIsNone(self.code1.check_yo())