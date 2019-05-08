from unittest import TestCase

import knot_pull.config as cfg

class TestConfig(TestCase):
    def test_variables_for_production(self):
        self.assertFalse(cfg.VERBOSE)
        self.assertTrue(int(cfg.NUMBER_PRECISION_FUNCTION(5)) == 5)
        self.assertTrue(int(cfg.NUMBER_PRECISION_FUNCTION(2)/2) == 1)

