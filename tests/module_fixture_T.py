# Creating a module fixture -- exercise from Python Testing book p. 113

from unittest import TestCase
from mocker import Mocker
from datetime import date

mocker = Mocker()

def setup():
    fake_date = mocker.replace(date)
    fake_date.today()
    mocker.result(date(year = 2009, month = 6, day = 12))
    mocker.count(1, None)
    mocker.replay()

def teardown():
    mocker.restore()
    mocker.verify()

class first_tests(TestCase):
    def test_year(self):
        self.assertEqual(date.today() .year, 2009)

    def test_month(self):
        self.assertEqual(date.today() .month, 6)

    def test_day(self):
        self.assertEqual(date.today() .day, 12)

class second_tests(TestCase):
    def test_isoformat(self):
        self.assertEqual(date.today() .isoformat(), '2009-06-12')
        