'''
FilePath: monitor.py
Author: wang yu
Date: 2024-11-03 11:01:43
LastEditTime: 2024-11-20 16:54:14
'''
"""
Name: Progress Monitor

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Function(s):
Get the progress and the time left
"""

from time import time


class Monitor:

    def __init__(self, t0):
        """
        introduction: Restore monitor settings.
        """
        self.position = 0
        self.last_time = t0

    # def restore(self):
    #     """
    #     introduction: Restore monitor settings.
    #     """
    #     self.__init__(self.last_time)

    def output(self, current_length, total_length, extra_informs=None):
        """
        introduction: Print the progress bar for required "for" sentence.

        :param current_length: Current position of the "for" sentence. (int)
        :param total_length: Total length of the "for" sentence. (int)
        :param extra_informs: extra information in the specific tasks. ()
        """
        if round(current_length / total_length * 100000) % 10 == 0:
            position = round(current_length / total_length * 100)

            self.position = position
            string = "["
            for index in range(100):
                if position > index:
                    string += "|"
                else:
                    string += " "
            string += "]"

            time_left = time() - self.last_time

            if extra_informs is None:
                string += "(detect = " + str(current_length) + ", total = " + str(total_length) + ") "
            else:
                string += "(detect = " + str(current_length) + ", {"
                for index, extra_inform in enumerate(extra_informs):
                    string += extra_inform[0] + " = " + extra_inform[1]
                    if index < len(extra_informs) - 1:
                        string += ", "
                string += "}, total = " + str(total_length) + ") "
            if self.position < 100:
                string += str(position) + "%, will be completed in " + str(
                    round(time_left,4)) + " seconds."
            else:
                string += str(position) + "%, was spent " + str(round(time_left, 4)) + " seconds."

            print(" " + string, end="\n")

            if current_length == total_length:
                print('\n')
