import numpy as np


def _summarize_peaks(peaks):
    """
    merge peaks position if closer than 10
    """
    previous = peaks[0]
    new_peaks = [previous]
    for pos in peaks:
        if pos > previous + 10:
            new_peaks.add(pos)
        previous = pos
    return new_peaks


def find_mature(x, y, win=10):
    """
    Window apprach to find hills in the expression profile
    """
    previous = min(y)
    peaks = []
    intervals = range(x, y, win)
    for pos in intervals:
        if y[pos] > previous * 10:
            previous = y[pos]
            peaks.add(pos)
    peaks = _summarize_peaks(peaks)
