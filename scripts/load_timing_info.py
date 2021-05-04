"""Module for loading timing info gathered via TSCTimer class"""
import itertools
import pandas
import json
import statistics
import sys


def load_timings(filename):
    with open(filename) as f:
        timing_data = [json.loads(line) for line in f.readlines()]
    grouped_timing_data = list()
    for _, group in itertools.groupby(timing_data, lambda j: j['step']):
        grouped_timing_data.append(list(group))
    return grouped_timing_data


def get_durations(timing_data_group):
    return {t['name']: t['duration'] for t in timing_data_group}


def get_duration_stats(timing_data):
    durations = list(map(get_durations, timing_data))
    duration_moments = dict()
    for key in durations[0]:
        dur_vals = [duration[key] for duration in durations]
        duration_moments[key] = {
            'mean': statistics.mean(dur_vals),
            'median': statistics.median(dur_vals),
            'std': statistics.pstdev(dur_vals)
        }
    return duration_moments


if __name__ == '__main__':
    filename = sys.argv[1] if len(sys.argv) > 1 else 'timings.json'
    print(f"Loading timing info from {filename}...")
    t = load_timings(filename)
    stats = get_duration_stats(t)
    print("Timing information:")
    print(f"{'Section' : >25}\t Mean/Median/Std")
    print(60*"=")
    for section, s in stats.items():
        print(f"{section : >25}\t", '/'.join(map("{:g}".format, (s['mean'], s['median'], s['std']))))
