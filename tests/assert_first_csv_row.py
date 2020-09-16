#!/usr/local/bin/python

import pandas as pd
import argparse

def assert_first_csv_row(filename, field, value):
	df = pd.read_csv(filename)

	if field not in df:
		raise Exception('CSV does not contain field: %s (%s)'%(field, filename))

	if not len(df):
		raise Exception('No rows found in %s'%filename)

	if not df[field][0] == value:
		raise Exception("First entry in %s field mismatch (Got: %s \t Expected: %s)" % (field, df[field][0], value))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('filename', help='path to csv')
	parser.add_argument('field', help='name of column to check')
	parser.add_argument('value', help='what the column should read')

	args = parser.parse_args()

	assert_first_csv_row(args.filename, args.field, args.value)