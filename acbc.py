#!/usr/bin/env python
#
# Utilities/support routines for the ACBC pipeline code.
#
# David L. Swofford, 5 Feb 2019

#from __future__ import print_function	# for python 2 vs 3 compatibility

import os
import logging
import sys

def make_set_from_genelist_file(filename):
	"""
	Returns a set containing gene names contained in 'filename' (one gene per line).
	"""
	list = []
	if os.path.isfile(filename):
		list = [line.split()[0] for line in open(filename)]
		sys.stderr.write("make_set_from_genelist_file with filename={}: list --> {}\n".format(filename, list))
	else:
		sys.stderr.write("make_set_from_genelist_file: file {} doesn't exist!\n".format(filename))
	return set(list)

def make_genelist_file_from_set(gene_set, filename):
	"""
	Creates a file ('filename') containing elements of 'gene_set', one gene per line.
	"""
	with open(filename, 'w') as file:
		file.write("\n".join(sorted(gene_set)))

def delete_file(filename):
	"""
	Deletes file 'filename' (ignoring errors if it doesn't exist)/
	"""
	try:
		os.remove(filename)
	except OSError:
		pass

def py_which(cmd, mode=os.F_OK | os.X_OK, path=None):
	"""
	Given a command, mode, and a PATH string, return the path which conforms to the given mode
	on the PATH, or None if there is no such file.

	`mode` defaults to `os.F_OK | os.X_OK`. `path` defaults to the result of
	`os.environ.get("PATH")`, or can be overridden with a custom search path.
	"""
	# Check that a given file can be accessed with the correct mode.  Additionally check that
	# `file` is not a directory, as on Windows directories pass the os.access check.
	def _access_check(fn, mode):
		return (os.path.exists(fn) and os.access(fn, mode)
				and not os.path.isdir(fn))

	# If we're given a path with a directory part, look it up directly rather than referring to
	# PATH directories. This includes checking relative to the current directory, e.g. ./script
	if os.path.dirname(cmd):
		if _access_check(cmd, mode):
			return cmd
		return None

	if path is None:
		path = os.environ.get("PATH", os.defpath)
	if not path:
		return None
	path = path.split(os.pathsep)

	if sys.platform == "win32":
		# The current directory takes precedence on Windows.
		if not os.curdir in path:
			path.insert(0, os.curdir)

		# PATHEXT is necessary to check on Windows.  See if the given file matches any of the
		# expected path extensions.  This will allow us to short circuit when given "python.exe".
		# If it does match, only test that one, otherwise we have to try others.
		pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
		if any([cmd.lower().endswith(ext.lower()) for ext in pathext]):
			files = [cmd]
		else:
			files = [cmd + ext for ext in pathext]
	else:
		# On other platforms you don't have things like PATHEXT to tell you what file suffixes
		# are executable, so just pass on cmd as-is.
		files = [cmd]

	seen = set()
	for dir in path:
		normdir = os.path.normcase(dir)
		if not normdir in seen:
			seen.add(normdir)
			for thefile in files:
				name = os.path.join(dir, thefile)
				if _access_check(name, mode):
					return name

	return None

def set_logger(logger_name, debug=False):
	logger = logging.getLogger(logger_name)
	logger.setLevel(logging.DEBUG if debug else logging.INFO)
	consoleHandler = logging.StreamHandler()
	logger.addHandler(consoleHandler)
	format = logging.Formatter("#################### (%(name)s) %(levelname)s: %(message)s")
	consoleHandler.setFormatter(format)
	return logger
