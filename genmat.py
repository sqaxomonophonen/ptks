#!/usr/bin/env python

import sys
import struct
import Image

def open_texture(d):
	diffuse_file, emission_file = d.split(":")
	diffuse = None
	emission = None
	if len(diffuse_file) > 0:
		diffuse = Image.open(diffuse_file).convert("RGB")

	if len(emission_file) > 0:
		emission = Image.open(emission_file).convert("RGB")
	return diffuse, emission

for d in sys.argv[1:]:
	diffuse, emission = open_texture(d)

	width = 128
	height = 128

	for y in range(height):
		for x in range(width):
			d = [0,0,0]
			if diffuse:
				d = list(diffuse.getpixel((x,y)))
			e = [0,0,0]
			if emission:
				e = list(emission.getpixel((x,y)))
			sys.stdout.write(struct.pack("BBBBBB", *(d + e)))

