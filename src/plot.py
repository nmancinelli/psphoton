#!/usr/bin/env python
#
def main():
	"""
	A simple program that plots phonon output.
	"""
	from numpy import zeros, log10
	infile=raw_input('Infile name: ')

	fin=open(infile,'r')
	hdr=fin.readline().strip('\n').split()
	x1=float(hdr[0])
	x2=float(hdr[1])
	nxdim=int(hdr[2])
	t1=float(hdr[3])
	t2=float(hdr[4])
	ntdim=int(hdr[5])
	
	rbin=zeros(ntdim*nxdim).reshape(ntdim,nxdim)

	for ix in range(nxdim):
		rbin[:,ix]=fin.readline().strip('\n').split()
		
	lbin=log10(rbin)
		
	import matplotlib
	matplotlib.use('TkAGG')
	from matplotlib import pylab as plt
	plt.imshow(lbin,extent=[x1,x2,t1,t2],aspect='auto',origin='lower')
	cbar=plt.colorbar()
	cbar.set_label('$\log_{10}$[ Amplitude ]')
	plt.xlabel('Distance (degrees)')
	plt.ylabel('Time (s)')
	plt.savefig('mypost.eps')
		

main()
